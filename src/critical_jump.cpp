// Continuous-time Gillespie simulation of bosonic number and phase transitions.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <regex>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
// Configuration schema matches qtp.io. Metadata describes recorded runs and
// does not change rates. Unknown keys are rejected before defaults are applied.
const std::set<std::string> integer_keys = {
    "nx", "ny", "mode", "phaseinit", "phase_draw_order", "seed", "disseed"
};
const std::set<std::string> float_keys = {
    "eps", "U", "delta", "D2", "D4", "zeta", "nu", "b", "c", "r", "eta",
    "leak", "repair", "Db", "smooth", "Dcl", "disorder", "T", "save", "m0",
    "t1", "t2", "eta_mid"
};
const std::set<std::string> metadata_keys = {
    "family", "group", "id", "parent", "block", "rep", "stream", "tf", "E0",
    "initseed", "sigma_init", "engine_sha256", "initial_sha256",
    "configuration_source_sha256", "configuration_source_format", "source_id"
};

// Parse an exact, nonnegative decimal integer without passing through double.
// Integral decimal/scientific forms (e.g. 5.0, 5e0) are accepted. Fractional
// values, overflow and negative nonzero values are rejected, including seeds.
uint64_t parse_integer(const std::string& text, const std::string& key,
                       uint64_t limit) {
    static const std::regex syntax(
        R"([+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?)");
    auto invalid = [&]() -> std::runtime_error {
        return std::runtime_error("invalid integer parameter: " + key + "=" + text);
    };
    if (!std::regex_match(text, syntax)) throw invalid();
    bool negative = text[0] == '-';
    size_t begin = (text[0] == '-' || text[0] == '+') ? 1 : 0;
    size_t exp_pos = text.find_first_of("eE", begin);
    std::string mantissa = text.substr(begin, exp_pos - begin);
    size_t dot = mantissa.find('.');
    size_t decimals = dot == std::string::npos ? 0 : mantissa.size() - dot - 1;
    std::string digits = mantissa;
    if (dot != std::string::npos) digits.erase(dot, 1);
    size_t first = digits.find_first_not_of('0');
    if (first == std::string::npos) return 0;
    if (negative) throw invalid();
    digits.erase(0, first);
    long long exponent = 0;
    if (exp_pos != std::string::npos) {
        try {
            exponent = std::stoll(text.substr(exp_pos + 1));
        } catch (const std::exception&) {
            throw invalid();
        }
    }
    if (exponent > 100000 || exponent < -100000 || decimals > 100000) {
        throw invalid();
    }
    long long scale = exponent - static_cast<long long>(decimals);
    if (scale < 0) {
        size_t trim = static_cast<size_t>(-scale);
        if (trim >= digits.size()) throw invalid();
        for (size_t i = digits.size() - trim; i < digits.size(); ++i) {
            if (digits[i] != '0') throw invalid();
        }
        digits.resize(digits.size() - trim);
    } else {
        if (digits.size() + static_cast<size_t>(scale) > 20) throw invalid();
        digits.append(static_cast<size_t>(scale), '0');
    }
    uint64_t value = 0;
    for (char ch : digits) {
        uint64_t digit = static_cast<uint64_t>(ch - '0');
        if (value > (limit - digit) / 10) throw invalid();
        value = 10 * value + digit;
    }
    return value;
}

double parse_number(const std::string& text, const std::string& key) {
    try {
        size_t used = 0;
        double value = std::stod(text, &used);
        if (used == text.size() && std::isfinite(value)) return value;
    } catch (const std::exception&) {
        // All conversion errors use the same parameter-specific diagnostic.
    }
    throw std::runtime_error("invalid numeric parameter: " + key + "=" + text);
}

// Fenwick tree for cumulative event rates and O(log M) event selection.
struct Tree{
    int n;
    std::vector<double> a,v;
    Tree(int N):n(N),a(N+1,0),v(N,0){
    }
    void set(int i,double x){
        if(!(x>=0)||!std::isfinite(x))throw std::runtime_error("invalid rate");
        double d=x-v[i];
        v[i]=x;
        for(int j=i+1;j<=n;j+=j&-j)a[j]+=d;
    }
    double total()const{
        double x=0;
        for(int j=n;j;j-=j&-j)x+=a[j];
        return x;
    }
    int choose(double x)const{
        int i=0,bit=1;
        while(2*bit<=n)bit*=2;
        for(;bit;bit/=2){
            int j=i+bit;
            if(j<=n&&a[j]<=x){
                i=j;
                x-=a[j];
            }
        }
        return std::min(i,n-1);
    }
};
int main(int argc,char**argv){
    try{
        if(argc!=3)throw std::runtime_error("usage: critical_jump CONFIG OUTPUT_PREFIX");
        const std::filesystem::path cfg_path=std::filesystem::absolute(argv[1]);
        std::ifstream cfg(cfg_path);
        if(!cfg)throw std::runtime_error("cannot open configuration: "+cfg_path.string());
        std::map<std::string,std::string> par;
        std::string line;
        size_t line_number = 0;
        while(std::getline(cfg,line)){
            ++line_number;
            std::istringstream row(line);
            std::string key,val;
            if(!(row>>key)||key[0]=='#')continue;
            if (!integer_keys.count(key) && !float_keys.count(key) &&
                !metadata_keys.count(key) && key != "init") {
                throw std::runtime_error(cfg_path.string() + ":" +
                    std::to_string(line_number) + ": unknown configuration key: " + key);
            }
            std::getline(row>>std::ws,val);
            const auto end=val.find_last_not_of(" \t\r");
            if(end==std::string::npos)throw std::runtime_error("missing value for "+key);
            val.erase(end+1);
            if(par.count(key))throw std::runtime_error("duplicate configuration key: "+key);
            par[key]=val;
        }
        for (const auto& key : {"seed", "disseed", "init"}) {
            if (!par.count(key)) {
                throw std::runtime_error(std::string("missing required parameter: ") + key);
            }
        }
        std::map<std::string, uint64_t> integers;
        std::map<std::string, double> numbers;
        for (const auto& entry : par) {
            const auto& key = entry.first;
            if (integer_keys.count(key)) {
                const uint64_t limit = (key == "seed" || key == "disseed")
                    ? std::numeric_limits<uint64_t>::max()
                    : static_cast<uint64_t>(std::numeric_limits<int>::max());
                integers[key] = parse_integer(entry.second, key, limit);
            } else if (float_keys.count(key)) {
                numbers[key] = parse_number(entry.second, key);
            }
        }
        auto get = [&](const std::string& key, double fallback) {
            if (integers.count(key)) return static_cast<double>(integers.at(key));
            if (numbers.count(key)) return numbers.at(key);
            return fallback;
        };
        if(get("nx",24)<5||get("ny",28)<5||get("nx",24)*get("ny",28)>std::numeric_limits<int>::max()/2) throw std::runtime_error("require nx, ny >= 5 and a supported array size");
        int nx=get("nx",24),ny=get("ny",28),M=nx*ny,mode=get("mode",0);
        if(nx<5||ny<5)throw std::runtime_error("use nx,ny >=5");
        double eps=get("eps",.0005),U=get("U",1),delta=get("delta",.3),
            D2=get("D2",1.6),D4=get("D4",1),zeta=get("zeta",0),nu=get("nu",0),
            b=get("b",1),c=get("c",.1),r=get("r",20),eta=get("eta",1),leak=get("leak",.0002),
            repair=get("repair",2),Db=get("Db",.00001),smooth=get("smooth",.002),
            Dcl=get("Dcl",.001),dis=get("disorder",.05),T=get("T",250),
            save=get("save",2);
        if(eps<=0||U<=0||save<=0||T<0||eta<0||Db<0||leak<0||repair<0||r<0||D2<0||D4<=0||b<=0||c<0||zeta<0||smooth<=0||Dcl<0||delta<=0||dis<0||dis>=1||(mode!=0&&mode!=1)) throw std::runtime_error("invalid physical rates, time, disorder, or mode");
        if(U/eps>std::numeric_limits<int>::max()/2||std::abs(U/eps-std::round(U/eps))>1e-8) throw std::runtime_error("U/eps must be a supported positive integer");
        double t1=get("t1",T+1),t2=get("t2",T+2),etamid=get("eta_mid",eta),eta0=eta;
        int stage=0;
        if(!std::isfinite(t1)||!std::isfinite(t2)||t1<0||t2<t1||etamid<0)
            throw std::runtime_error("invalid dephasing schedule");
        uint64_t seed=integers.at("seed"),disseed=integers.at("disseed");
        int ncap=std::lround(U/eps);
        std::vector<int32_t> n(M);
        std::vector<int8_t>s(M);
        std::vector<double> u(M),base(M),up(M),down(M),dd(M),bx(M),by(M);
        std::mt19937_64 rng(seed),drng(disseed);
        std::uniform_real_distribution<double>uni(0,1);
        auto randu=[&](){
            return std::max(uni(rng),1e-300);
        };
        double m0=get("m0",0);
        if(m0 < -1 || m0 > 1) throw std::runtime_error("invalid m0");
        // Select the uniform-to-phase initialization map for this dataset.
        // Both conventions have the same phase probabilities and transition rates.
        int phase_draw_order=get("phase_draw_order",0);
        if(phase_draw_order!=0&&phase_draw_order!=1)throw std::runtime_error("phase_draw_order must be 0 or 1");
        int phaseinit=get("phaseinit",0);
        if(phaseinit!=0&&phaseinit!=1)throw std::runtime_error("phaseinit must be 0 or 1");
        std::filesystem::path init_path(par.at("init"));
        if(init_path.is_relative())init_path=cfg_path.parent_path()/init_path;
        std::ifstream ini(init_path);
        if(!ini)throw std::runtime_error("cannot open initial configuration: "+init_path.string());
        for(int i=0;i<M;i++){
            int ns,ss;
            if(!(ini>>ns>>ss))throw std::runtime_error("invalid or truncated initial configuration: "+init_path.string());
            if(ns < -ncap || ns > ncap || (ss!=-1&&ss!=0&&ss!=1&&ss!=9)) throw std::runtime_error("initial occupation or label outside allowed range");
            n[i]=ns;
            s[i]=ss==9?(phase_draw_order==0?(randu()<(1+m0)/2?1:-1):(randu()<(1-m0)/2?-1:1)):ss;
            if(phaseinit==1)s[i]=1;
            if(mode==1&&s[i]!=0)s[i]=1;
            u[i]=eps*(n[i]-(s[i]==0?.5:0));
            dd[i]=delta*(1+dis*(2*uni(drng)-1));
            bx[i]=D2*(1+dis*(2*uni(drng)-1));
            by[i]=D2*(1+dis*(2*uni(drng)-1));
        }
        std::string extra_initial_value;
        if (ini >> extra_initial_value) {
            throw std::runtime_error("extra entries in initial configuration: " + init_path.string());
        }
        auto index=[&](int x,int y){
            return ((y%ny+ny)%ny)*nx+(x%nx+nx)%nx;
        };
        std::vector<std::vector<int>>nb(M),supp(M);
        std::vector<std::vector<double>>co(M),edge(M);
        int dx[13]={
            0,1,-1,0,0,1,1,-1,-1,2,-2,0,0
        },dy[13]={
            0,0,0,1,-1,1,-1,1,-1,0,0,2,-2
        };
        for(int i=0;i<M;i++){
            int x=i%nx,y=i/nx;
            nb[i]={
                index(x+1,y),index(x-1,y),index(x,y+1),index(x,y-1)
            };
            edge[i]={
                bx[i],bx[index(x-1,y)],by[i],by[index(x,y-1)]
            };
            double cf[13]={
                dd[i]+20*D4+2*zeta-4*Dcl,-8*D4+Dcl,-8*D4+Dcl,-8*D4-zeta+Dcl,-8*D4-zeta+Dcl,2*D4,2*D4,2*D4,2*D4,D4,D4,D4,D4
            };
            for(int a=0;a<13;a++){
                int j=index(x+dx[a],y+dy[a]);
                auto it=std::find(supp[i].begin(),supp[i].end(),j);
                if(it==supp[i].end()){
                    supp[i].push_back(j);
                    co[i].push_back(cf[a]);
                }
                else{
                    co[i][it-supp[i].begin()]+=cf[a];
                }
            }
            for(size_t a=0;a<supp[i].size();a++)base[i]+=co[i][a]*u[supp[i][a]];
        }
        Tree tree(2*M);
        double background=Db/(eps*eps);
        // Field jump rates, Eq. (9); mode 1 uses the exact sign average in Appendix C.
        auto rate_update=[&](int i){
            double v=u[i],f=-base[i]+nu*v*v-b*v*v*v-c*v*v*v*v*v,ru=0,rd=0;
            if(mode==0){
                for(int a=0;a<4;a++)f+=edge[i][a]*.5*(int(s[i])+int(s[nb[i][a]]))*(v-u[nb[i][a]]);
                double q=std::hypot(f,smooth);
                ru=(q+f)/(2*eps)+background;
                rd=(q-f)/(2*eps)+background;
            }
            else {
                for(int mask=0;mask<32;mask++){
                    double fz=f;
                    int si=s[i]==0?0:((mask&1)?1:-1);
                    for(int a=0;a<4;a++){
                        int sj=s[nb[i][a]]==0?0:((mask>>(a+1)&1)?1:-1);
                        fz+=edge[i][a]*.5*(si+sj)*(v-u[nb[i][a]]);
                    }
                    double q=std::hypot(fz,smooth);
                    ru+=(q+fz)/(64*eps);
                    rd+=(q-fz)/(64*eps);
                }
                ru+=background;
                rd+=background;
            }
            if(n[i]>=ncap)ru=0;
            if(n[i]<=-ncap)rd=0;
            up[i]=ru;
            down[i]=rd;
            tree.set(i,ru+rd);
        };
        // Local coherence source, dephasing, leakage and repair rates.
        auto phase_update=[&](int i){
            double z=s[i]==0?repair:leak+(mode==1?0:(s[i]<0?r+eta:eta));
            tree.set(M+i,z);
        };
        for(int i=0;i<M;i++){
            rate_update(i);
            phase_update(i);
        }
        std::string prefix=argv[2];
        auto output_parent=std::filesystem::path(prefix).parent_path();
        if(!output_parent.empty())std::filesystem::create_directories(output_parent);
        for(const auto& ext:{
            ".bin",".csv","_engine.json"
        }){
            if(std::filesystem::exists(prefix+ext))throw std::runtime_error("output already exists: "+prefix+ext);
        }
        std::ofstream out(prefix+".bin",std::ios::binary),log(prefix+".csv");
        if(!out||!log)throw std::runtime_error("cannot create outputs for prefix: "+prefix);
        int32_t hd[3]={
            nx,ny,M
        };
        out.write((char*)hd,sizeof hd);
        log<<std::setprecision(17)<<"time,energy,mean,coherence,leak_fraction,transport_power,weighted_leak,field_events,phase_events,leak_events,repair_events,max_abs\n";
        uint64_t nf=0,np=0,nl=0,nr=0;
        double time=0,next=0,maxerr=0;
        uint64_t boundary=0;
        // Save centered occupation labels and number/phase observables.
        auto write=[&](double t){
            double E=0,av=0,cs=0,ls=0,pow=0,wl=0,mx=0;
            for(int i=0;i<M;i++){
                E+=u[i]*u[i];
                av+=u[i];
                cs+=mode==0?s[i]:0;
                ls+=s[i]==0;
                mx=std::max(mx,std::abs(u[i]));
                for(int a=0;a<4;a++){
                    double v=u[i]-u[nb[i][a]];
                    pow+=edge[i][a]*(mode==0?s[i]:0)*v*v;
                    wl+=(s[i]==0)*v*v;
                }
            }
            out.write((char*)&t,8);
            out.write((char*)n.data(),M*4);
            out.write((char*)s.data(),M);
            log<<t<<","<<E/M<<","<<av/M<<","<<cs/M<<","<<ls/M<<","<<pow/M<<","<<wl/M<<","<<nf<<","<<np<<","<<nl<<","<<nr<<","<<mx<<"\n";
        };
        // Optional externally controlled dephasing window; the field generator is unchanged.

        write(0);
        next=save;
        while(time<T){
            double rate=tree.total();
            if(!(rate>0))throw std::runtime_error("unexpected absorbing state");
            double event=time-std::log(randu())/rate;
            double sw=stage==0?t1:(stage==1?t2:T+10);
            double stop=std::min(T,std::min(event,sw));
            while(next<=stop+1e-10){
                write(std::min(next,T));
                next+=save;
            }
            if(sw<=event&&sw<T){
                time=sw;
                stage++;
                eta=stage==1?etamid:eta0;
                for(int i=0;i<M;i++)phase_update(i);
                continue;
            }
            if(event>T)break;
            time=event;
            int j=tree.choose(randu()*rate);
            double du=0;
            bool changed_u=false;
            if(j<M){
                int step=randu()*(up[j]+down[j])<up[j]?1:-1;
                n[j]+=step;
                du=eps*step;
                u[j]=eps*(n[j]-(s[j]==0?.5:0));
                nf++;
                changed_u=true;
                if(std::abs(n[j])>=ncap)boundary++;
            }
            else{
                j-=M;
                int old=s[j];
                if(old==0){
                    s[j]=mode==1?1:(randu()<.5?-1:1);
                    du=eps/2;
                    nr++;
                    changed_u=true;
                }
                else{
                    double flip=mode==1?0:(old<0?r+eta:eta);
                    if(randu()*(flip+leak)<leak){
                        s[j]=0;
                        du=-eps/2;
                        nl++;
                        changed_u=true;
                    }
                    else{
                        s[j]=-s[j];
                        np++;
                    }
                }
                u[j]=eps*(n[j]-(s[j]==0?.5:0));
                phase_update(j);
            }
            if(changed_u){
                for(size_t a=0;a<supp[j].size();a++){
                    int i=supp[j][a];
                    base[i]+=co[j][a]*du;
                }
                for(int i:supp[j])rate_update(i);
            }
            else{
                rate_update(j);
                for(int i:nb[j])rate_update(i);
            }
            if((nf+np+nl+nr)%1000000==0){
                for(int i=0;i<M;i++){
                    double z=0;
                    for(size_t a=0;a<supp[i].size();a++)z+=co[i][a]*u[supp[i][a]];
                    maxerr=std::max(maxerr,std::abs(z-base[i]));
                    base[i]=z;
                    rate_update(i);
                }
            }
        }
        if(next-save<T-1e-9)write(T);
        if(!out||!log)throw std::runtime_error("failed while writing trajectory: "+prefix);
        std::ofstream meta(prefix+"_engine.json");
        if(!meta)throw std::runtime_error("cannot write event diagnostics: "+prefix);
        meta<<std::setprecision(17)<<"{\"field_events\":"<<nf<<",\"phase_events\":"<<np<<",\"leak_events\":"<<nl<<",\"repair_events\":"<<nr<<",\"boundary_hits\":"<<boundary<<",\"stencil_error\":"<<maxerr<<"}\n";
        return 0;
    }
    catch(const std::exception&e){
        std::cerr<<"critical_jump: "<<e.what()<<"\n";
        return 1;
    }
}
