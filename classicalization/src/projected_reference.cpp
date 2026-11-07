// Independent exact physical-number-projected control of the autonomous model.
// Nonlinear rates are averaged over inactive relative-phase signs before sampling.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Fenwick {
    int n; std::vector<double> bit, val;
    explicit Fenwick(int n): n(n), bit(n+1,0), val(n,0) {}
    void set(int i,double v) {
        if(v<0 || !std::isfinite(v)) throw std::runtime_error("Invalid event rate");
        const double d=v-val[i]; val[i]=v;
        for(int j=i+1;j<=n;j+=j&-j) bit[j]+=d;
    }
    double total() const { double a=0;for(int j=n;j;j-=j&-j)a+=bit[j];return a; }
    int choose(double v) const {
        int i=0,k=1;while(2*k<=n)k*=2;
        for(;k;k/=2) {int j=i+k;if(j<=n && bit[j]<=v) {i=j;v-=bit[j];}}
        return std::min(i,n-1);
    }
};

int main(int argc,char**argv) {
 try {
    if(argc!=3) throw std::runtime_error("Usage: autonomous_jump CONFIG OUTPUT_PREFIX");
    const auto cfg=std::filesystem::absolute(argv[1]);
    std::ifstream input(cfg); if(!input)throw std::runtime_error("Cannot read configuration");
    const std::set<std::string> allowed={"nx","ny","eps","U","delta","D2","D4","zeta","nu","b","c","r","eta","leak","repair","Db","smooth","Dcl","disorder","T","save","seed","disseed","init","gamma","beta","activity_scale","classical_efficiency","quantum_enabled","environment_initial","m0","verify_every"};
    std::map<std::string,std::string> p;std::string line;
    while(std::getline(input,line)) {
        std::istringstream row(line);std::string key,value;
        if(!(row>>key)||key[0]=='#')continue;
        std::getline(row>>std::ws,value);
        if(!allowed.count(key)||value.empty()||p.count(key))throw std::runtime_error("Invalid/duplicate key: "+key);
        p[key]=value;
    }
    auto get=[&](const std::string& k,double def){if(!p.count(k))return def;size_t used;double v=std::stod(p.at(k),&used);if(used!=p.at(k).size()||!std::isfinite(v))throw std::runtime_error("Invalid value: "+k);return v;};
    auto integer=[&](const std::string& k,int def){double v=get(k,def);if(v!=std::round(v)||v<0||v>1000000000)throw std::runtime_error("Invalid integer: "+k);return int(v);};
    auto seed_value=[&](const std::string& k){if(!p.count(k))throw std::runtime_error("Missing "+k);size_t j;auto v=std::stoull(p.at(k),&j);if(j!=p.at(k).size()||p.at(k)[0]=='-')throw std::runtime_error("Invalid seed");return v;};
    const int nx=integer("nx",48),ny=integer("ny",56),M=nx*ny;
    const bool quantum=false; // Exact number projection is imposed in this control only.
    const int initial_h=integer("environment_initial",0);
    const double eps=get("eps",.01),U=get("U",1),delta=get("delta",.3),D2=get("D2",1.6),D4=get("D4",1),zeta=get("zeta",0),nu=get("nu",0),b=get("b",1),c=get("c",.1),r=get("r",20),eta=get("eta",1),leak=get("leak",.0002),repair=get("repair",2),Db=get("Db",.00001),smooth=get("smooth",.002),Dcl=get("Dcl",.001),dis=get("disorder",.05),T=get("T",1200),save=get("save",5),gamma=get("gamma",.002),beta=get("beta",.08),scale=get("activity_scale",.10),eff=get("classical_efficiency",.91),m0=get("m0",0);
    if(nx<5||ny<5||M>1000000||eps<=0||U<=0||D4<=0||b<=0||delta<=0||smooth<=0||scale<=0||T<0||save<=0||dis<0||dis>=1||initial_h>1||std::abs(m0)>1)throw std::runtime_error("Invalid geometry or positive parameter");
    for(double x:{D2,zeta,c,r,eta,leak,repair,Db,Dcl,gamma,beta,eff})if(x<0)throw std::runtime_error("Negative rate");
    if(std::abs(U/eps-std::round(U/eps))>1e-8)throw std::runtime_error("U/eps must be integral");
    const int cap=std::lround(U/eps),verify_every=integer("verify_every",1000000);
    std::mt19937_64 rng(seed_value("seed")),drng(seed_value("disseed"));
    std::uniform_real_distribution<double> uniform(0,1);
    auto random=[&](){return std::max(uniform(rng),1e-300);};
    std::vector<int32_t> n(M); std::vector<int8_t> s(M),h(M,initial_h);
    std::vector<double> u(M),base(M),up(M),down(M),env_local(M),env_spread(M),dd(M),bx(M),by(M);
    std::vector<std::array<int,4>> nb(M);
    std::vector<std::array<int,13>> support(M);
    std::vector<std::array<double,13>> coeff(M);
    std::vector<std::array<double,4>> edge(M);
    if(!p.count("init"))throw std::runtime_error("Missing init");
    std::filesystem::path ip=p.at("init");if(ip.is_relative())ip=cfg.parent_path()/ip;
    std::ifstream ini(ip);if(!ini)throw std::runtime_error("Cannot read initial state");
    for(int i=0;i<M;i++) {
        int ns,ss;if(!(ini>>ns>>ss)||ns< -cap||ns>cap||(ss!=-1&&ss!=0&&ss!=1&&ss!=9))throw std::runtime_error("Invalid initial state");
        n[i]=ns;s[i]=ss==9?(random()<(1+m0)/2?1:-1):ss;
        if((!quantum||h[i])&&s[i]!=0)s[i]=1;
        u[i]=eps*(n[i]-(s[i]==0?.5:0));
        dd[i]=delta*(1+dis*(2*uniform(drng)-1));
        bx[i]=D2*(1+dis*(2*uniform(drng)-1));by[i]=D2*(1+dis*(2*uniform(drng)-1));
    }
    std::string extra;if(ini>>extra)throw std::runtime_error("Extra initialization entries");
    auto index=[&](int x,int y){return ((y%ny+ny)%ny)*nx+(x%nx+nx)%nx;};
    const int dx[13]={0,1,-1,0,0,1,1,-1,-1,2,-2,0,0};
    const int dy[13]={0,0,0,1,-1,1,-1,1,-1,0,0,2,-2};
    for(int i=0;i<M;i++) {
        int x=i%nx,y=i/nx;nb[i]={index(x+1,y),index(x-1,y),index(x,y+1),index(x,y-1)};
        edge[i]={bx[i],bx[index(x-1,y)],by[i],by[index(x,y-1)]};
        coeff[i]={dd[i]+20*D4+2*zeta-4*Dcl,-8*D4+Dcl,-8*D4+Dcl,-8*D4-zeta+Dcl,-8*D4-zeta+Dcl,2*D4,2*D4,2*D4,2*D4,D4,D4,D4,D4};
        for(int k=0;k<13;k++){support[i][k]=index(x+dx[k],y+dy[k]);base[i]+=coeff[i][k]*u[support[i][k]];}
    }
    Fenwick tree(3*M);
    auto update_field=[&](int i){
        const double v=u[i];
        double f=-base[i]+nu*v*v-b*v*v*v-c*v*v*v*v*v;
        double ai=0; std::array<double,5> a{}; int count=0;
        for(int k=0;k<4;k++) {
            const int j=nb[i][k];const double w=edge[i][k]*.5*(v-u[j]);
            f+=w*eff*(h[i]+h[j]);ai+=w;
            if(!h[j]&&s[j]!=0)a[count++]=w;
        }
        if(!h[i]&&s[i]!=0)a[count++]=ai;
        double mean_root=0;
        for(int mask=0;mask<(1<<count);mask++) {
            double ff=f;for(int k=0;k<count;k++)ff+=((mask>>k)&1)?a[k]:-a[k];
            mean_root+=std::hypot(ff,smooth);
        }
        mean_root/=double(1<<count);
        up[i]=n[i]>=cap?0:(mean_root+f)/(2*eps)+Db/(eps*eps);
        down[i]=n[i]<=-cap?0:(mean_root-f)/(2*eps)+Db/(eps*eps);
        tree.set(i,up[i]+down[i]);
    };
    auto update_phase=[&](int i){double flip=quantum&&!h[i]?(s[i]<0?r+eta:eta):0;tree.set(M+i,s[i]==0?repair:leak+flip);};
    auto update_environment=[&](int i){
        if(h[i]){env_local[i]=env_spread[i]=0;tree.set(2*M+i,0);return;}
        double a=u[i]*u[i],activated=0;
        for(int k:nb[i]){a+=u[k]*u[k];activated+=h[k];}
        const double z=a/(5*scale),z2=z*z,z4=z2*z2;
        env_local[i]=gamma*z4/(1+z4);env_spread[i]=beta*activated/4;
        tree.set(2*M+i,env_local[i]+env_spread[i]);
    };
    for(int i=0;i<M;i++){update_field(i);update_phase(i);update_environment(i);}
    const std::string prefix=argv[2];auto parent=std::filesystem::path(prefix).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
    for(auto ext:{".bin",".csv","_events.csv","_engine.json"})if(std::filesystem::exists(prefix+ext))throw std::runtime_error("Output exists: "+prefix+ext);
    std::ofstream out(prefix+".bin",std::ios::binary),log(prefix+".csv"),events(prefix+"_events.csv");
    if(!out||!log||!events)throw std::runtime_error("Cannot write outputs");
    int32_t hdr[3]={nx,ny,M};out.write((char*)hdr,sizeof hdr);
    log<<std::setprecision(17)<<"t,mean_squared_contrast,mean_contrast,coherence,leakage,activated_fraction,inactive_nonleaked_fraction,quantum_transport,classical_transport,max_abs,field_events,phase_events,leak_events,repair_events,activation_events\n";
    events<<std::setprecision(17)<<"t,site,cause,activity_rate,propagation_rate\n";
    uint64_t nf=0,np=0,nl=0,nr=0,na=0,nbound=0;int remaining=initial_h?0:M;
    double time=0,first=-1,complete=initial_h?0:-1,maxerr=0,maxrateerr=0;
    auto record=[&](double t){
        double E=0,av=0,cs=0,ls=0,hs=0,qn=0,tq=0,tc=0,mx=0;
        for(int i=0;i<M;i++) {
            E+=u[i]*u[i];av+=u[i];cs+=quantum&&!h[i]?s[i]:0;ls+=s[i]==0;hs+=h[i];qn+=!h[i]&&s[i]!=0;
            mx=std::max(mx,std::abs(u[i]));
            for(int k=0;k<4;k++){double v=u[i]-u[nb[i][k]],w=edge[i][k]*v*v;tq+=w*(quantum&&!h[i]?s[i]:0);tc+=w*(h[i]?eff:0)+Dcl*v*v;}
        }
        out.write((char*)&t,8);out.write((char*)n.data(),4*M);out.write((char*)s.data(),M);out.write((char*)h.data(),M);
        log<<t<<","<<E/M<<","<<av/M<<","<<cs/M<<","<<ls/M<<","<<hs/M<<","<<qn/M<<","<<tq/M<<","<<tc/M<<","<<mx<<","<<nf<<","<<np<<","<<nl<<","<<nr<<","<<na<<"\n";
    };
    std::vector<double> saves={0,.05,.1,.2,.5,1,2};for(double t=save;t<=T+1e-8;t+=save)saves.push_back(std::min(t,T));saves.push_back(T);
    std::sort(saves.begin(),saves.end());saves.erase(std::unique(saves.begin(),saves.end()),saves.end());while(!saves.empty()&&saves.back()>T)saves.pop_back();
    size_t save_idx=0;record(0);save_idx=1;
    const auto start=std::chrono::steady_clock::now();
    while(time<T) {
        const double total=tree.total();if(total<=0) {while(save_idx<saves.size())record(saves[save_idx++]);break;}
        const double next=time-std::log(random())/total;
        while(save_idx<saves.size()&&saves[save_idx]<=std::min(next,T))record(saves[save_idx++]);
        if(next>T)break;time=next;
        int j=tree.choose(random()*total);double du=0;bool changed_u=false,changed_h=false;
        if(j<M) {
            int step=random()*(up[j]+down[j])<up[j]?1:-1;n[j]+=step;du=eps*step;changed_u=true;nf++;if(std::abs(n[j])>=cap)nbound++;
        } else if(j<2*M) {
            j-=M;
            if(s[j]==0){s[j]=!quantum||h[j]?1:(random()<.5?-1:1);du=eps/2;nr++;changed_u=true;}
            else {double flip=quantum&&!h[j]?(s[j]<0?r+eta:eta):0;if(random()*(leak+flip)<leak){s[j]=0;du=-eps/2;nl++;changed_u=true;}else{s[j]=-s[j];np++;}}
            update_phase(j);
        } else {
            j-=2*M;if(h[j])throw std::runtime_error("Double environmental activation");
            const bool local=random()*(env_local[j]+env_spread[j])<env_local[j];events<<time<<","<<j<<","<<(local?"activity":"propagation")<<","<<env_local[j]<<","<<env_spread[j]<<"\n";
            h[j]=1;if(s[j]!=0)s[j]=1;na++;remaining--;changed_h=true;
            if(first<0)first=time;if(!remaining)complete=time;update_phase(j);
        }
        if(changed_u) {
            u[j]=eps*(n[j]-(s[j]==0?.5:0));
            for(int k=0;k<13;k++)base[support[j][k]]+=coeff[j][k]*du;
            for(int i:support[j])update_field(i);
            update_environment(j);for(int i:nb[j])update_environment(i);
        } else {
            update_field(j);for(int i:nb[j])update_field(i);
            if(changed_h){update_environment(j);for(int i:nb[j])update_environment(i);}
        }
        if(verify_every>0&&(nf+np+nl+nr+na)%verify_every==0) {
            for(int i=0;i<M;i++){double v=0;for(int k=0;k<13;k++)v+=coeff[i][k]*u[support[i][k]];maxerr=std::max(maxerr,std::abs(v-base[i]));base[i]=v;double old=up[i]+down[i];update_field(i);maxrateerr=std::max(maxrateerr,std::abs(old-up[i]-down[i]));}
        }
    }
    while(save_idx<saves.size())record(saves[save_idx++]);
    if(!out||!log||!events)throw std::runtime_error("Output stream error");
    std::ofstream meta(prefix+"_engine.json");
    meta<<std::setprecision(17)<<"{\"field_events\":"<<nf<<",\"phase_events\":"<<np<<",\"leak_events\":"<<nl<<",\"repair_events\":"<<nr<<",\"activation_events\":"<<na<<",\"boundary_hits\":"<<nbound<<",\"first_activation_time\":"<<first<<",\"all_activated_time\":"<<complete<<",\"remaining_inactive\":"<<remaining<<",\"max_stencil_error\":"<<maxerr<<",\"max_rate_recomputation_error\":"<<maxrateerr<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
