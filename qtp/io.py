"""Load parameter files and exact saved number/phase trajectories."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json
import numpy as np
from .observables import occupation_contrast

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'


# Defaults are the defaults of src/critical_jump.cpp. Explicit configuration files
# are authoritative: no value is taken from a neighbouring parameter JSON.
ENGINE_DEFAULTS: dict[str, int | float] = {
    "nx": 24, "ny": 28, "mode": 0, "phaseinit": 0, "phase_draw_order": 0,
    "eps": 0.0005, "U": 1.0, "delta": 0.3, "D2": 1.6, "D4": 1.0,
    "zeta": 0.0, "nu": 0.0, "b": 1.0, "c": 0.1, "r": 20.0,
    "eta": 1.0, "leak": 0.0002, "repair": 2.0, "Db": 0.00001,
    "smooth": 0.002, "Dcl": 0.001, "disorder": 0.05,
    "T": 250.0, "save": 2.0, "m0": 0.0,
}
INTEGER_KEYS = {"nx", "ny", "mode", "phaseinit", "phase_draw_order", "seed", "disseed"}
FLOAT_KEYS = (set(ENGINE_DEFAULTS) - INTEGER_KEYS) | {"t1", "t2", "eta_mid"}
METADATA_KEYS = {
    "family", "group", "id", "parent", "block", "rep", "stream", "tf", "E0",
    "initseed", "sigma_init", "engine_sha256", "initial_sha256",
    "configuration_source_sha256", "configuration_source_format", "source_id",
}


def validate_parameter_keys(parameters: dict[str, Any], source: str | Path) -> None:
    """Reject unknown settings before engine defaults are applied.

    CFG and JSON accept exactly the same dynamical and provenance keys.
    Suggestions identify likely spelling errors without changing any input.
    """
    from difflib import get_close_matches
    allowed = INTEGER_KEYS | FLOAT_KEYS | METADATA_KEYS | {"init"}
    unknown = sorted(set(parameters) - allowed)
    if unknown:
        details = []
        for key in unknown:
            nearest = get_close_matches(key, sorted(allowed), n=1, cutoff=0.6)
            hint = f" (did you mean {nearest[0]!r}?)" if nearest else ""
            details.append(repr(key) + hint)
        raise ValueError(f"{source}: unknown configuration key(s): " + ", ".join(details))


def parameter_file(run: str | Path) -> Path:
    """Resolve an explicit file or a recorded run ID to one parameter source.

    An explicit .cfg or .json is read as supplied. A run ID/prefix or .npz path
    uses ID_params.json. Ambiguous IDs and missing explicit files raise errors.
    """
    p = Path(run).expanduser()
    if p.suffix.lower() in {".cfg", ".json"}:
        if not p.is_file():
            raise FileNotFoundError(f"Parameter file not found: {p}")
        return p.resolve()
    if p.suffix.lower() == ".npz":
        p = p.with_suffix("")
    candidate = p.with_name(p.name + "_params.json")
    if candidate.is_file():
        return candidate.resolve()
    if p.parent != Path("."):
        raise FileNotFoundError(f"Parameter JSON not found: {candidate}")
    found = [DATA / group / (p.name + "_params.json")
             for group in ("primary", "preparations", "dephasing")
             if (DATA / group / (p.name + "_params.json")).is_file()]
    if len(found) != 1:
        raise FileNotFoundError(f"Run ID {run!r} has {len(found)} matches; give an explicit file.")
    return found[0].resolve()


def _prefix(source: Path) -> Path:
    name = source.name
    suffix = "_params.json" if name.endswith("_params.json") else source.suffix
    return source.with_name(name[:-len(suffix)])


def run_prefix(run: str | Path) -> Path:
    """Resolve the trajectory prefix associated with a parameter source."""
    return _prefix(parameter_file(run))


def parse_config(path: str | Path) -> dict[str, Any]:
    """Read the engine's ``key value`` configuration syntax.

    Blank lines and full-line # comments are accepted. The complete value of
    ``init`` is a path, including any spaces. Numeric values must be finite;
    seeds are parsed as exact unsigned 64-bit integers. Duplicates and unknown
    keys raise errors so that misspelled simulation settings cannot be ignored.
    """
    from decimal import Decimal, InvalidOperation
    path = Path(path)
    result: dict[str, Any] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        pair = text.split(maxsplit=1)
        if len(pair) != 2:
            raise ValueError(f"{path}:{number}: missing value for {pair[0]!r}")
        key, value = pair
        if key in result:
            raise ValueError(f"{path}:{number}: duplicate configuration key {key!r}")
        validate_parameter_keys({key: value}, f"{path}:{number}")
        if key in INTEGER_KEYS:
            try:
                d = Decimal(value)
                if not d.is_finite() or d != d.to_integral_value():
                    raise ValueError
                parsed = int(d)
                if key in {"seed", "disseed"} and not 0 <= parsed < 2**64:
                    raise ValueError
            except (InvalidOperation, ValueError, OverflowError):
                raise ValueError(f"{path}:{number}: invalid integer {key}={value!r}") from None
        elif key in FLOAT_KEYS:
            try:
                parsed = float(value)
                if not np.isfinite(parsed):
                    raise ValueError
            except ValueError:
                raise ValueError(f"{path}:{number}: invalid numeric value {key}={value!r}") from None
        else:
            # Metadata is not used to fill in missing dynamical parameters.
            parsed = value
            if key in {"E0", "sigma_init", "tf"}:
                parsed = float(value)
            elif key in {"initseed", "block", "rep", "stream"}:
                parsed = int(value)
        result[key] = parsed
    return result


def load_parameters(run: str | Path) -> tuple[Path, dict[str, Any]]:
    """Read exactly the selected source; IDs/prefixes select parameter JSON."""
    source = parameter_file(run)
    if source.suffix.lower() == ".cfg":
        p = parse_config(source)
    else:
        def unique_keys(pairs):
            result = {}
            for key, value in pairs:
                if key in result:
                    raise ValueError(f"{source}: duplicate JSON key {key!r}")
                result[key] = value
            return result
        def reject_nonfinite(value):
            raise ValueError(f"{source}: non-finite JSON value {value!r}")
        p = json.loads(source.read_text(encoding="utf-8"), object_pairs_hook=unique_keys,
                       parse_constant=reject_nonfinite)
        if not isinstance(p, dict):
            raise ValueError(f"{source}: parameters must be a JSON object")
    validate_parameter_keys(p, source)
    return _prefix(source), p


def simulation_parameters(run: str | Path) -> tuple[Path, dict[str, Any], Path]:
    """Load run settings, apply engine defaults, and resolve the initialization.

    ``init`` is relative to the explicitly selected file. JSON files without
    ``init`` use ID_init.txt beside that JSON, as recorded datasets do. An
    explicit CFG must contain init; it never borrows initialization from JSON.
    """
    source = parameter_file(run)
    prefix, supplied = load_parameters(source)
    for key in ("seed", "disseed"):
        if key not in supplied:
            raise ValueError(f"{source}: missing required parameter {key!r}")
    if source.suffix.lower() == ".cfg" and "init" not in supplied:
        raise ValueError(f"{source}: missing required parameter 'init'")
    config = ENGINE_DEFAULTS | supplied
    for key in INTEGER_KEYS:
        value = config[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{source}: {key} must be an integer")
        if not np.isfinite(value) or int(value) != value:
            raise ValueError(f"{source}: {key} must be a finite integer")
        config[key] = int(value)
        if key in {"seed", "disseed"} and not 0 <= config[key] < 2**64:
            raise ValueError(f"{source}: {key} must fit an unsigned 64-bit integer")
    for key in FLOAT_KEYS & config.keys():
        try:
            config[key] = float(config[key])
        except (ValueError, TypeError):
            raise ValueError(f"{source}: {key} must be numeric") from None
        if not np.isfinite(config[key]):
            raise ValueError(f"{source}: {key} must be finite")
    initial_value = config.get("init", prefix.name + "_init.txt")
    if not isinstance(initial_value, str) or not initial_value.strip():
        raise ValueError(f"{source}: init must be a nonempty path string")
    initial = Path(initial_value).expanduser()
    if not initial.is_absolute():
        initial = source.parent / initial
    initial = initial.resolve()
    if not initial.is_file():
        raise FileNotFoundError(f"Initial configuration not found: {initial}")
    return prefix, config, initial


@dataclass(frozen=True)
class Trajectory:
    """Arrays with time first and spatial axes ordered as (ny, nx)."""
    t: np.ndarray
    n: np.ndarray
    s: np.ndarray
    eps: float

    @property
    def u(self) -> np.ndarray:
        return occupation_contrast(self.n, self.s, self.eps)


def load_trajectory(run: str | Path) -> tuple[dict[str, Any], Trajectory]:
    if Path(run).suffix.lower() == ".cfg":
        raise ValueError("Read a recorded trajectory by ID, NPZ path or parameter JSON; "
                         "run an edited CFG with run_simulation.py before analyzing its outputs.")
    prefix, p = load_parameters(run)
    with np.load(str(prefix)+'.npz', allow_pickle=False) as z:
        t, n, s = (z[key].copy() for key in ('t', 'n', 's'))
    if n.shape != s.shape or n.shape != (len(t), int(p['ny']), int(p['nx'])):
        raise ValueError(f'Array dimensions disagree with parameters for {prefix.name}')
    if len(t) == 0 or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError(f'Invalid saved times for {prefix.name}')
    if not np.issubdtype(n.dtype, np.integer) or not np.all(np.isin(s, [-1, 0, 1])):
        raise ValueError(f'Invalid occupation/phase arrays for {prefix.name}')
    return p, Trajectory(t=t, n=n, s=s, eps=float(p['eps']))


def binary_to_npz(binary: Path, output: Path, eps: float) -> dict[str, Any]:
    """Convert the C++ little-endian stream to t/n/s/eps NPZ without changing values.

    Header: int32 nx, ny, nx*ny. Each record: float64 t, int32 n[M], int8 s[M].
    A truncated record raises ValueError. Existing output files are not overwritten.
    """
    binary, output = Path(binary), Path(output)
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    if binary.stat().st_size < 12:
        raise ValueError('Binary trajectory is missing its 12-byte header.')
    with binary.open('rb') as stream:
        nx, ny, m = map(int, np.fromfile(stream, dtype='<i4', count=3))
        if nx < 1 or ny < 1 or m != nx*ny:
            raise ValueError('Invalid binary dimensions.')
        record = np.dtype([('t', '<f8'), ('n', '<i4', (m,)), ('s', 'i1', (m,))])
        if (binary.stat().st_size-12) % record.itemsize:
            raise ValueError('Truncated binary trajectory record.')
        values = np.fromfile(stream, record)
    if not len(values) or not np.all(np.isfinite(values['t'])) or np.any(np.diff(values['t']) <= 0):
        raise ValueError('Empty or invalid binary trajectory times.')
    if not np.all(np.isin(values['s'], [-1, 0, 1])):
        raise ValueError('Invalid labels in binary trajectory.')
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, t=values['t'].copy(),
                        n=values['n'].reshape(-1, ny, nx).copy(),
                        s=values['s'].reshape(-1, ny, nx).copy(), eps=float(eps))
    return dict(saved_frames=len(values), array_shape=[ny, nx], first_time=float(values['t'][0]),
                last_time=float(values['t'][-1]), output=str(output))
