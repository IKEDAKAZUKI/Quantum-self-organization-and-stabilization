"""Compile the C++ generator, run portable configurations, and convert trajectories."""
from __future__ import annotations
from pathlib import Path
from typing import Any
import hashlib
import json
import os
import shutil
import subprocess
import sys
from .io import ROOT, parameter_file, simulation_parameters, binary_to_npz
from .safety import output_directory


def compile_engine(output: Path, compiler: str | None = None) -> Path:
    """Build a C++17 executable in an explicit results directory."""
    output=output_directory(output)
    compiler=compiler or os.environ.get('CXX','g++')
    found=shutil.which(compiler)
    if found is None:
        raise FileNotFoundError(f'C++17 compiler not found: {compiler}. Set CXX or install a compiler.')
    engine=output/'critical_jump'
    result=subprocess.run([found,'-O3','-std=c++17',str(ROOT/'src/critical_jump.cpp'),'-o',str(engine)],
                          text=True,capture_output=True)
    if result.returncode:
        raise RuntimeError('C++ compilation failed:\n'+result.stderr)
    return engine


def simulate(run: str | Path, output: Path, until: float | None = None,
             engine: Path | None = None) -> dict[str, Any]:
    """Run the selected CFG/JSON settings and their chosen initialization.

    Omit ``until`` for the entire recorded duration. A shorter duration is useful
    for a prefix check. Outputs include a portable config, init, params, BIN, CSV,
    event diagnostics and converted NPZ. Existing results are never overwritten.
    """
    source=parameter_file(run)
    prefix,p,initial=simulation_parameters(source)
    output=output_directory(output)
    name=prefix.name
    if any((output/(name+suffix)).exists() for suffix in ('.cfg','.npz','.bin','.csv','_params.json','_init.txt')):
        raise FileExistsError('Output already contains this run; choose an empty directory.')
    if sys.byteorder != 'little':
        raise RuntimeError('The recorded binary format currently requires a little-endian platform.')
    config=dict(p)
    original_id = config.get('id')
    if original_id is not None and original_id != name:
        config['source_id'] = original_id
    config['id'] = name
    # The contrast budget follows the actual chosen initialization, including
    # its leakage half-step, even if a user changes init or eps.
    import numpy as np
    initial_values = np.loadtxt(initial, dtype=np.int64, ndmin=2)
    if initial_values.shape != (int(config['nx'])*int(config['ny']), 2):
        raise ValueError(f'Initial configuration dimensions do not match nx*ny: {initial}')
    leaked = (initial_values[:,1] == 0) if int(config['phaseinit']) == 0 else np.zeros(len(initial_values), bool)
    contrast = float(config['eps']) * (initial_values[:,0] - .5*leaked)
    config['E0'] = float(np.mean(contrast**2))
    if until is not None:
        if not 0 <= until <= float(p['T']):
            raise ValueError('until must lie between zero and the recorded endpoint T.')
        config['T']=float(until)
    shutil.copyfile(initial,output/(name+'_init.txt'))
    config['init']=name+'_init.txt'
    (output/(name+'.cfg')).write_text('\n'.join(f'{key} {val}' for key,val in config.items()
                                            if isinstance(val,(int,float,str)))+'\n')
    source_hash=hashlib.sha256((ROOT/'src/critical_jump.cpp').read_bytes()).hexdigest()
    params=config | {'engine_sha256':source_hash,'initial_sha256':hashlib.sha256(initial.read_bytes()).hexdigest(),
                     'configuration_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                     'configuration_source_format':source.suffix.lower()[1:]}
    (output/(name+'_params.json')).write_text(json.dumps(params,indent=2)+'\n')
    if engine is None:
        engine=compile_engine(output/'build')
    completed=subprocess.run([str(Path(engine).resolve()),str(output/(name+'.cfg')),str(output/name)],
                             cwd=output.parent,text=True,capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or 'The trajectory generator failed.')
    converted=binary_to_npz(output/(name+'.bin'),output/(name+'.npz'),float(config['eps']))
    return dict(id=name,generator_sha256=source_hash,status='passed',**converted,
                diagnostics=json.loads((output/(name+'_engine.json')).read_text()))
