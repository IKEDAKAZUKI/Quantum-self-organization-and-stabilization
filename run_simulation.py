#!/usr/bin/env python3
"""Run an explicit CFG/JSON or a recorded trajectory configuration and write BIN, CSV and NPZ outputs."""
import argparse
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from qtp.simulation import simulate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',help='Run ID/prefix (selects JSON), explicit .cfg, or parameter JSON')
    parser.add_argument('--output',type=Path,required=True,help='Empty results directory')
    parser.add_argument('--until',type=float,help='Shorter terminal time; omitted for the full trajectory')
    args=parser.parse_args()
    try:
        report=simulate(args.run,args.output,args.until)
    except (ValueError,FileNotFoundError,FileExistsError,RuntimeError) as exc:
        parser.exit(1,f'run_simulation: {exc}\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
