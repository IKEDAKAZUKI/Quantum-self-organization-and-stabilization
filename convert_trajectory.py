#!/usr/bin/env python3
"""Convert a C++ binary trajectory to NPZ without modifying the input."""
import argparse,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
from qtp.io import binary_to_npz
from qtp.safety import output_directory
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary',type=Path);p.add_argument('--eps',type=float,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    output_directory(args.output.parent)
    print(json.dumps(binary_to_npz(args.binary,args.output,args.eps),indent=2))
