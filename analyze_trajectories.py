#!/usr/bin/env python3
"""Recompute the article's observables and summaries from stored trajectories."""
import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
sys.dont_write_bytecode=True
os.environ.setdefault('NUMBA_CACHE_DIR',os.path.join(tempfile.gettempdir(),'qtp-numba-cache'))
from qtp.analysis import analyze_all,analyze_primary,analyze_preparations,analyze_dephasing


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',choices=['all','primary','preparations','dephasing'],default='all')
    p.add_argument('--output',type=Path,default=Path('results/analysis'))
    args=p.parse_args()
    action={'all':analyze_all,'primary':analyze_primary,'preparations':analyze_preparations,'dephasing':analyze_dephasing}[args.dataset]
    print(json.dumps(action(args.output),indent=2))
if __name__=='__main__':main()
