#!/usr/bin/env python3
"""Recreate Figures 1--4 from the stored configurations, ensemble summaries and model equations."""
import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
sys.dont_write_bytecode=True
os.environ.setdefault('MPLCONFIGDIR',os.path.join(tempfile.gettempdir(),'qtp-matplotlib'))
from qtp.figures import render_all


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis',type=Path,default=Path('results/analysis'))
    p.add_argument('--output',type=Path,default=Path('results/figures'))
    args=p.parse_args();report=render_all(args.analysis,args.output)
    print('Recreated all four numerical figures in',args.output)
    return report
if __name__=='__main__':
    main()
