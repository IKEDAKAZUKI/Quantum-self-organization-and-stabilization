#!/usr/bin/env python3
"""Verify input integrity, reported ensembles, generators and short trajectories."""
import os
import sys
import tempfile
sys.dont_write_bytecode=True
os.environ.setdefault('NUMBA_CACHE_DIR',os.path.join(tempfile.gettempdir(),'qtp-numba-cache'))
from qtp.verification import main
if __name__ == '__main__':
    raise SystemExit(main())
