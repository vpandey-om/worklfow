#!/usr/bin/env python3
import sys
from downstream_cli import main

if __name__ == "__main__":
    sys.argv.insert(1, "differential-expression")
    raise SystemExit(main())
