#!/usr/bin/env python3
from downstream_cli import main
import sys

if __name__ == "__main__":
    sys.argv.insert(1, "pca")
    raise SystemExit(main())
