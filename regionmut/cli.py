#!/usr/bin/env python3

import subprocess
import re
import os
import tempfile
import sys
from regionmut.helpers import *

## i am going for the simpler version first
## there might be a way to do this with click
## but i am not sure how to do it

## EXEC example
## python3 regionmut.py bin -m exec -p pkg -o out

# sys.argv is the list of command-line arguments.

def cli():
    if len(sys.argv) == 0:
        raise ValueError("No arguments passed")
    # 0 is script
    mode_sel = sys.argv[1]
    other_args = sys.argv[2:]
    args_str = " ".join(other_args)
    print("mode -> {}".format(mode_sel))
    print("args -> {}".format(args_str))
    path_to_file = get_exec_path_mode(mode_sel, "regionMut")
    print("path to file -> {}".format(path_to_file))
    print(other_args)
    other_args.insert(0, path_to_file)
    print("full_cmd -> {}".format(" ".join(other_args)))
    outprocess = subprocess.run(other_args, check=False, capture_output=True)
    if outprocess.returncode == 0:
        ## need to generate also the stdout?
        print(outprocess.stderr)
    else:
        print(outprocess.stderr)
        sys.exit(outprocess.returncode)

if __name__ == "__main__":
    cli()
