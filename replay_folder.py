#python exec.py -o resumeoffline OUTDIR FRAMENUM

import os
import sys
import subprocess
import argparse


# get command line arguments
ap = argparse.ArgumentParser()
ap.add_argument("folder", type=str)

args = ap.parse_args()
args = vars(args)

workdir = os.getcwd()
executable = os.path.join(os.getcwd(), "build-Release", "bin", "arcsim_0.2.1")

try:
    for conf in os.listdir(args["folder"]):
        print("STARTING:", conf)
        conffile = os.path.join(args["folder"],conf)
        subprocess.check_call([executable, "replay", conffile], cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
