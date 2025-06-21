#!/usr/bin/env python3
import subprocess
import sys

s = sys.argv[1]
if s not in ["b", "d", "e", "f"]:
    raise ValueError("Argument must be [b, d, e, f]")


if s == "b":
    nu = 1
    spot_theta = 90
    obs_theta = 90
elif s == "d":
    nu = 200
    spot_theta = 90
    obs_theta = 90
elif s == "e":
    nu = 400
    spot_theta = 60
    obs_theta = 30
elif s == "f":
    nu = 400
    spot_theta = 20
    obs_theta = 80

subprocess.run("make sd1bdef", shell=True)
subprocess.run(f"./build/sd1bdef {nu} {spot_theta} {obs_theta}", shell=True)
subprocess.run(f"python3 p_sd1.py {s}", shell=True)
