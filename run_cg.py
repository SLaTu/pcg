import subprocess
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os


os.environ["OMP_SCHEDULE"]="auto"

subprocess.call(["perf", "stat", "./cg.out", sys.argv[1], sys.argv[2], "10000", sys.argv[3], "3", "3", sys.argv[4], sys.argv[5]])

