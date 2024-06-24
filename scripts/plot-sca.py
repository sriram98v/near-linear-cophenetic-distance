#!/usr/bin/env python
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import os, sys
import tqdm
import matplotlib as mpl
import pandas as pd
import math

mpl.rcParams["figure.figsize"] = (10, 8)
mpl.rcParams["figure.dpi"] = 300
file = "./sca-study"
data = dict()
with open(file, 'r') as f:
    for line in f.readlines():
        run_data = line.split(":")
        run_type = run_data[0].split("_")
        run_time = np.array([float(i) for i in run_data[1].split(",")])
        data["_".join(run_type)] = run_time

fig, ax = plt.subplots(1,1)
for key,val in data.items():
    ax.set_title(f"Scalability analysis for multiple norms")
    ax.plot(range(200, 10001, 200), val, label=f"{key.split("_")[0]} norm={key.split("_")[1]}")
    ax.grid(color='gray', linewidth=0.5, linestyle="--")
fig.legend()
fig.tight_layout()

plt.subplots_adjust(right=0.8) 
plt.savefig("sca_p.png")




