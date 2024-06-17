#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["figure.figsize"] = (10, 8)
mpl.rcParams["figure.dpi"] = 300


dist_types = ["uniform", "yule"]
zeta_types = ["depth", "height", "size"]
n_bins=100

file = "./emp-study"
dists = {dist: {zeta_type: dict() for zeta_type in zeta_types} for dist in dist_types}
with open(file, 'r') as f:
    for line in f.readlines():
        run_data = line.split(":")
        run_info = run_data[0].split(",")
        tree_dist = run_info[0]
        zeta_type = run_info[1]
        norm = run_info[2]
        dist = np.array([float(i) for i in run_data[1].split(",")])
        dists[tree_dist][zeta_type][int(norm)] = dist

fig, ax = plt.subplots(3, 2)
for col,model in enumerate(dist_types):
    for row,zeta in enumerate(zeta_types):
        ax[row,col].hist(dists[model][zeta][1], bins=n_bins)
        ax[row,col].set_title(f"Sampled distibution for {zeta} under {model} model")
        ax[row,col].tick_params(axis='x', labelsize=10)
        ax[row,col].tick_params(axis='y', labelsize=10)
        ax[row,col].grid(color='gray', linewidth=0.5, linestyle="--")
fig.tight_layout()
plt.savefig("distribs.png")

fig, ax = plt.subplots(3, 2)
for col,model in enumerate(dist_types):
    for row,zeta in enumerate(zeta_types):
        for norm in range(1, max(dists[model][zeta].keys())):
            ax[row,col].hist(dists[model][zeta][norm], bins=n_bins)
        ax[row,col].set_title(f"Sampled distibution for {zeta} \n under {model} model for multiple norms")
        ax[row,col].tick_params(axis='x', labelsize=10)
        ax[row,col].tick_params(axis='y', labelsize=10)
        ax[row,col].set_xscale('log', base=10)
        ax[row,col].grid(color='gray', linewidth=0.5, linestyle="--")
fig.legend([f"p={i}" for i in range(1, max(dists[model][zeta].keys()))], loc="upper right")
fig.tight_layout()

plt.subplots_adjust(right=0.9) 
plt.savefig("distribs_p.png")




