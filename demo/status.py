#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from ase.io import read
iters = [2,3]
trajs = {}
trajl = {}
for it in iters:
    trajs[it] = [read(glob(f'../trajs/iter{it}/seed{i}/*.xyz')[0], index=':')
                 for i in range(1,4)]
    trajl[it] = read(f'../trajs/iter{it}/label.xyz', index=':')

fig, axs = plt.subplots(3,1, figsize=[5,5], sharex=True)
t0 = 0
dt = 0.005
getE = lambda atoms: atoms.get_potential_energy()
getF = lambda atoms: np.abs(atoms.get_forces()).max()
getD = lambda atoms: atoms.get_masses().sum()/6.022e23/atoms.get_volume()*1e24

for it in iters:
    t=t0+np.arange(len(trajs[it][0]))*dt
    for i, traj in enumerate(trajs[it]):
        style = {'color':'blue' if i==0 else 'grey',
                 'label':f'S{i+1} {"[md]" if i==0 else "[qbc]"}'}
        axs[0].plot(t, list(map(getE, traj)), **style)
        axs[1].plot(t, list(map(getF, traj)), **style)
        axs[2].plot(t, list(map(getD, traj)), **style)

    style = {'color':'red', 'label':'lammps'}
    t=t0+np.arange(len(trajl[it]))*dt
    axs[0].plot(t, list(map(getE, trajl[it])), **style)
    axs[1].plot(t, list(map(getF, trajl[it])), **style)
    axs[2].plot(t, list(map(getD, trajl[it])), **style)
    t0+=len(trajs[it][0])*dt
    if it==2: axs[0].legend()

axs[0].set_ylabel('Energy [eV]')
axs[1].set_ylabel('Force [eV/$\AA$]')
axs[1].set_ylim(0,20)
axs[2].set_ylabel('Density [g/$cm^3$]')
fig.align_ylabels(axs)
plt.tight_layout()
plt.savefig('status.png')
