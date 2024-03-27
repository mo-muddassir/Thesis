#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

#Never change these
pi = np.pi
r = np.logspace(-2.9,3)
vr_2 = (1/6)*(1/(1+(r**2)))**0.5
den = (3/(4*pi))*pow((1+(r**2)),-5/2)

#read data
den_eq_smbh = np.loadtxt('/home/moe/research/Thesis/analysis/equal_small_bh_final.den', unpack=True)
vel_eq_smbh = np.loadtxt('/home/moe/research/Thesis/analysis/equal_small_bh_final.vel', unpack=True)
ic_den_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.den', unpack=True)
ic_vel_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.vel', unpack=True)

#plot
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den))
axs[1].plot(np.log10(r), np.log10(vr_2))
axs[0].scatter(ic_den_eq[0],ic_den_eq[1], color = 'blue', marker = '.',label = 'No BH')
axs[1].scatter(ic_vel_eq[0],ic_vel_eq[1], color='blue', marker = '.', label = 'No BH')
axs[0].scatter(den_eq_smbh[0],den_eq_smbh[1], color='red', marker='.', label='Small BH (Mbh/M = 0.001)')
axs[1].scatter(vel_eq_smbh[0],vel_eq_smbh[1], color='red', marker='.', label='Small BH (Mbh/M = 0.001')


#plot details
axs[0].legend(frameon=False, fontsize=12, prop={'family':'serif'})
axs[1].legend(frameon=False, fontsize=12, prop = {'family':'serif'})
axs[0].set_xlim(-3,3)
axs[1].set_xlim(-3,3)
axs[0].set_ylim(-15,2)
axs[1].set_ylim(-4.1,0)
axs[1].set_xlabel(r"$log_{10}(r)$", fontsize=16,fontweight='bold', family='serif')
axs[0].set_ylabel(r"$log_{10}(rho)$", fontsize=16, fontweight = 'bold', family='serif')
axs[1].set_ylabel(r"$log_{10}(<v_r^2>)$", fontsize=16, family ='serif')
plt.savefig("Small_BH_Profiles.pdf")
plt.show()

