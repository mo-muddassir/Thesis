#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

#Never change these
pi = np.pi
r = np.logspace(-2.9,3)
vr_2 = (1/6)*(1/(1+(r**2)))**0.5
den = (3/(4*pi))*pow((1+(r**2)),-5/2)

#read data
ic_den_un = np.loadtxt('/home/moe/research/Thesis/analysis/unequal_ic.den', unpack=True)
ic_vel_un = np.loadtxt('/home/moe/research/Thesis/analysis/unequal_ic.vel', unpack=True)
ic_den_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.den', unpack=True)
ic_vel_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.vel', unpack=True)

#plot
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den))
axs[1].plot(np.log10(r), np.log10(vr_2))
axs[0].scatter(ic_den_eq[0],ic_den_eq[1], color = 'blue', marker = '^',label = 'Equal Mass')
axs[1].scatter(ic_vel_eq[0],ic_vel_eq[1], color='blue', marker = '^', label = 'Equal Mass')
axs[0].scatter(ic_den_un[0],ic_den_un[1], color='red', marker=',', label='Unequal Mass')
axs[1].scatter(ic_vel_un[0],ic_vel_un[1], color='red', marker=',', label='Unequal Mass')


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
plt.savefig("Comparing_IC_Profiles.pdf")
plt.show()

