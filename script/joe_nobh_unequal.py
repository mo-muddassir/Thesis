#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('xtick', labelsize=13)
matplotlib.rc('ytick', labelsize=13)

#Never change these
pi = np.pi
r = np.logspace(-2.9,3)
vr_2 = (1/6)*(1/(1+(r**2)))**0.5
den = (3/(4*pi))*pow((1+(r**2)),-5/2)

#read data
ic_den_joe = np.loadtxt('/home/moe/research/Thesis/analysis/joe_IC.den', unpack=True)
ic_vel_joe = np.loadtxt('/home/moe/research/Thesis/analysis/joe_IC.vel', unpack=True)
joe_nobh_den = np.loadtxt('/home/moe/research/Thesis/analysis/joe_nobh_unequal.den', unpack=True)
joe_nobh_vel = np.loadtxt('/home/moe/research/Thesis/analysis/joe_nobh_unequal.vel', unpack=True)


#plot
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den), zorder=1, linewidth=2.5, color='black', label='Theoretical ICs')
axs[1].plot(np.log10(r), np.log10(vr_2), zorder=1, linewidth=2.5, color='black', label='Theoretical ICs')


axs[0].scatter(ic_den_joe[0],ic_den_joe[1], zorder=2, marker='.', color='blue', label='Sim ICs')
axs[1].scatter(ic_vel_joe[0],ic_vel_joe[1], zorder=2, marker='.', color='blue', label='Sim ICs')

axs[0].scatter(joe_nobh_den[0],joe_nobh_den[1], zorder=3, marker='.', color='red', label = 'Sim noBH')
axs[1].scatter(joe_nobh_vel[0],joe_nobh_vel[1], zorder=3, marker='.', color='red', label='Sim noBH')



#plot details
#axs[0].legend(frameon=False, fontsize=12, prop={'family':'serif'})
#axs[1].legend(frameon=False, fontsize=12, prop = {'family':'serif'})
axs[0].set_xlim(-3,3)
axs[1].set_xlim(-3,3)
axs[0].set_ylim(-15,2)
axs[1].set_ylim(-4.1,0)
axs[0].legend()
axs[1].legend()

axs[1].set_xlabel(r"$log_{10}(r)$", fontsize=16,fontweight='bold', family='serif')
axs[0].set_ylabel(r"$log_{10}(rho)$", fontsize=16, fontweight = 'bold', family='serif')
axs[1].set_ylabel(r"$log_{10}(<v_r^2>)$", fontsize=16, family ='serif')
plt.savefig("Joe_noBH_Unequal.pdf")
plt.show()

