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
r_ad01, rho_ad01,phi_ad01, vr2_ad01, vt2_ad01, rhostar_ad01, phistar_ad01, vr2star_ad01, vt2star_ad01 = np.loadtxt('/home/moe/research/Thesis/outputs/pred_01m.txt', unpack = True)

r_ad001, rho_ad001,phi_ad001, vr2_ad001, vt2_ad001, rhostar_ad001, phistar_ad001, vr2star_ad001, vt2star_ad001 = np.loadtxt('/home/moe/research/Thesis/outputs/pred_001m.txt', unpack = True)

r_ad1, rho_ad1,phi_ad1, vr2_ad1, vt2_ad1, rhostar_ad1, phistar_ad1, vr2star_ad1, vt2star_ad1 = np.loadtxt('/home/moe/research/Thesis/outputs/pred_1m.txt', unpack = True)

#plot
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den))
axs[1].plot(np.log10(r), np.log10(vr_2))

axs[0].plot(np.log10(r_ad001), np.log10(rhostar_ad001), label = 'Small BH (Mbh/Mtot = 0.001)',color='blue',zorder=1)
axs[1].plot(np.log10(r_ad001), np.log10(vr2star_ad001), label ='Small BH (Mbh/Mtot = 0.001)',color='blue',zorder=1)


axs[0].plot(np.log10(r_ad01), np.log10(rhostar_ad01), label = 'Medium BH (Mbh/Mtot = 0.01)',color='green',zorder=2)
axs[1].plot(np.log10(r_ad01), np.log10(vr2star_ad01), label ='Medium BH (Mbh/Mtot = 0.01)',color='green',zorder=2)

axs[0].plot(np.log10(r_ad1), np.log10(rhostar_ad1), label = 'Large BH (Mbh/Mtot = 0.1)',color='red',zorder=3)
axs[1].plot(np.log10(r_ad1), np.log10(vr2star_ad1), label ='Large BH (Mbh/Mtot = 0.1)',color='red',zorder=3)

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
plt.savefig("Adiabatic_Predictions.pdf")
plt.show()
