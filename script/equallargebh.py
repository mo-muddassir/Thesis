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
den_eq_lgbh = np.loadtxt('/home/moe/research/Thesis/analysis/equal_large_bh_final1.den', unpack=True)
vel_eq_lgbh = np.loadtxt('/home/moe/research/Thesis/analysis/equal_large_bh_final1.vel', unpack=True)

ic_den_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.den', unpack=True)
ic_vel_eq = np.loadtxt('/home/moe/research/Thesis/analysis/equal_ic.vel', unpack=True)

r_ad, rho_ad,phi_ad, vr2_ad, vt2_ad, rhostar_ad, phistar_ad, vr2star_ad, vt2star_ad = np.loadtxt('/home/moe/research/Thesis/outputs/pred_1m.txt', unpack = True)

#plot
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den))
axs[1].plot(np.log10(r), np.log10(vr_2))


axs[0].scatter(ic_den_eq[0],ic_den_eq[1], color = 'blue', marker = '.',label = 'No BH')
axs[1].scatter(ic_vel_eq[0],ic_vel_eq[1], color='blue', marker = '.', label = 'No BH')


axs[0].scatter(den_eq_lgbh[0],den_eq_lgbh[1], color='red', marker='.', label='Large BH (Mbh/M = 0.1)')
axs[1].scatter(vel_eq_lgbh[0],vel_eq_lgbh[1], color='red', marker='.', label='Large BH (Mbh/M = 0.1')

axs[0].plot(np.log10(r_ad), np.log10(rhostar_ad), label = 'Adiabatic Prediction',color='blue',zorder=2,linestyle='dashed')
axs[1].plot(np.log10(r_ad), np.log10(vr2star_ad), label ='Adiabatic Prediction',color='blue',zorder=2,linestyle='dashed')


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
plt.savefig("UMD_Large_BH_Profiles.pdf")
plt.show()

