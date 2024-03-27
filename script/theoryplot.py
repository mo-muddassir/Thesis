#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
r = np.logspace(-2.9,3)
vr_2 = (1/6)*(1/(1+(r**2)))**0.5

den = (3/(4*pi))*pow((1+(r**2)),-5/2)
fig, axs = plt.subplots(2, figsize = (12,12), constrained_layout=True)
axs[0].plot(np.log10(r),np.log10(den))
axs[1].plot(np.log10(r), np.log10(vr_2))
#axs[0].legend(frameon=False, fontsize=12, prop={'family':'serif'})
#axs[1].legend(frameon=False, fontsize=12, prop = {'family':'serif'})
axs[0].set_xlim(-3,3)
axs[1].set_xlim(-3,3)
axs[0].set_ylim(-15,2)
axs[1].set_ylim(-4.1,0)
axs[1].set_xlabel(r"$log_{10}(r)$", fontsize=24,fontweight='bold', family='serif')
axs[0].set_ylabel(r"$log_{10}(rho)$", fontsize=24, fontweight = 'bold', family='serif')
axs[1].set_ylabel(r"$log_{10}(<v_r^2>)$", fontsize=24, family ='serif')
plt.savefig("Theoretical Plummer Profiles.pdf")
plt.show()

