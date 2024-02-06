import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

#rc('text.latex',preamble='\\usepackage{libertine}\n\\usepackage[libertine]{newtxmath}')
#rc('font',**{'family':'serif','serif':['Linux Libertine O']}, size=18)
#rc('text', usetex=True)

data = []
for f in ["data_m0.0000.txt", "data_m0.0001.txt", "data_m0.0010.txt",  "data_m0.0100.txt", "0.01m_e0.05.den", "0.01m_e0.05.vel"]:
	data.append(np.loadtxt(f))

fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(2,2,1)
ax.loglog(data[0][:,0], data[0][:,5], color='k')
ax.loglog(data[1][:,0], data[1][:,5], color='k')
ax.loglog(data[2][:,0], data[2][:,5], color='k')
ax.loglog(data[3][:,0], data[3][:,5], color='k')
ax.loglog(data[4][:,0],data[4][:,1],color = 'blue', marker='o')
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\rho(r)$")
ax.set_xlim(1e-4, 300)
ax.set_ylim(1e-13, 5000)

ax = fig.add_subplot(2,2,2)
ax.semilogx(data[0][:,0], data[0][:,6], color='k')
ax.semilogx(data[1][:,0], data[1][:,6], color='k')
ax.semilogx(data[2][:,0], data[2][:,6], color='k')
ax.semilogx(data[3][:,0], data[3][:,6], color='k')
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\Phi(r)$")
ax.set_xlim(1e-4, 300)
ax.set_ylim(-4, 0.2)

ax = fig.add_subplot(2,2,3)

ax.loglog(data[0][:,0], data[0][:,7], color='k')
ax.loglog(data[1][:,0], data[1][:,7], color='k')
ax.loglog(data[2][:,0], data[2][:,7], color='k')
ax.loglog(data[3][:,0], data[3][:,7], color='k')
ax.loglog(data[5][:,0],data[5][:,1],color = 'blue', marker = 'o')
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\langle v_r^2 \rangle (r)$")
ax.set_xlim(1e-4, 300)
ax.set_ylim(1e-3, 100)

ax = fig.add_subplot(2,2,4)
ax.loglog(data[0][:,0], data[0][:,8], color='k')
ax.loglog(data[1][:,0], data[1][:,8], color='k')
ax.loglog(data[2][:,0], data[2][:,8], color='k')
ax.loglog(data[3][:,0], data[3][:,8], color='k')
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\langle v_t^2 \rangle (r)$")
ax.set_xlim(1e-4, 300)
ax.set_ylim(1e-3, 100)

plt.tight_layout()
plt.savefig("fig_results.pdf", bbox_inches="tight")
