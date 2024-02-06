#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from numba import njit, vectorize, float64

rc('text.latex',preamble='\\usepackage{libertine}\n\\usepackage[libertine]{newtxmath}')
rc('font',**{'family':'serif','serif':['Linux Libertine O']}, size=18)
rc('text', usetex=True)

# we'll need this later, and it's a little better than trapz
@njit
def int_simpson(f, x):
	
	return np.trapz(f, x)
	
	"""
	N = len(x) - 1
	h = np.diff(x)
		   
	result = 0.0
	for i in range(1, N, 2):
		if h[i] != 0.0:
			hph = h[i] + h[i - 1]
			result += f[i] * ( h[i]**3 + h[i - 1]**3 + 3. * h[i] * h[i - 1] * hph ) / ( 6 * h[i] * h[i - 1] )
			result += f[i - 1] * ( 2. * h[i - 1]**3 - h[i]**3 + 3. * h[i] * h[i - 1]**2) / ( 6 * h[i - 1] * hph)
			result += f[i + 1] * ( 2. * h[i]**3 - h[i - 1]**3 + 3. * h[i - 1] * h[i]**2) / ( 6 * h[i] * hph )

	if (N + 1) % 2 == 0:
		if h[N] != 0.0:
			result += f[N] * ( 2 * h[N - 1]**2 + 3. * h[N - 2] * h[N - 1]) / ( 6 * ( h[N - 2] + h[N - 1] ) )
			result += f[N - 1] * ( h[N - 1]**2 + 3*h[N - 1]* h[N - 2] ) / ( 6 * h[N - 2] )
			result -= f[N - 2] * h[N - 1]**3 / ( 6 * h[N - 2] * ( h[N - 2] + h[N - 1] ) )
	return result
	"""
	
# remember M = a = 1!
def density(r):
	return 3 / 4 / np.pi * np.power(1 + r**2, -5/2)

def potential(r):
	return -np.power(1 + r**2, -1/2)

def average_radial_velocity_squared(r):
	return 1.0/6.0 * np.power(1 + r**2, -1/2)

def average_tangential_velocity_squared(r):
	return 1.0/3.0 * np.power(1 + r**2, -1/2)

@vectorize([float64(float64, float64)])
def distribution_function(E, J):
	return 24 * np.sqrt(2) / 7 / np.pi**3 * np.power(-E, 7/2)

@njit
def calc_density(r, phi, f):
	"""
	Calculate the density of the system on a grid based on the potential and 
	distribution function.
	"""
	
	N = f.shape[0]
	Nj = f.shape[1]
	
	E = phi
	
	rho = np.zeros(N)
	
	# loop over radius
	for i in range(N):
		
		# placeholder for energy integral
		h = np.zeros(N)
		
		# loop over energy
		for j in range(N):
			
			jmax2 = 2 * r[i]**2 * (E[j] - phi[i])
			if jmax2 > 0:
				jmax = np.sqrt(jmax2)
			else:
				jmax = 0
			
			J = np.linspace(0, jmax, Nj)
			
			g = np.zeros(Nj)
			for k in range(Nj):
				vr = 1.0 / r[i] * np.sqrt(jmax**2 - J[k]**2)
				if vr > 0:
					g[k] = J[k] / r[i] / r[i] / vr * f[j, k]
				else:
					g[k] = 0.0
					
			# do the angular momentum integral
			h[j] = int_simpson(g, J)#np.trapz(g, J)
		
		# then the energy integral
		rho[i] = 4.0 * np.pi * int_simpson(h, E)
	
	return rho
	
@njit
def calc_potential(r, rho, mass):
	"""
	Calculate the density of the system on a grid based on the potential and 
	distribution function.
	"""

	N = len(r)

	rho_r = rho * r
	rho_r2 = rho * r * r

	phi = np.zeros(N)
	for i in range(N):
		g1 = int_simpson(rho_r2[0:i], r[0:i])
		g2 = int_simpson(rho_r[i:N-1], r[i:N-1])
				
		phi[i] = -4.0 * np.pi * (g1 / r[i] + g2) - mass / r[i]

	return phi
	
@njit
def calc_radial_action(r, phi):
	"""
	Calculate the radial action of the system on a grid.
	"""   

	E = phi
	vr = np.zeros(N)
	Ir = np.zeros((N, Nj))

	# energy loop
	for j in range(N):

		# first find the largest value of Jmax with this energy by looping
		# over all radii
		jmax_max = 0.0
		for i in range(N):
			jmax2 = 2 * r[i]**2 * (E[j] - phi[i])
			if jmax2 > jmax_max**2:
				jmax_max = np.sqrt(jmax2)

		# now set the angular momentum to go from 0 to jmax_max
		J = np.linspace(0, jmax_max, Nj)

		# angular momentum loop
		for k in range(Nj):

			# radius loop to perform integral
			for i in range(N):
				vr2 = 2 * (E[j] - phi[i]) - J[k]**2 / r[i]**2
				if vr2 > 0:
					vr[i] = np.sqrt(vr2)
				else:
					vr[i] = 0
			Ir[j, k] = 2 * int_simpson(vr, r)

	return Ir

@njit
def calc_velocity_moment(r, rho, phi, f, m, n):
	"""
	Calculate the velocity moment <v_r^m v_t^n> of the system on a grid.
	"""

	N = f.shape[0]
	Nj = f.shape[1]

	E = phi
	vel = np.zeros(N)

	# loop over radius
	for i in range(N):

		# placeholder for energy integral
		h = np.zeros(N)

		# loop over energy
		for j in range(N):

			jmax2 = 2 * r[i]**2 * (E[j] - phi[i])
			if jmax2 > 0:
				jmax = np.sqrt(jmax2)
			else:
				jmax = 0

			J = np.linspace(0, jmax, Nj)

			g = np.zeros(Nj)
			for k in range(Nj):
				vr = 1.0 / r[i] * np.sqrt(jmax**2 - J[k]**2)
				vt = J[k] / r[i]
				if vr > 0:
					g[k] = J[k] / r[i] / r[i] / vr * f[j, k] * np.power(vr, m) * np.power(vt, n)
				else:
					g[k] = 0.0


			# do the angular momentum integral
			h[j] = int_simpson(g, J)

		# then the energy integral
		if rho[i] > 0.0:
			vel[i] = 4.0 * np.pi * int_simpson(h, E) / rho[i]
		else:
			vel[i] = 0.0

	return vel

@njit
def calc_dfstar(phi, ir, irstar):

	N = ir.shape[0]
	Nj = ir.shape[1]

	E = phi

	fstar = np.zeros((N, Nj))

	for i in range(N):
		for k in range(Nj):
			a = irstar[i, k]
			b = np.interp(a, ir[:,k], E)

			fstar[i, k] = distribution_function(b, 0)

	return fstar


#
# do the calculations on a grid
#

# takes 13 minutes to run with N = 2000
N = 2001
Nj = 2001

# radial grid
rmin = -4.0
rmax = 4.0
r_grid = np.logspace(rmin, rmax, N)
phi_grid = potential(r_grid)
rho_grid = density(r_grid)
ir_grid = calc_radial_action(r_grid, phi_grid)
f_grid = np.zeros((N, Nj))
for i in range(N):
    for j in range(Nj):
        f_grid[i, j] = distribution_function(phi_grid[i], 0)

vr2_grid = calc_velocity_moment(r_grid, rho_grid, phi_grid, f_grid, 2, 0)
vt2_grid = calc_velocity_moment(r_grid, rho_grid, phi_grid, f_grid, 0, 2)
		
mbh = 0.01
error = 5e-2

rhostar_grid = rho_grid
delta = 1e6
while delta > error:
	print("Calculating potential with BH mass", mbh)
	phistar_grid = calc_potential(r_grid, rhostar_grid, mbh)
	print("Calculating radial action ...")
	irstar_grid = calc_radial_action(r_grid, phistar_grid)
	print("Calculating DF ...")
	fstar_grid = calc_dfstar(phi_grid, ir_grid, irstar_grid)
	print("Calculating density ...")
	rhostar_grid_temp = calc_density(r_grid, phistar_grid, fstar_grid)
	delta = np.max(np.abs(rhostar_grid_temp - rhostar_grid))
	rhostar_grid = rhostar_grid_temp
	print(f"Current delta = {delta}")
		

# calculate velocity moments
vr2star_grid = calc_velocity_moment(r_grid, rhostar_grid, phistar_grid, fstar_grid, 2, 0)
vt2star_grid = calc_velocity_moment(r_grid, rhostar_grid, phistar_grid, fstar_grid, 0, 2)

# save it all
np.savetxt(f"pred0.01m_e0.05.txt", np.column_stack((r_grid, rho_grid, phi_grid, vr2_grid, vt2_grid, rhostar_grid, phistar_grid, vr2star_grid, vt2star_grid)))

