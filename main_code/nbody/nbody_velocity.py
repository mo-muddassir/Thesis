#!/usr/bin/env python3

from nbody import *
from random import *
import numpy as np
import argparse
from sys import argv
from sys import stderr
from os.path import splitext
import matplotlib.pyplot as plt
from math import log10


parser = argparse.ArgumentParser(description='Calculate the binned density profile of the N-body system.')
parser.add_argument('filename', type=str, help='file name to read in')
parser.add_argument('-n', type=int, dest='num_bins', default=100, help='Number of bins')
parser.add_argument('-f', dest = 'outfile', default = 'nbody-velocity-plot.pdf', help = 'Output file name')
parser.add_argument('-a', dest = 'R', default = 1, help = 'plummer radius scale factor')
parser.add_argument('-equal', dest = 'equal', action='store_true',help = 'use equal mass bins')
parser.add_argument('-e', dest = 'softening', default = 0.0005 ,help = 'softening')
args = parser.parse_args()


#Read system details from source file
file = System.read(args.filename)



def radial_spacing(s):
	rmax = s.rmax + 10**-3
	rmin = s.rmin
	dr = log10(rmax / rmin) / args.num_bins

	print("  r_min = {0:6.3f}, r_max = {1:6.3f}, dr = {2:6.3f}".format(rmin, rmax,dr), file=stderr)


	# Calculate and store mean radial velocities and rms velocities

	mean_vr = np.zeros(args.num_bins)

	mean_vr_2 = np.zeros(args.num_bins)

	bin_pop = np.zeros(args.num_bins)

	bin_vr = np.zeros(args.num_bins)

	bin_vr_2 = np.zeros(args.num_bins)
	
	c = 0

	for p in s: 
	
		c+=1
		r = np.sqrt(np.sum(p['position']**2))

		vr = ((p['position'][0]*p['velocity'][0])+(p['position'][1]*p['velocity'][1])+(p['position'][2]*p['velocity'][2]))/r
	
		pos = int( (log10(r) - log10(rmin)) / dr)
	
		bin_pop[pos] += 1
	
		bin_vr[pos] += vr
	
		bin_vr_2[pos] += vr**2
	
	radii = np.zeros(args.num_bins)

	for i in range(args.num_bins):
		radii[i] = np.power(10.0, np.log10(rmin) + 0.5 * dr + dr*i)
	
		mean_vr[i] = bin_vr[i] / bin_pop[i]
	
		mean_vr_2[i] = bin_vr_2[i] / bin_pop[i]
		
	rows_to_delete = []
	for i in range(len(radii)):
		if radii[i]<float(args.softening):
			rows_to_delete.append(i)
	radii = np.delete(radii, rows_to_delete)
	
	mean_vr_2 = np.delete(mean_vr_2, rows_to_delete)
	

	return radii,mean_vr_2

def equal_mass(s):
	rmax = s.rmax + 1e-3
	rmin = s.rmin
	dr = log10(rmax / rmin) / args.num_bins

	print("  r_min = {0:6.3f}, r_max = {1:6.3f}, dr = {2:6.3f}".format(rmin, rmax,dr), file=stderr)
	
	print("Binning in radial shells with equal mass spacing", file=stderr)
    
	num_per_bin = len(s) / args.num_bins + 1
	s.sort_by_radius()
    
	rmax = s.rmax + 1e-10
	rmin = s.rmin
	print(f"  r_min = {rmin:6.3f}, r_max = {rmax:6.3f}", file=stderr)
    
	m_enc = 0.0
	r_last = rmin
	c = 0
	vr_2 = 0 
	bin_vr_2 = np.zeros(args.num_bins)
	radii = np.zeros(args.num_bins)
	density = np.zeros(args.num_bins)
	mean_vr_2 = np.zeros(args.num_bins)
    
	for i in range(len(s)):
		m_enc += s[i]['mass']
		r = np.sqrt(np.sum(s[i]['position']**2))

		vr_2 += np.power((((s[i]['position'][0]*s[i]['velocity'][0])+(s[i]['position'][1]*s[i]['velocity'][1])+(s[i]['position'][2]*s[i]['velocity'][2]))/r),2)
        
		if (i != 0 and i % num_per_bin == 0) or (i == len(s) - 1):        
        
			r = s.radius(i)
			radii[c] = r_last + 0.5 * (r - r_last)
			bin_vr_2[c] = vr_2
			r_last = r
			c += 1
			m_enc = 0.0
			vr_2 = 0 
			
	for i in range(args.num_bins):
	
		mean_vr_2[i] = bin_vr_2[i] / num_per_bin
	
	return radii, mean_vr_2
	


if args.equal:
	print("With equal mass bins")
	radii, mean_vr_2 = equal_mass(file)
else:
	print("Binning in radial shells with logarithmic spacing", file=stderr)
	radii, mean_vr_2= radial_spacing(file)
	
	
np.savetxt(splitext(args.filename)[0] + ".vel", np.transpose([np.log10(radii), mean_vr_2])) 


