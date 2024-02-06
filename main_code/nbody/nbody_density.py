#!/usr/bin/env python3

from nbody import *
from random import *
import numpy as np
import argparse
from sys import argv
from sys import stderr
from os.path import splitext
from math import log10

parser = argparse.ArgumentParser(description='Calculate the binned density profile of the N-body system.')
parser.add_argument('filename', type=str, help='file name to read in')
parser.add_argument('-f', dest = 'outfile', default = 'Density_profile',required = True, help = 'output file name base')
parser.add_argument('-n', type=int, dest='num_bins', default=100, help='Number of bins')
parser.add_argument('-e', type=float, dest='softening', default=None, help='Softening length:  drop data less than this radius')
parser.add_argument('-c', action='store_true', dest='cm', help='Adjust for centre of mass')
parser.add_argument('-equal', action='store_true', dest='equal', help='Use equal-mass bins')

args = parser.parse_args()

s = System.read(args.filename)

if args.cm:
    s.translate_to(s.centre_of_mass())

if args.equal:

    print("Binning in radial shells with equal mass spacing", file=stderr)
    
    num_per_bin = len(s) / args.num_bins + 1
    s.sort_by_radius()
    
    rmax = s.rmax + 1e-10
    rmin = s.rmin
    print(f"  r_min = {rmin:6.3f}, r_max = {rmax:6.3f}", file=stderr)
    
    m_enc = 0.0
    r_last = rmin
    c = 0
    
    radii = np.zeros(args.num_bins)
    density = np.zeros(args.num_bins)
    
    for i in range(len(s)):
        m_enc += s[i]['mass']
        
        if (i != 0 and i % num_per_bin == 0) or (i == len(s) - 1):        
        
            r = s.radius(i)
            radii[c] = r_last + 0.5 * (r - r_last)
            volume = 4.0/3.0 * np.pi * (r**3 - r_last**3)
            density[c] = m_enc / volume
            r_last = r
            c += 1
            m_enc = 0.0   
    

else:

    print("Binning in radial shells with logarithmic spacing", file=stderr)

    rmax = s.rmax + 1e-3
    rmin = s.rmin
    dr = log10(rmax / rmin) / args.num_bins
    print("  r_min = {0:6.3f}, r_max = {1:6.3f}".format(rmin, rmax), file=stderr)

    density = np.zeros(args.num_bins)
    for p in s:
        r = np.sqrt(np.sum(p['position']**2))
        pos = int( (log10(r) - log10(rmin)) / dr)
        
        density[pos] += p['mass']

    radii = np.zeros(args.num_bins)
    for i in range(args.num_bins):
        radii[i] = np.power(10.0, log10(rmin) + 0.5 * dr + dr*i)
        r = np.power(10.0, dr*i) * rmin
        rp1 = np.power(10.0, dr) * r
        volume = 4.0 / 3.0 * np.pi * ( np.power(rp1, 3.0) - np.power(r, 3.0) )
        
        density[i] /= volume

    rows_to_delete = []
    for i in range(len(radii)):
    	if radii[i]<0.0005:
    		rows_to_delete.append(i)
    #for i in range(len(radii)):
        #if density[i] == 0 or (args.softening != None and radii[i] < 0.1):
            #rows_to_delete.append(i)
     
            

            
    radii = np.delete(radii, rows_to_delete)
    density = np.delete(density, rows_to_delete)

radii = np.log10(radii)
density = np.log10(density)
    
np.savetxt(f"/home/moe/research/Thesis/analysis/{args.outfile}.den", np.transpose([radii, density]))  

