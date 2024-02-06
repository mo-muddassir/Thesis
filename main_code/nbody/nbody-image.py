#!/usr/bin/env python3

from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import mpl_toolkits.mplot3d.axes3d
import matplotlib.animation as animation
import matplotlib as mpl
import argparse
from sys import argv
from time import sleep
import scipy.spatial

parser = argparse.ArgumentParser(description='Create an image of an N-body file.')
parser.add_argument('filename', type=str, help='the list filename of the snapshot')
parser.add_argument('-save', dest='out_file', default=None, help='Save as a file')
parser.add_argument('-G', dest="gadget", action='store_true', help="Use gadget file format")
parser.add_argument('-O', dest="old", action='store_true', help="Use old file format")
parser.add_argument('-pixels', dest='pixels', default=256, type=int, help="Number of pixels")
parser.add_argument('-density', dest='density', action='store_true', help="Use density of particles (slow)")
args = parser.parse_args()

if args.gadget:
    s = System.read_gadget(args.filename)
elif args.old:
    s = System.read_old(args.filename)
else:
    s = System.read(args.filename)
    
xmax = s.xmax
xmin = s.xmin
delta = (xmax-xmin) / args.pixels

grid = np.zeros((args.pixels, args.pixels))

# should split this out into a separate file and program
weights = np.ones(len(s))
if args.density:

    print("----building KD tree ...")
    tree = scipy.spatial.cKDTree(s['position'])
    print("    ... done")
    
    max_parts = 3000000
    if len(s) > max_parts:  # need to split it up for memory reasons
        pos = np.array_split(s['position'], np.ceil(len(s)/max_parts))
        
        ntot = 0
        for j in range(len(pos)):
            print("QUERYING", j, " of", len(pos), ": ", len(pos[j])," particles")
            d,i = tree.query(pos[j], k=64)
            mtot = np.sum(s['mass'][i], axis=1)
            vol = 4.0 / 3.0 * np.pi * d[:, -1]**3
            den = mtot/vol
            
            weights[ntot:ntot + len(pos[j])] = den*den
            ntot += len(pos[j])
    else:
        pos = s['position']
        d, i = tree.query(pos, k = 64)

        mtot = np.sum(s['mass'][i], axis=1)
        vol = 4.0 / 3.0 * np.pi * d[:, -1]**3
        den = mtot/vol
        weights = den*den
    

h, x, y = np.histogram2d(s.all_x(), s.all_y(), bins=args.pixels, weights=weights)
#h2, x, y = np.histogram2d(s.all_x(), s.all_y(), bins=args.pixels)
#h2 += 0.00001
#h1 += 0.00001
#h = h1/h2
#print(h)

h+=0.00001
h = np.log(h)

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
plt.imshow(h.T, origin='lower', extent=[x[0],x[-1],y[0],y[-1]], interpolation='bicubic', cmap=cm.magma)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.colorbar()

#ax = fig.add_subplot(111, aspect='equal')
#X, Y = np.meshgrid(x, y)
#ax.pcolormesh(X, Y, h, cmap=cm.afmhot)

#ax = fig.add_subplot(111, title='NonUniformImage: interpolated', aspect='equal', xlim=x[[0, -1]], ylim=y[[0, -1]])
#im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
#xcenters = (x[:-1] + x[1:]) / 2
#ycenters = (y[:-1] + y[1:]) / 2
#im.set_data(xcenters, ycenters, h)
#ax.images.append(im)

if args.out_file == None:
    plt.show()
else:
    plt.savefig(args.out_file)
