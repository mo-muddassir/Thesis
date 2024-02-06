#!/usr/bin/env python3

import argparse
import numpy as np
from nbody import System
from random import random, seed
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = "Create a Plummer System of N particles of equal mass")

parser.add_argument('-N', dest='N', default=100000, type=int, help='Number of particles')

parser.add_argument('-R', dest = 'plum_rad', default = 1, type = int, help ='Plummer radius')

parser.add_argument('-G', dest='G', default=1, type=int, help='Gravitational Constant')

parser.add_argument('-f', dest='outfile', default=None, help='File to output to')

args = parser.parse_args()

def rand(N):
	return np.random.uniform(0,1,N)


m = 1/args.N
pi=np.pi
i=0
q =[]

X1 = rand(args.N)
X2 = rand(args.N)
X3 = rand(args.N)
X6 = rand(args.N)
X7 = rand(args.N)

r = (pow(X1,-2/3)-1)**(-1/2)

M = (r**3)*pow((1+(r**2)),-3/2)


# Set positions for each particle
	
z = (1-2*X2)*r

x = pow((r**2 - z**2),0.5)*((np.cos(2*pi*X3)))

y = pow((r**2 - z**2),0.5)*((np.sin(2*pi*X3)))



def prop_check(X4,X5):
	g = (X4**2)*((1-X4**2)**(7/2))
	
	if 0.1*X5 < g:
		return True
	else:
		return False
		
while i<args.N:
	var = np.random.uniform(0,0.99)
	var2 = np.random.uniform(0,0.99)
	if prop_check(var,var2) == True:
		q.append(var)
		i+=1

g = np.zeros(args.N)
for i in range(args.N):
	g[i] = (q[i]**2)*((1-q[i]**2)**(7/2))

#escape velocity
Ve = (2**0.5)*pow((1+r**2),-1/4)

#individual total velocity of each particle
V = q*Ve


#Calculate and set velocity components to each particle
w = (1-2*X6)*V

u = pow(abs(V**2 - w**2),0.5)*((np.cos(2*pi*X7)))

v = pow(abs(V**2 - w**2),0.5)*((np.sin(2*pi*X7)))	



#Create system with N particles, mass 1/N, positions [x,y,z], and velocity[u,v,w]	
parts = np.zeros(args.N, dtype=System._particle_type)
parts['mass'] = m

parts['position'] = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
parts['velocity'] = np.vstack((u.flatten(), v.flatten(), w.flatten())).T

plummer_system = System(particles = parts)

plummer_system.write(f"/home/moe/research/Thesis/init_files/{args.N}.dat")
