#!/usr/bin/env python3

import os

import argparse
from os.path import splitext

f = '/home/moe/research/Thesis/init_files/1M_init.dat'
outfile = '1M_den'
os.system(f"/home/moe/research/Thesis/main_code/nbody/nbody_density.py {f} -f {outfile}")
