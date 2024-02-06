from struct import unpack, pack
from array import array
from sys import stderr
from random import random, seed
import time
import copy

import numpy as np

"""
An N-body system
"""

class System:

    _particle_type = np.dtype([('mass', 'f4'), ('position', '3f4'), ('velocity', '3f4'), ('id', 'i4')])
    _particle_type_dnc = np.dtype([('mass', 'f4'), ('position', '3f4'), ('velocity', '3f4')])
    _header_type = np.dtype([('N', 'i4'), ('time', 'f8')])


#
# N-body system creation
#


    def __init__(self, time = 0.0, particles = [], G = 1.0):
        
        self.time = time
        self._particles = particles
        
        self.T = 0.0
        self.U = 0.0
        
        self.G = G
        
    @classmethod
    def cube_grid(cls, Ng, L, Mtot, time = 0.0):
        N = Ng**3
        dx = L / Ng
        m = Mtot / N
        
        print("  Making grid: {0} particles per side;  side length {1}; total mass {2}".format(Ng, L, Mtot), file=stderr)
        
        parts = np.zeros(N, dtype=cls._particle_type)
        parts['mass'] = m
        x, y, z = np.mgrid[0:L:dx, 0:L:dx, 0:L:dx]
        parts['position'] = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
        
        #
        # Have not implemented perturbing the positions by some displacements Sx, Sy, Sz
        #
        
        s = cls(time = time, particles = parts)
        
        return s
        
    @classmethod
    def cube_random(cls, N, L, Mtot, t = 0.0, random_seed = None):
        print("  Making random cube: {0} particles;  side length {1}; total mass {2}".format(N, L, Mtot), file=stderr)
    
        seed(random_seed)
        m = Mtot / N
        
        parts = np.zeros(N, dtype=cls._particle_type)
        
        parts['mass'] = m
        parts['position'] = L * np.random.rand(N, 3)            

        return cls(time = t, particles = parts)

    @classmethod
    def sphere_random(cls, N, R, Mtot, t = 0.0, random_seed = None):
        print("  Making random sphere: {0} particles;  radius {1}; total mass {2}".format(N, R, Mtot), file=stderr)
        seed(random_seed)
        m = Mtot / N
        
        parts = np.zeros(N, dtype=cls._particle_type)
        
        parts['mass'] = m
        
        phi = np.random.uniform(0, 2.0*np.pi, N)
        cosphi = np.cos(phi)
        sinphi = np.sin(phi)
        costheta = np.random.uniform(-1.0, 1.0, N)
        sintheta = np.sin(np.arccos(costheta))
        r = R * np.cbrt(np.random.uniform(0, 1, N))
        
        x = r * sintheta * cosphi
        y = r * sintheta * sinphi
        z = r * costheta
        
        p = np.transpose(np.array([x, y, z]))
        #print(p)
        
        parts['position'] = p
        
        return cls(time = t, particles = parts)
   
    def add_particles(self, p):
        p_new = np.append(self._particles, p) 
        self._particles = p_new  # slow!
        
    def remove_particle(self, pi):
        np.delete(self._particles, pi)   # slow!
        
    def __getitem__(self, index):
        return self._particles[index]
        
    def __len__(self):
        return len(self._particles)
    
    def __str__(self):
        return "N-body system: {0} particles at time {1}".format(len(self), self.time)
        
#
# Reading/writing routines
#        

        
    def write(self, filename):
        print("Writing nbody system to file {0} (as binary)".format(filename), file=stderr)
        f = open(filename, 'wb')
        
        h = np.array([(len(self._particles), self.time)], dtype=self._header_type)
        h.tofile(f)
        
        self._particles.tofile(f)
        
        print("--Wrote {0} particles at time {1}".format(len(self._particles), h['time'][0]), file=stderr)
        
        f.close()
        
        
    @classmethod
    def read(cls, filename):
        print("Reading nbody system from file {0}".format(filename), file=stderr)
        f = open(filename, 'rb')
        h = np.fromfile(f, count=1, dtype=cls._header_type)
        parts = np.fromfile(f, dtype=cls._particle_type)
        f.close()
        
        print("--Read {0} particles at time {1}".format(len(parts), h['time'][0]), file=stderr)

        s = cls(time = h['time'][0], particles = parts)
        
        return s
        
    @classmethod
    def read_dnc(cls, filename):
        print("Reading nbody system from file {0}".format(filename), file=stderr)
        f = open(filename, 'rb')
        h = np.fromfile(f, count=1, dtype=cls._header_type)
        parts = np.fromfile(f, dtype=cls._particle_type_dnc)
        f.close()
        
        print("--Read {0} particles at time {1}".format(len(parts), h['time'][0]), file=stderr)

        s = cls(time = h['time'][0], particles = parts)
        
        return s
        
    @classmethod
    def read_gadget(cls, filename):
    
        # this function will read in a Gadget file and return a system object.
        # However:  it assumes there is only one file and that the mass is specified in mass_arr.
    
        f = open(filename, 'rb')
        print("Reading Gadget N-body system from {0}".format(filename), file=stderr)
        block_size = unpack('<i', f.read(4))[0]
        
        # read in particle numbers; there are six types.
        npart = [0,0,0,0,0,0]
        N = 0
        for i in range(6):
            npart[i] = unpack('<I', f.read(4))[0]
            if npart[i] > 0:
                print("--Number of particle of type {0}: {1}".format(i, npart[i]), file=stderr)
            N += npart[i]
            
        # read in masses for all 6 types
        mass_arr = [0,0,0,0,0,0]
        variable_mass_arr = [0,0,0,0,0,0]
        for i in range(6):
            mass_arr[i] = unpack('<d', f.read(8))[0]
            if mass_arr[i] > 0.0:
                print("--Mass of particle of type {0}: {1}".format(i, mass_arr[i]), file=stderr)
            if mass_arr[i] == 0.0 and npart[i] > 0:
                print("----Particle type {0} has no mass specified, will read it individually".format(i), file=stderr)
                variable_mass_arr[i] = npart[i]
            
        # I'll just get time; the rest can come later if we need it
        t = unpack('<d', f.read(8))[0]
        print("--Time:  {0}".format(t))        
        
        
        f.seek(block_size + 8)        
        
        print("--Reading {0} total positions and velocities".format(N), file=stderr)        
        
        # now read in particle positions
        block_size = unpack('<i', f.read(4))[0]        
        
        pos_arr = array('f')
        pos_arr.fromfile(f, N * 3)
        block_size = unpack('<i', f.read(4))[0]
        
        # now velocities ...
        block_size = unpack('<i', f.read(4))[0]
        vel_arr = array('f')
        vel_arr.fromfile(f, N * 3)
        block_size = unpack('<i', f.read(4))[0]
        
        # now particles ids ...
        block_size = unpack('<i', f.read(4))[0]
        id_arr = array('I')
        id_arr.fromfile(f, N)
        block_size = unpack('<i', f.read(4))[0]
                        
        # now masses if necessary
        mass_arr_indiv = array('f')
        tot_masses = np.sum(variable_mass_arr)
        if tot_masses > 0:
            print("  Reading {0} masses for individual particles".format(tot_masses), file=stderr)       
            block_size = unpack('<i', f.read(4))[0]
            mass_arr_indiv.fromfile(f, tot_masses)
            block_size = unpack('<i', f.read(4))[0]
        
        f.close()
        
        parts = np.zeros(N, dtype=cls._particle_type)
        c = 0
        c_mass = 0
        for i in range(6):
            for j in range(npart[i]):
                if mass_arr[i] > 0.0:
                    parts[c]['mass'] = mass_arr[i]
                else:
                    parts[c]['mass'] = mass_arr_indiv[c_mass]
                    c_mass += 1
                parts[c]['position'] = [pos_arr[3 * c + 0], pos_arr[3 * c + 1], pos_arr[3 * c + 2]]
                parts[c]['velocity'] = [vel_arr[3 * c + 0], vel_arr[3 * c + 1], vel_arr[3 * c + 2]]
                parts[c]['id'] = id_arr[c]
                c += 1

        s = cls(time = t, particles = parts)
        return s
    
    
    def write_gadget(self, filename, redshift=0.0, box_size=0.0, omega0=0.0, omega_lambda=0.0, hubble_param=0.0):
        print("Writing nbody system to file {0} (as Gadget format)".format(filename), file=stderr)
        
        # this function will make the folowing assumptions:
        # * only dark matter particles
        # * all particles have the same mass
        # * only one file
        
        f = open(filename, 'wb')
        
        # header start
        
        block_size = pack('<i', 256)
        f.write(block_size)
        
        Npart = [0, len(self), 0, 0, 0, 0]
        f.write(pack('<6i', *Npart))
        mass_arr = [0.0, self[0]['mass'], 0.0, 0.0, 0.0, 0.0]
        f.write(pack('<6d', *mass_arr))
        
        f.write(pack('<d', self.time))
        f.write(pack('<d', redshift))  # redshift
        f.write(pack('<i', 0))    # FlagSfr
        f.write(pack('<i', 0))    # FlagFeedback
        f.write(pack('<6i', *Npart))   # Nall[6]
        f.write(pack('<i', 0))    # FlagCooling
        f.write(pack('<i', 1))    # NumFiles
        f.write(pack('<d', box_size))  # BoxSize
        f.write(pack('<d', omega0))  # Omega0
        f.write(pack('<d', omega_lambda))  # OmegaLambda
        f.write(pack('<d', hubble_param))  # HubbleParam
        #f.write(pack('<i', 0))    # FlagAge
        #f.write(pack('<i', 0))    # FlagMetals
        #f.write(pack('<6i', 0, 0, 0, 0, 0, 0))    # NallHW[6]
        #f.write(pack('<i', 0))    # flag_entr_ics
        f.write(pack('<96x'))
        
        f.write(block_size)
        
        # header done!
        
        # next is particle positions
        block_size = pack('<i', 3 * 4 * len(self))
        f.write(block_size)
        
        pos_arr = array('f')
        for p in self._particles:
            pos_arr.fromlist([p['position'][0], p['position'][1], p['position'][2] ])
        pos_arr.tofile(f)
        
        f.write(block_size)
        
        # next is velocities
        f.write(block_size)
        
        vel_arr = array('f')
        for p in self._particles:
            vel_arr.fromlist([p['velocity'][0], p['velocity'][1], p['velocity'][2]  ])
        vel_arr.tofile(f)
        
        f.write(block_size)
        
        # particle IDs
        block_size = pack('<i', 4 * len(self))
        f.write(block_size)
        id_arr = array('I')
        id_arr.fromlist([ i for i in range(len(self)) ])
        id_arr.tofile(f)
        f.write(block_size)
        
        f.close()

#
# Manipulating and info routines
#
    def set_particle_mass(self, m):
        self._particles['mass'] = m
        
    def scale_positions(self, scale_factor):
        self._particles['position'] *= scale_factor
        
    def scale_velocities(self, scale_factor):
        self._particles['velocity'] *= scale_factor
      
    def all_x(self):
        return self._particles['position'][:,0]
        
    def all_y(self):
        return self._particles['position'][:,1]   
        
    def all_z(self):
        return self._particles['position'][:,2]
    
    @property
    def mass(self):
        return np.sum(self._particles['mass'])
    
    @property
    def rmax(self):
        return np.sqrt(np.amax(np.sum(self._particles['position']**2, axis=1)))
        
    @property
    def rmin(self):
        return np.sqrt(np.amin(np.sum(self._particles['position']**2, axis=1)))
        
    @property
    def xmax(self):
        return np.amax(self._particles['position'][:,0])
        
    @property
    def xmin(self):
        return np.amin(self._particles['position'][:,0])
        
    @property
    def ymax(self):
        return np.amax(self._particles['position'][:,1])
        
    @property
    def ymin(self):
        return np.amin(self._particles['position'][:,1])
        
    @property
    def zmax(self):
        return np.amax(self._particles['position'][:,2])
        
    @property
    def zmin(self):
        return np.amin(self._particles['position'][:,2])

    def translate_to(self, new_position):
        print("--Translating system to new origin {0}".format(new_position), file=stderr) 
        
        self._particles['position'] -= new_position
        
    def extract_sphere(self, radius):
        r2 = radius**2
        
        radii2 = np.sum(self._particles['position']**2, axis=1)
        
        print("--Extracting sphere of radius {0}".format(radius), file=stderr) 
        parts = np.extract(radii2 < r2, self._particles)
        
        return System(self.time, parts)
        
    def extract_cube(self, L):
        # extent of cube is -L/2 to L/2
        
        pos = self['position']
        inn = (pos[:,0] > -0.5 * L) & (pos[:,0] < 0.5 * L) & (pos[:,1] > -0.5 * L) & (pos[:,1] < 0.5 * L) & (pos[:,2] > -0.5 * L) & (pos[:,2] < 0.5 * L)
        #print(inn)
        
        print("--Extracting cube of length {0}".format(L), file=stderr) 
        parts = np.extract(inn, self._particles)
        
        return System(self.time, parts)
        
    def _find_global_centre_of_mass(self):
        
        # here's the global c/m
        
        mtot = self.mass
        xm = self._particles['mass'] * self._particles['position'][:,0]
        ym = self._particles['mass'] * self._particles['position'][:,1]
        zm = self._particles['mass'] * self._particles['position'][:,2]
        
        xcm = np.sum(xm) / mtot
        ycm = np.sum(ym) / mtot
        zcm = np.sum(zm) / mtot
        
        return np.array([xcm, ycm, zcm])
    
    def centre_of_mass(self, tolerance = None, cm = None):
        
        if cm == None:
            cm = self._find_global_centre_of_mass()
            
        if tolerance == None:
            return cm
        
        print("----Finding the centre of mass iteratively, hold on ...", file=stderr) 
        s = copy.deepcopy(self)
                
        error = 100.0
        cm_static = cm
        while error > tolerance:
            s.translate_to(cm)
            s = s.extract_sphere(0.75 * s.rmax)
            cm_new = s._find_global_centre_of_mass()
            error = np.sqrt(np.sum((cm - cm_new)**2)) / np.sqrt(np.sum((cm_static)**2))
            
            cm = cm_new
            cm_static += cm_new
       
        return cm_static
        
    def radius2(self, i):
        p = self[i]
        return np.sum(np.square(p['position']))
        
    def radius(self, i):
        p = self[i]
        return np.sqrt(np.sum(np.square(p['position'])))
 
    def radii(self):
        return np.sqrt(np.sum(self._particles['position']**2, axis=1))
        
    def kinetic_energy(self):
        ke = 0.5 * self[0]['mass'] * np.sum(np.sum(self['velocity']**2, axis=1))
        return ke
        
    def potential_energy(self):
        pe = 0.0
        for i in range(len(self)):
            print("-- ", i)
            for j in range(i+1, len(self)):
                x = self[i]['position'][0] - self[j]['position'][0]
                y = self[i]['position'][1] - self[j]['position'][1]
                z = self[i]['position'][2] - self[j]['position'][2]
        
                dr = np.sqrt(x**2 + y**2 + z**2)
        
                pe += -self[i]['mass'] * self[j]['mass'] / dr
        
        return pe
        
    def potential_energy_shell(self):
        pe = 0.0
        out_shell = 0.0
        r = self.radii()
        ind = np.argsort(r)
        for i in range(len(self)):    
            out_shell += self[ind[i]]['mass'] / r[ind[i]]

        enc_mass = 0.0        
        for i in range(len(self)):
            enc_mass += self[ind[i]]['mass']
            out_shell -= self[ind[i]]['mass'] / r[ind[i]]
            pe -= (enc_mass / r[ind[i]] + out_shell)
        pe *= 0.5 * self[0]['mass'] 
        return pe

    # returns a list of AM for all particles
    def angular_momenta(self):
        J = np.cross(self['position'], self['velocity'])
        return J
        
    # returns a list of the magnitude of AM for all particles
    def angular_momenta_mag(self):
        J = np.cross(self['position'], self['velocity'])
        J = np.sum(J**2, axis=1)
        return J
    
    # total AM    
    def angular_momentum(self):
        return np.sum(self.angular_momenta_mag())
    
    def sort_by_radius(self):
        radii = np.sum(np.square(self._particles['position']), axis=1)
        key = np.argsort(radii)
        
        self._particles = self._particles[key]
                
    
if __name__ == '__main__':
    
    
    s = System.read("plum.dat")
    print(s.radius(0))
    s.sort_by_radius()
    print(s.radius(0))    
    
