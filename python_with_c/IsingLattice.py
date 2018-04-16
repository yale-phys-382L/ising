import numpy as np
from ctypes import *
libc = CDLL("./ising_lattice_lib.so")

libc.newMatrix.restype  = c_void_p
libc.newMatrix.argtypes = [c_int, c_int]

libc.step.argtypes        = [ c_void_p, c_float, c_float ]
libc.nsteps.argtypes      = [ c_void_p, c_int, c_float, c_float ]

libc.get_E.restype = c_float
libc.get_M.restype = c_float
libc.auto_correlation.restype = c_float

libc.auto_correlation.restype = c_float

class IsingLattice(object):
    def __init__(self, N, flip_prop):
        n_flip = N*N*flip_prop
        self.N = N
        if n_flip != int(n_flip):
            n_flip = int(n_flip) + 1
        self.matrix_ptr = c_void_p(libc.newMatrix(int(N),int(n_flip)))
        self.auto_correlation = []
    def free_memory(self):
        # print ('deleting ising_lattice_lib object')
        libc.delMatrix( self.matrix_ptr ) 
    def step(self, T, B):
        return libc.step(self.matrix_ptr, c_float(T), c_float(B))
    def nsteps(self, T, B, n):
        return libc.nsteps(self.matrix_ptr,n,T,B)
    def get_E(self):
        return float( libc.get_E(self.matrix_ptr) )
    def get_M(self):
        return float( libc.get_M(self.matrix_ptr) )
    def calc_auto_correlation(self):
        libc.calc_auto_correlation(self.matrix_ptr)
        for i in range(1,int(self.N/2)):
            val = libc.auto_correlation(self.matrix_ptr, i)
            self.auto_correlation.append([i,val])
        return self.auto_correlation
    def set_Nflip(self,npick):
        return  libc.set_Npick(self.matrix_ptr, c_int(npick))
    def set_flip_prop(self,flip_prop):
        n_flip = self.N*flip_prop
        if n_flip != int(n_flip):
            n_flip = int(n_flip) + 1
        else:
            n_flip = int(n_flip)
        return (libc.set_Npick(self.matrix_ptr, c_int(n_flip)))
    def get_Nspin(self):
        return libc.get_Nspin(self.matrix_ptr)
    def get_Nalign(self):
        return libc.get_Nbond(self.matrix_ptr)
    def get_spin(self, i, j):
        return libc.get_spin(self.matrix_ptr, i, j)
    def set_spin(self, i, j, k):
        return libc.set_spin(self.matrix_ptr, i, j, k)
    def print_spins(self):
        return libc.print_spins(self.matrix_ptr)
    def print_aligned(self):
        return libc.print_bonds(self.matrix_ptr)
    def randomize_spins(self):
        return libc.initialize_spins(self.matrix_ptr)
    def get_numpy_spin_matrix(self):
        matrix = np.ones((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                matrix[i,j] = libc.get_spin(self.matrix_ptr,i,j)
        return matrix

if __name__ == '__main__':
    from sys import argv
    n = 10
    if len(argv) == 1:
        pass
    elif argv[1] == "1":
        lattice = IsingLattice(n,.1)
        lattice.print_spins()
        print("----")
        for i in range(n):
            for j in range(n):
                continue
                lattice.set_spin(i,j,1)
        lattice.print_spins()
        print("Nspin : %i"%lattice.get_Nspin())
        print(" This is mag: %f"%lattice.get_M())
        print("auto_corr:")
        x = lattice.calc_auto_correlation()
        for val in x:
            print(val)
        lattice.free_memory()
        
    elif argv[1] == "0":
        lattice = IsingLattice(5,.1)
        print(argv[0])
        lattice.print_spins()
        print('---')
        lattice.set_spin(0,0,1)
        lattice.print_spins()
        print('---')
        lattice.set_spin(0,0,-1)
        lattice.print_spins()
        print('---')

        flip = True
        for i in range(n):
            for j in range(n):
                flip = not flip
                if flip:
                    lattice.set_spin(i,j,1)
                else:
                    lattice.set_spin(i,j,-1)
        lattice.print_spins()
        print('---')

        lattice.free_memory()
    else:
        lattice = IsingLattice(10,.01)
        lattice.print_spins()
        print('---')
        lattice.print_aligned()
        lattice.nsteps(8.,0.1,10000)
        print('---')
        lattice.print_spins()
        print('---')
        lattice.print_aligned()

        lattice.randomize_spins()
        # lattice.nsteps(100,2.2,3.0)
        print('E: ',lattice.get_E())
        print('M: ',lattice.get_M())
        print('autocorrelation ', lattice.calc_auto_correlation())
        print("set_Nflip:      ", lattice.set_Nflip(20))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        print("nsteps:         " , lattice.nsteps(2.9,1.0,10000))
        print("set_flip_prop   " , lattice.set_flip_prop(0.20))
        print("nsteps:         " , lattice.nsteps(2.9,1.0,20000))
        print("step:         " , lattice.step(2.9,1.0))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        print('E: ',lattice.get_E())
        print("N               " , lattice.N)
        lattice.print_aligned()
        lattice.free_memory()
    exit(1)
