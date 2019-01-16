from scipy import signal
import sys
import numpy as np
import random
class IsingLattice(object):
    def __init__(self, N, flip_prop):
        self.J = 1
        self.N = N
        self.NN = N*N
        self.B = 0
        self.flip_prop = flip_prop
        self.neighbors = None
        # self.B = 0
        # if self.n_flip != int(self.n_flip):
            # self.n_flip = int(self.n_flip) + 1

        self.spin = np.random.choice([-1,1],(N,N))
        self.conv_mat = np.matrix('0 1 0; 1 0 1; 0 1 0')

        # self.matrix_ptr = c_void_p(libc.newMatrix(int(N),int(n_flip)))
        self.auto_correlation = []
    def free_memory(self):
        pass
    def step(self, T, B):
        #Calculating the total spin of neighbouring cells
        self.B = B
        # print('B: ',self.B)
        self.neighbors = signal.convolve2d(self.spin,self.conv_mat,mode='same',boundary='wrap')
        M = float(np.sum(self.spin))/float(self.NN)
        E = float(-self.J*(np.sum((self.spin*self.neighbors)))/2.0)/float(self.NN) - float(B)*M
        DeltaE = 2.0 * (self.J*(self.spin*self.neighbors) + float(self.B)*self.spin)

        #Calculate the transition
        p_trans = np.where(DeltaE >= 0.0, np.exp(-1.0*DeltaE/float(T)),1.0)

        #Decide which transitions will occur
        transitions = [[-1 if (cell>random.random() and self.flip_prop>random.random()) 
                        else 1 for cell in row] for row in p_trans]

        #Perform the transitions
        self.spin *= transitions
        return None
    def nsteps(self, T, B, n):
        for i in range(n):
            self.step(T,B)
        return None
        # return libc.nsteps(self.matrix_ptr,n,T,B)
    def get_E(self):
        if type(self.neighbors) == type(None):
            self.neighbors = signal.convolve2d(self.spin,self.conv_mat,mode='same',boundary='wrap')
        return float(-self.J*(np.sum((self.spin*self.neighbors)))/2.0)/float(self.NN) - float(self.B)*self.get_M()
    def get_M(self):
        return float(np.sum(self.spin))/float(self.NN)
    def calc_auto_correlation(self):
        # n = len(spin)
        corr_array = []
        for k in range(1,int(self.N/2)):
            col_mean, row_mean = self.spin.mean(axis=0),self.spin.mean(axis=1)
            #compute r values for rows and cols
            r_col = [np.multiply(self.spin[j,:]-col_mean,self.spin[(j+k)%self.N,:]-col_mean) for j in range(1,self.N)]
            r_row = [np.multiply(self.spin[:,j]-row_mean,self.spin[:,(j+k)%self.N]-row_mean) for j in range(1,self.N)]
            #normalize r values
            r_col = np.divide(r_col,float(self.N))
            r_row = np.divide(r_row,float(self.N))
            #calculate corr for k and add it to array
            corr = (r_col.mean() + r_row.mean())/2.0
            corr_array.append([k,corr])
        return corr_array
    def set_Nflip(self,npick):
        self.flip_prop = npick / self.NN
    def set_flip_prop(self,flip_prop):
        self.flip_prop = flip_prop
    def get_Nspin(self):
        return np.sum(self.spin)
    def get_Nalign(self):
        return np.sum(self.neighbors)/2
        # return libc.get_Nbond(self.matrix_ptr)
    def get_spin(self, i, j):
        return self.spin[i,j]
    def set_spin(self, i, j, k):
        self.spin[i,j] = k
        self.neighbors = None
    def print_spins(self):
        print(self.spin)
    def print_aligned(self):
        if type(self.neighbors) == type(None):
            self.neighbors = signal.convolve2d(self.spin,self.conv_mat,mode='same',boundary='wrap')
        print(self.neighbors)
    def randomize_spins(self):
        self.spin = np.random.choice([-1,1],(self.N,self.N))
    def get_numpy_spin_matrix(self):
        return self.spin

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
        lattice.print_aligned()
        lattice.nsteps(4.,2.1,10000)

        # lattice.randomize_spins()
        print('E: ',lattice.get_E())
        print('M: ',lattice.get_M())
        print('autocorrelation ', lattice.calc_auto_correlation())
        print("set_Nflip:      ", lattice.set_Nflip(20))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        # print("step:           " , lattice.step(2.9,1.0))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        print('E: '              , lattice.get_E())
        print("N               " , lattice.N)
        print('aligned')
        # lattice.print_aligned()
        print('spins')
        # lattice.print_spins()
        lattice.free_memory()
    print('exit')
    sys.exit()
