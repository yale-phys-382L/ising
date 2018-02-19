# Ising Python with C implementation

## Differences from normal python version
There are two new files:
- `ising_lattice_lib.cxx`: This is a library written in C++ and C
- `IsingLattice.py`: This defines the Python class `IsingLattice` that is basically wrapper to use the `ising_lattice_lib` library

The files `main_c.py`, `ising_c.py`, and `main-multiprocessing_c.py` are minimally modified from the files `main.py`, `ising.py`, and `main-multiprocessing.py` in the ising-starter-python respository in order to use the C++ library.

## Why use this?
Speed: The C++/C library implements (with one exception, see below) the same algorithm for the lattice simulation as in the pure python code. The code base is much more verbose, but it is also faster (in informal testing somewhere from 7 to more than 20 times as fast, depending on the lattice size of the simulation).

## Why not use this?
1. In order to use this you have to compile the `ising_lattice_lib.cxx` file into a shared library.
2. The code works by generating `IsingLattice` python objects which make calls to the C++ library. Exception handling between the library and Python is not coded. If something fails the resulting Python message probably won't be very informative. (Don't be too scared by this, the C++ library has been tested, and shouldn't fail...)
3. You must remember to release the memory conumed by each instance of IsingLattice that you generate in Python, otherwise you will have a memory leak.
4. The results of `main-multiprocessing_c.py` are suspiciously fast. (This is not meant facetiously, it is *much* faster than anticipated). While brief debugging hasn't uncovered any errors, and it does appear to generate full results, take it with a grain of salt, or better yet, email either a bug report or a justification for its performance.

## How to use this?
1. Compile ising_lattice_lib.cxx into a shared library. On OSX the command is: `g++ --std=c++11 -shared -o ising_lattice_lib.so ising_lattice_lib.cxx`. 
2. Use the python object `IsingLattice`, which is basically a wrapper for calls to the shared library. This is already done in `main_c.py`, `ising_c.py` and `main-multiprocessing_c.py`. (See #4 above under "Why not use this?" for a disclainer on use of the `main-multiprocessing_c.py`)

Generally, use the library though the `IsingLattice` class in Python. For example:
  `lattice = IsingLattice(<N>,<flip_prop>)` will generate a class instance `lattice`. This actually talks to the C-library which then generates a pointer to a C++ class which generates the Ising Lattice. C++ was used because it has better support for a thread safe random number generator. C was used because it works with the Python module ctypes.
  
  Various member functions of `IsingLattice` call C++ to perform operations and return results. Python has a garbage collector that deletes objects as soon as there are no longer any pointers to them. However, it will note automatically call the C++ library to delete an object to which it is making calls. Therefore, when you are done with any instance of `IsingLattice`, call the `.free_memory()` function. (Obviously, this is the *last* call to be made any given instance of `IsingLattice`).
  
  ## Some useful function calls to IsingMatrix
  
  - `__init__(self, N, flip_prop)`: initializes the ising lattice object in C++
  - `free_memory(self)`: frees the memory in the C++ library. If this isn't called, then when the `IsingLattice` python object is destroyed, the object in the C++ library no longer is pointed to by anything and becomes leaked memry.
  - `step(self, T, B)`: run on step of the ising simulation
  - `nsteps(self, T, B, n)`: run n-steps of the ising simulation
  - `get_E(self)`: return the energy per lattice site value for the current lattice configuation
  - `get_M(self)`: return the magnetization for the current lattice configuration
  - `set_flip_prop`(self, flip_prop): reset the flip_prop. Useful during a simulated annealing process.
  - `calc_auto_correlation(self)`: returns a Python list of the autocorrelation for the current lattice configuration.
  - `randomize_spins(self)`: re-start the lattice to a randomized spin state.
  
  It is likely that those are the only functions which will be used. For convenience,
  there are the following (and a few more besides, see `IsingLattice.py` for full details.
  - `print_spings(self)`: print out a matrix of 0's and 1's for the lattice site (1's are aligned with the magnetic field and 0's anti-aligned), example:  
     `1 0 0 0 0`  
     `0 0 1 0 1`  
     `1 0 0 0 0`    
     `1 1 0 0 1`      
     `0 1 0 1 1`  
 - `print_aligned(self)`: print out a matrix of how many parallel aligned neighbors each lattice site has (+1 for aligned, -1 for anti-alligned). (Looks like above, but with -4, -2, 0, 2, and 4 entries).
 - `get_numpy_spin_matrix()`: returns a numpy matrix of the spins (like the one above, but with -1s and 1s, instead of 1's and 0's.
 
## Difference in algorithms

The one place where the algorithm has changed is that currently in the pure python implementation the flip proportion probability is applied to all lattice sites. I.e. if flip prop = 10%, then in each step all lattice sites have a 10% change of being chosen for the possibliity to flip. Therefore in each step there is a binomial distrubtion centered around (flip_prop * n) for the number of lattice sites that could flip. In the C++ library, exactly (flip_prop * n) random lattice sites are selected to possibly flip in each step. The change is for performance reasons (many less calls to the RNG) and do not change the resulting validity of the model. 
