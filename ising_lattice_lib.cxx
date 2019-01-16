#include <random>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
 
using namespace std;

struct lattice_site{
    /* ~lattice_site(){ */ 
        /* cout << "del lattice" << endl; */
        /* delete left; delete right; delete up; delete down; }; */
    /* ~lattice_site(){ delete left; delete right; delete up; delete down; }; */
	bool  spin;       // whether it is spin up (true) or spin down (false)
	short n_bonds;    // number of aligned spins with neighbor (+1 aligned, -1 anti-aligned)
	bool  bond_right; // true if spin aligned with site to right
	bool  bond_down;  // true is spin aligned with site down
    lattice_site  *left;
    lattice_site  *right;
    lattice_site  *up;
	lattice_site  *down;
};

class IsingMatrix{
    // Make it of lattice sites (ls), each of which is linked to four other 
    // sites: up/down/left/rigth: U/D/L/R.

    // Each site "owns" the boolean of opposing flip with the site next to
    // it on the D/R sides

	// Each LS has a pointer to each four adjacent sites


	// Each site also has pointers to the each of the sites 
    public:
    IsingMatrix(int set_N = 10, int set_Npick = 10); // same as allocate
    ~IsingMatrix();
    
    float get_E() { return -(J*Nbond+B*Nspin)/NN; };
    float get_M() { return (float)Nspin/NN; };
    int   get_N() { return N; };
    int   set_Npick(int in_val) { Npick = in_val; return 0; };
    int   set_J (float in_val) { J = in_val; return 0;} ;
    int   print_E() { cout << " Value of E: " << get_E() << endl; return 0; };
    int   print_Nspin() { cout << " #spin: " << Nspin << endl; return 0; };
    int   get_Nspin() { return Nspin; };
    int   get_Nbond() { return Nbond; };
    int   print_Nbond() { cout << " #bonds: " << Nbond << endl; return 0; };
    int   get_spin(int i, int j) { return (site[i*N + j].spin ? 1 : -1); };
    int   set_spin(int i, int j, int k) ;
    void  flip_site(int i);
          /* { site[i*N + j].spin = ((k > 0) ? true : false); return 0; }; */
    int   print_spins();     // defined
    int   print_bonds();     // defined
    int   initialize_spins(); // defined
    int   step(float T, float B);
    int   nsteps(int nsteps, float T, float B);
    int   calc_auto_correlation();
    float *auto_correlation;

    /* bool *spin; */
    /* short int *pairs; */

    private:
	lattice_site* site;

    int N;
    int NN;
    int Nbond;
    int Nspin;
    float J = 1.0;
    float flip_prob = 0.1;
    float E;
    float M;
    float B;
    /* std::random_device rd; */
    std::mt19937 gen;
    std::uniform_real_distribution<float> dis;
    int flip_spin(int i);

    vector<int> indices; // used to pick the order of latice sights
    int Npick;
};

int IsingMatrix::set_spin(int ii, int jj, int kk){
    int i = ii*N + jj;
    bool set_spin_state = kk > 0;
    if (site[i].spin != set_spin_state) flip_site(i);
    return 1;
};

void IsingMatrix::flip_site(int i){
    Nbond -= 2*site[i].n_bonds;
    Nspin += ( site[i].spin ? -2 : +2 );
  
    site[i].spin    ^= true;
    site[i].n_bonds *= -1;
    
    site[i].right->n_bonds += ( site[i].bond_right ? -2 : 2);
    site[i].bond_right     ^= true;
  
    site[i].down->n_bonds += ( site[i].bond_down ? -2 : 2);
    site[i].bond_down     ^= true;
  
    site[i].up->n_bonds   += ( site[i].up->bond_down ? -2 : 2);
    site[i].up->bond_down ^= true;
  
    site[i].left->n_bonds    += ( site[i].left->bond_right ? -2: 2);
    site[i].left->bond_right ^= true;
}

int IsingMatrix::print_spins() {
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            cout << " " << site[i*N+j].spin;
        }
        cout << endl;
    };
    return 0;
}
int IsingMatrix::print_bonds() {
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            printf("%3i ", site[i*N+j].n_bonds);
            /* cout << " " << site[i*N+j].n_bonds; */
        }
        cout << endl;
    };
    return 0;
}

int IsingMatrix::initialize_spins(){
	Nspin = 0;
    for (int i = 0; i < NN; ++i){
        if (dis(gen) > 0.5){
            Nspin += 1;
            site[i].spin = true;
        } else {
            Nspin -= 1;
            site[i].spin = false;
        }
    };

    for (unsigned int i = 0; i < NN; ++i){
        site[i].n_bonds = 0;
    }

	Nbond = 0;
	for (unsigned int i = 0; i < NN; ++i){
		if (site[i].spin ^ site[i].right->spin)	{
			site[i].bond_right = false;
			site[i].n_bonds -= 1;
			site[i].right->n_bonds -= 1;
			Nbond -= 1;	
		} else {
			site[i].bond_right = true;
			site[i].n_bonds += 1;
			site[i].right->n_bonds += 1;
			Nbond += 1;
		}

		if (site[i].spin ^ site[i].down->spin)	{
			site[i].bond_down = false;
			site[i].n_bonds -= 1;
			site[i].down->n_bonds -= 1;
			Nbond -= 1;	
		} else {
			site[i].bond_down = true;
			site[i].n_bonds += 1;
			site[i].down->n_bonds += 1;
			Nbond += 1;
		}
	}		
    return 0;
}


IsingMatrix::IsingMatrix(int set_N, int set_Npick) 
  : N{set_N}, gen{std::mt19937(time(0))}, dis {0.0, 1.0}, Npick{set_Npick}
{ 
    NN = N*N;
    auto_correlation = new float[N/2 - 1];
    for (int i = 0; i < (N/2-1); ++i) auto_correlation[i] = 0;
	site = new lattice_site[NN];
	for (int i = 0; i < NN; ++i){
		indices.push_back(i);
	}

	// set up the links between the sites 
	//link left edge to right edge
	for (int i = 0; i < N; ++i){
		site[i*N].left        = &site[(i+1)*N-1];
		site[(i+1)*N-1].right = &site[i*N];
	}
	//link left to right, all other columns
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N - 1; ++j){
			site[i*N+j].right  = &site[i*N+j+1];
			site[i*N+j+1].left = &site[i*N+j];
		}
	}
	//link top edge to bottom edge
	for (int j = 0; j < N; ++j){
		site[j].up           = &site[N*(N-1) + j];
		site[N*(N-1)+j].down = &site[j];
	}
	//link bottom to top, all other rows
	for (int i = 0; i < N - 1; ++i){
		for (int j = 0; j < N; ++j){
			site[i*N+j].down = &site[(i+1)*N + j];
			site[(i+1)*N+j].up    = &site[i*N+j];
		}
	}
	initialize_spins(); 
}

IsingMatrix::~IsingMatrix(){ 
    delete[] site;
    delete[] auto_correlation;
}

int IsingMatrix::calc_auto_correlation(){
    // get the row and column means
    vector<float> col_mean;
    vector<float> row_mean;
    for (int i = 0; i < N; ++i){
        vector<float> row_spin;
        vector<float> col_spin;
        for (int j = 0; j < N; ++j){
            row_spin.push_back( (float)get_spin(i,j) );
            col_spin.push_back( (float)get_spin(j,i) );
        }
        row_mean.push_back( accumulate(row_spin.begin(), row_spin.end(), 0.0)/N );
        col_mean.push_back( accumulate(col_spin.begin(), col_spin.end(), 0.0)/N );
    }
    /* cout << "row mean "; for (auto i : row_mean) cout << " " << i; cout << endl; */
    /* cout << "col mean "; for (auto i : col_mean) cout << " " << i; cout << endl; */
    for (int k = 1; k < N/2; ++k){
        float sum = 0;
        int n_entries = 0;
        for (int j = 1; j < N; ++j){
            for (int i = 0; i < N; ++i){
                int t = (j+k)%N;
                sum += (get_spin(j,i) - col_mean[i])*(get_spin(t,i)-col_mean[i]);
                sum += (get_spin(i,j) - row_mean[i])*(get_spin(i,t)-row_mean[i]);
                ++n_entries;
            }
        }
        sum /= (2*N*n_entries);
        /* sum /= (2*N); */
        auto_correlation[k-1] = sum;
    }
    return 0;
}

int IsingMatrix::nsteps(int n, float T, float B){
    for (int i = 0; i < n; ++i){
        step(T, B);
    }
    return 0;
}

int IsingMatrix::step(float T, float B_in){
    B = B_in;
    
	// randomly pick Npick sites to potentially flip
	int n_flipped = 0;
	/* int num_random = Npick; */
	auto begin = indices.begin();
	auto end   = indices.end();
	size_t left = std::distance(begin, end);
	/* cout << " left " << left << endl; */
	/* cout << "begin" << endl; */
	/* cout << indices[0] << endl; */
	for (int i = 0; i < Npick; ++i){
	/* while(num_random--){ */
		auto r = begin;	
		float temp = (int)(dis(gen)*left);
		/* cout << " " << temp; */
		/* advance(r, (int)(dis(gen)*left)); // find a random number past begin(), */
		advance(r, temp); // find a random number past begin(),

		/* cout << " picked " << *r; */

		// the site r is allowed to flip if:
        //    1. deltaE < 0 or 2. exp(-1.0*deltaE/T) > rand()
		float deltaE = 2*(J*site[*r].n_bonds+B*(2*site[*r].spin-1));
		/* cout << " " << *r << " spin: " << site[*r].spin << " " << deltaE << endl; */
		/* cout << " deltaE: " << deltaE << " -> exp(-1*deltaE/T)=" << exp(-1.0*deltaE/T) << " "; */
		if (deltaE < 0.0 || dis(gen) < exp(-1.0*deltaE/T)){
			/* cout << " YES swap " << endl; */
			//move r to the front
			swap(*begin, *r);
			++begin;
			++n_flipped;
		} else {
			/* cout << " NO swap " << endl; */
		    // move r to the end, where it cannot get picked
			swap(*(end-1), *r);
		}
		--left;
		/* cout << " begin-end " << *begin << " " << *(begin + left -1) << endl; */
		/* cout << " v: "; for (auto k : indices) cout << " " << k; cout << endl; */
	}
	/* cout << endl; */
	/* cout << " indices: "; for (auto i : indices) cout << " " << i; cout << endl; */
	// now the first Npick numbers in indeces are
	// the sites that are allowed to flip.
	// check them all at once to see if they should flip,
	// and then flip them.
	/* cout << " index flipped: "; */
	for (int index = 0; index < n_flipped; ++index){	
		int i = indices[index];
		flip_site(i);
	}
		
	return 0;	
}

// make a C-library
extern "C" // Tells the compiler to use C-linkage for this scope.
           // These are the functions available to ctypes in python
{
    /* void* c_new_IsingMatrix( void ) { return new(std::nothrow) IsingMatrix; } */
    void* newMatrix(int N, int n_flip ) { 
        return new IsingMatrix(N, n_flip);
    };
    int delMatrix( IsingMatrix* im) { delete im; return 0; }
    float get_E( IsingMatrix* im )  { return im->get_E(); }
    float get_M( IsingMatrix* im )  { return im->get_M(); }
    int   set_Npick( IsingMatrix* im, int val ) {
        return im->set_Npick( val ); }
    int   print_E( IsingMatrix* im )   { return im->print_E(); };
	int   print_Nspin( IsingMatrix* im ) { return im->print_Nspin(); };
	int   get_Nspin( IsingMatrix* im ) { return im->get_Nspin(); };
	int   get_Nbond( IsingMatrix* im ) { return im->get_Nbond(); };
	int   print_Nbond( IsingMatrix* im ) { return im->print_Nbond(); };
    int   get_spin( IsingMatrix* im, int i, int j) { return im->get_spin(i,j); };
    int   set_spin( IsingMatrix* im, int i, int j, int k ) { return im->set_spin(i,j,k); };
    int   print_spins( IsingMatrix* im ) { return im->print_spins(); };     // defined
    int   print_bonds( IsingMatrix* im ) { return im->print_bonds(); };     // defined
    int   initialize_spins( IsingMatrix* im ){ return im->initialize_spins(); }; // defined
    int   step( IsingMatrix* im, float T, float B) { return im->step(T, B); };
    int   nsteps( IsingMatrix* im, int i, float T, float B) { return im->nsteps(i, T, B); };
    int   calc_auto_correlation( IsingMatrix* im ) { return im->calc_auto_correlation(); };
    float auto_correlation( IsingMatrix* im, int i ){
        if (i >= im->get_N()/2+1){
            cout << "warning, there are only " 
                 << im->get_N()/2 
                 << " autocorrelation values available. " << endl;
            return 0;
        }
        return im->auto_correlation[i-1]; 
    };
};


int main(){
    cout << "starting main" << endl;
    IsingMatrix* z = reinterpret_cast<IsingMatrix*>(newMatrix(2,2));
    step(z, 2.2, 1.0); 
    print_spins(z);
    delMatrix(z);
    /* delete z; */
    return 0;
}
