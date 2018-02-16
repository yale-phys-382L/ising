import sys
import os
import time
import csv
import click
import numpy as np
import logging
import matplotlib.pyplot as plt
from ising_c import run_ising #import run_ising function from ising.py
import multiprocessing as mp
from IsingLattice import IsingLattice

def run_simulation(index,temp,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,data_listener,corr_listener):
    print("Working on Temp {0}".format(round(temp,4)))
    try:
        #run the Ising model
        lattice = IsingLattice(n, flip_prop)
        time_start = time.time()
        Msamp, Esamp = run_ising(lattice,temp,num_steps,num_burnin,j,b,disable_tqdm=True)

        try:
            #calculate statistical values
            M_mean = np.average(Msamp[-num_analysis:])
            E_mean = np.average(Esamp[-num_analysis:])
            M_std = np.std(Msamp[-num_analysis:])
            E_std = np.std(Esamp[-num_analysis:])

            data_array = [np.abs(M_mean),M_std,E_mean,E_std]
            data_listener.put([temp]+data_array)

            corr = lattice.calc_auto_correlation()
            lattice.free_memory()
            [corr_listener.put([temp]+corr_value) for corr_value in corr]

            print("Done with Temp {0} in {1} seconds".format(round(temp,4), round(time.time()-time_start,2)))
            return True

        except:
            logging.error("Temp="+str(temp)+": Statistical Calculation Failed. No Data Written.")
            return False

    except KeyboardInterrupt:
        print("\n\nProgram Terminated. Good Bye!")
        data_listener.put('kill')
        corr_listener.put('kill')
        sys.exit()

    except:
        logging.error("Temp="+str(temp)+": Simulation Failed. No Data Written")

#simulation options (enter python main.py --help for details)
@click.command()
@click.option('--t_min', default=2.0, prompt='Minimum Temp', help='Minimum Temperature (inclusive)', type=float)
@click.option('--t_max', default=2.6, prompt='Maximum Temp', help='Maximum Temperature (inclusive)', type=float)
@click.option('--t_step', default=0.1, prompt='Temp Step Size', help='Temperature Step Size', type=float)

@click.option('--n', prompt='Lattice Size', help='Lattice Size (NxN)',type=int)
@click.option('--num_steps', default=100000, help='Total Number of Steps',type=int)
@click.option('--num_analysis', default=50000, help='Number of Steps used in Analysis',type=int)
@click.option('--num_burnin', default=25000, help='Total Number of Burnin Steps',type=int)

@click.option('--j', default=1.0, help='Interaction Strength',type=float)
@click.option('--b', default=0.0, help='Applied Magnetic Field',type=float)
@click.option('--flip_prop', default=0.1, help='Proportion of Spins to Consider Flipping per Step',type=float)

@click.option('--output', default='data', help='Directory Name for Data Output',type=str)

@click.option('--processes', default=1, help='',type=int)

def main(t_min,t_max,t_step,n,num_steps,num_analysis,num_burnin,j,b,flip_prop,output,processes):
    simulation_start_time = time.time()
    data_filename, corr_filename = initialize_simulation(n,num_steps,num_analysis,num_burnin,output,j,b,flip_prop)
    run_processes(processes,t_min,t_max,t_step,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename)
    simulation_duration = round((time.time() - simulation_start_time)/60.0,2)
    print('\n\nSimulation finished in {0} minutes. Data written to {1}.'.format(simulation_duration,data_filename))

def initialize_simulation(n,num_steps,num_analysis,num_burnin,output,j,b,flip_prop):
    check_step_values(num_steps, num_analysis, num_burnin)
    data_filename, corr_filename = get_filenames(output)
    write_sim_parameters(data_filename,corr_filename,n,num_steps,num_analysis,num_burnin,j,b,flip_prop)
    print('\nSimulation Started! Data will be written to ' + data_filename + '\n')
    return data_filename, corr_filename

def check_step_values(num_steps,num_analysis,num_burnin): #simulation size checks and exceptions
    if (num_burnin > num_steps):
        raise ValueError('num_burning cannot be greater than available num_steps. Exiting simulation.')

    if (num_analysis > num_steps - num_burnin):
        raise ValueError('num_analysis cannot be greater than available num_steps after burnin. Exiting simulation.')

def get_filenames(dirname): #make data folder if doesn't exist, then specify filename
    try:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        data_filename = os.path.join(dirname,'data_'+str(time.strftime("%Y%m%d-%H%M%S"))+".csv")
        corr_filename = os.path.join(dirname,'corr_'+str(time.strftime("%Y%m%d-%H%M%S"))+".csv")
        #Write simulation parameters to file
        return data_filename, corr_filename
    except:
        raise ValueError('Directory name not valid. Exiting simulation.')
        sys.exit()

def get_temp_array(t_min,t_max,t_step):
    if (t_min > t_max):
        raise ValueError('T_min cannot be greater than T_max. Exiting Simulation')
        sys.exit()
    try:
        T = np.arange(t_min,t_max,t_step).tolist()
        return T
    except:
        raise ValueError('Error creating temperature array. Exiting simulation.')
        sys.exit()

def write_sim_parameters(data_filename,corr_filename,n,num_steps,num_analysis,num_burnin,j,b,flip_prop):
    try:
        with open(data_filename,'w') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', lineterminator='\n')
            #Write simulations parameters to CSV file
            writer.writerow(['Lattice Size (NxN)','Total Steps','Steps Used in Analysis','Burnin Steps','Interaction Strength','Applied Mag Field','Spin Prop'])
            writer.writerow([n,num_steps,num_analysis,num_burnin,j,b,flip_prop])
            writer.writerow([])
        with open(corr_filename,'w') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', lineterminator='\n')
            #Write simulations parameters to CSV file
            writer.writerow(['Lattice Size (NxN)','Total Steps','Steps Used in Analysis','Burnin Steps','Interaction Strength','Applied Mag Field','Spin Prop'])
            writer.writerow([n,num_steps,num_analysis,num_burnin,j,b,flip_prop])
            writer.writerow([])
    except:
        logging.error('Could not save simulation parameters. Exiting simulation')
        sys.exit()

def compute_autocorrelation(spin):
    n = len(spin)
    corr_array = []
    for k in range(1,int(n/2)):
        col_mean, row_mean = spin.mean(axis=0),spin.mean(axis=1)
        #compute r values for rows and cols
        r_col = [np.multiply(spin[j,:]-col_mean,spin[(j+k)%n,:]-col_mean) for j in range(1,n)]
        r_row = [np.multiply(spin[:,j]-row_mean,spin[:,(j+k)%n]-row_mean) for j in range(1,n)]
        #normalize r values
        r_col = np.divide(r_col,float(n))
        r_row = np.divide(r_row,float(n))
        #calculate corr for k and add it to array
        corr = (r_col.mean() + r_row.mean())/2.0
        corr_array.append([k,corr])
    return corr_array

def listener(q, fn):
    '''listens for messages on the q, writes to file. '''
    f = open(fn, 'a') 
    writer = csv.writer(f, delimiter=',', lineterminator='\n')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        writer.writerow(m)
        f.flush()
    f.close()

def run_processes(processes,t_min,t_max,t_step,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename):
    
    T = get_temp_array(t_min, t_max, t_step)
    
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    data_listener = manager.Queue()
    corr_listener = manager.Queue()    
    pool = mp.Pool(mp.cpu_count() + 2)

    #put listener to work first
    data_watcher = pool.apply_async(listener, args=(data_listener, data_filename,))
    corr_watcher = pool.apply_async(listener, args=(corr_listener, corr_filename,))

    #fire off workers 
    jobs = [pool.apply_async(run_simulation, args=(index,temp,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,data_listener,corr_listener,)) for index,temp in enumerate(T)]

    # collect results from the workers through the pool result queue   
    [job.get() for job in jobs]

    #now we are done, kill the listener
    data_listener.put('kill')
    corr_listener.put('kill')
    pool.close()

if __name__ == "__main__":
   main()
