import sys
import os
import time
import csv
import click
import numpy as np
import logging
import matplotlib.pyplot as plt
from IsingLattice import IsingLattice
from tqdm import tqdm #fancy progress bar generator
from ising_c import run_ising #import run_ising function from ising.py

def calculate_and_save_values(lattice, Msamp,Esamp,num_analysis,index,temp,data_filename,corr_filename):

    try:
        #calculate statistical values
        M_mean = np.average(Msamp[-num_analysis:])
        E_mean = np.average(Esamp[-num_analysis:])
        M_std = np.std(Msamp[-num_analysis:])
        E_std = np.std(Esamp[-num_analysis:])
        data_array = [np.abs(M_mean),M_std,E_mean,E_std]
    
        #write data to CSV file
        header_array = ['Temperature','Magnetizatio n Mean','Magnetization Std Dev','Energy Mean','Energy Std Dev']
        append_data_to_file(data_filename, header_array) if index == 0 else None
        append_data_to_file(data_filename, data_array, temp)

        #get correlation function
        # corr = compute_autocorrelation(spin)
        corr = lattice.calc_auto_correlation()

        #write correlation function to CSV file
        header_array = ['Temperature','K','Spatial Spin Correlation']
        append_data_to_file(corr_filename, header_array) if index == 0 else None
        [append_data_to_file(corr_filename, corr_value, temp) for corr_value in corr]

        return True

    except:
        logging.error("Temp="+str(temp)+": Statistical Calculation Failed. No Data Written.")
        return False

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
@click.option('--plots', default=True, help='Turn Automatic Plot Creation Off or On',type=bool)

def run_simulation(t_min,t_max,t_step,n,num_steps,num_analysis,num_burnin,j,b,flip_prop,output,plots):

    time_start = time.time()

    check_step_values(num_steps, num_analysis, num_burnin)

    lattice = IsingLattice(n, flip_prop)
    # lattice.print_spins()

    T = get_temp_array(t_min, t_max, t_step)

    data_filename, corr_filename = get_filenames(output)

    write_sim_parameters(data_filename,corr_filename,n,num_steps,num_analysis,num_burnin,j,b,flip_prop)

    if plots:
        #initialize vars for plotting values
        temp_arr, M_mean_arr, E_mean_arr, M_std_arr, E_std_arr = [],[],[],[],[]

    print('\nSimulation Started! Data will be written to ' + data_filename + '\n')

    temp_range = tqdm(T) #set fancy progress bar range
    for index, temp in enumerate(temp_range):

        #show current temperature
        temp_range.set_description("Simulation Progress");

        try:
            #run the Ising model
            Msamp, Esamp = run_ising(lattice, temp,num_steps,num_burnin,j,b)
            #get and save statistical values
            if calculate_and_save_values(lattice, Msamp,Esamp,num_analysis,index,temp,data_filename,corr_filename):
                if plots:
                    #for plotting
                    M_mean, E_mean, M_std, E_std = get_plot_values(temp,Msamp,Esamp,num_analysis)
                    temp_arr.append(temp)
                    M_mean_arr.append(M_mean)
                    E_mean_arr.append(E_mean)
                    M_std_arr.append(M_std)
                    E_std_arr.append(E_std)

        except KeyboardInterrupt:
            print("\n\nProgram Terminated. Good Bye!")
            sys.exit()

        except:
            logging.error("Temp="+str(temp)+": Simulation Failed. No Data Written")


    print('\n ------ Time start(%f) Time finished(%f), total time(%f)'%(time_start, time.time(), time.time()-time_start))
    print('\n\nSimulation Finished! Data written to '+ data_filename)

    lattice.free_memory()
    if plots:
        plot_graphs(temp_arr, M_mean_arr, M_std_arr, E_mean_arr, E_std_arr)

def get_plot_values(temp,Msamp,Esamp,num_analysis): #only for plotting at end
    try:
        M_mean = np.average(Msamp[-num_analysis:])
        E_mean = np.average(Esamp[-num_analysis:])
        M_std = np.std(Msamp[-num_analysis:])
        E_std = np.std(Esamp[-num_analysis:])
        return M_mean, E_mean, M_std, E_std
    except:
        logging.error("Temp={0}: Error getting plot values".format(temp))
        return False

def plot_graphs(temp_arr,M_mean_arr,M_std_arr,E_mean_arr,E_std_arr): #plot graphs at end
    plt.figure(1)
    plt.ylim(0,1)
    plt.errorbar(temp_arr, np.absolute(M_mean_arr), yerr=M_std_arr, uplims=True, lolims=True,fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Magnetization')
    plt.figure(2)
    plt.errorbar(temp_arr, E_mean_arr, yerr=E_std_arr, fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Energy')
    plt.show()

def check_step_values(num_steps,num_analysis,num_burnin): #simulation size checks and exceptions
    if (num_burnin > num_steps):
        raise ValueError('num_burning cannot be greater than available num_steps. Exiting simulation.')
        sys.exit()

    if (num_analysis > num_steps - num_burnin):
        raise ValueError('num_analysis cannot be greater than available num_steps after burnin. Exiting simulation.')
        sys.exit()

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

def append_data_to_file(filename,data_array,temp=False):
    try:
        with open(filename,'a') as csv_file: #appends to existing CSV File
            writer = csv.writer(csv_file, delimiter=',', lineterminator='\n')
            if temp:
                writer.writerow([temp]+data_array)
            else:
                writer.writerow(data_array)

    except:
        logging.error("Temp={0}: Error Writing to File".format(temp))

if __name__ == "__main__":
    print("\n2D Ising Model Simulation\n")
    run_simulation()
