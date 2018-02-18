# 2D Ising Model

Python code written by Surya Dutta

C library written by Dave Stewart 

## Getting Started

### GUI-based approach (Windows/Mac/Linux)

If you prefer to use GUIs as opposed to the command line, this section is for you!

#### Managing your Git Repository (optional, but recommended)

If you would like to use Git Version Control in your team to collaborate on and keep backups of your code, great! If not, no worries - just download the files here and follow the instructions below. Version control is always highly recommended.

There are plenty of great GUIs for Git. My personal favorite is [Github Desktop](https://desktop.github.com/).

If you are new to Git and want to learn more about version control, visit [this website](https://programminghistorian.org/lessons/getting-started-with-github-desktop) for a great primer on version control, Git, and Github Desktop.

The first step is to make a Github account and fork this repository (click on `Fork` in the top right). This will create a copy of this code onto your own account. Now you can follow the instructions for your respective GUI to clone this repository (download the files locally), and start working with the simulation!

### Command line (Mac/Linux)

1. Fork this repository to your own account

2. Navigate to the folder you would like to use, then use:
  ```bash
  git clone git@github.com:{{your-github-username}}/ising && cd ising
  ```

3. If you are using conda (recommended), use this to install required packages:
  ```bash
  conda install --yes --file requirements.txt
  ```

  If you have a standalone version of Python 3 installed and are using pip, use this instead (you may need to be a superuser to install Pip packages):
 ```bash
  pip install -r requirements.txt
  ```

4. The code should be ready to run! Use this to run the simulation, and it should save the results automatically to an auto-generated data folder:
  ```bash
  python main.py
  ```

5. If you want to change other parameters of the simulation, you can use:

  ```
  python main.py --help

  2D Ising Model Simulation
  Usage: main.py [OPTIONS]

  Options:
    --t_min FLOAT           Minimum Temperature (inclusive)
    --t_max FLOAT           Maximum Temperature (inclusive)
    --t_step FLOAT          Temperature Step Size
    --n INTEGER             Lattice Size (NxN)
    --num_steps INTEGER     Total Number of Steps
    --num_analysis INTEGER  Number of Steps used in Analysis
    --num_burnin INTEGER    Total Number of Burnin Steps
    --j FLOAT               Interaction Strength
    --b FLOAT               Applied Magnetic Field
    --flip_prop FLOAT       Proportion of Spins to Consider Flipping per Step
    --help                  Show this message and exit.
  ```

  This will list all of the parameters you can change. For example, if you run `python main.py --b=0.5 --flip_prop=0.2`, the simulation will add a magnetic field of 0.5T and increase the flip proportion to 0.2. You can also edit the default parameters directly in the `main.py` file.

  ### Some general benchmarks for code performance

  Measured in average execution time for each temperature step (lower the better):

  * Regular Python = 123.65 secs

  * Python - 4 processes = 30.92 secs

  * Regular Python with C = 12.38 secs

  * Python with C - 4 processes = 1.23 secs