# 2D Ising Model in Python

## New Features
* Command line interface to input simulation parameters
* Better interpretability with pythonic features (e.g. list comprehension)
* Modular codebase for easy changes and experimentation
* Fancy progress bars with time estimations
* Complete error handling with progress save
* ...and more coming!

## Getting Started

#### Installing and using Python

If you don't have Python installed yet, I would highly recommend using the **Anaconda distribution** to install Python 3. You can find the installation instructions [here](https://docs.anaconda.com/anaconda/install/)

Once this is installed on your computer, you will have Python 3 ready to go, as well as important packages like NumPy and SciPy. You can view these packages and install new ones using the [Anaconda Navigator](https://docs.anaconda.com/anaconda/navigator/) (need to install this separately).

In order to edit and run your code, I would recommend [Spyder (Scientific PYthon Development EnviRonment)](https://pythonhosted.org/spyder/) (I know, horrible acronym, but the IDE makes up for it). It should be really easy to edit your code and run it through this environment. **NOTE**: The IPython shell in Spyder does not support nested progress bars, so you will only see one when you run the simulation. In order to see both, you will need to change your run configuration to run in a normal Python shell.

Another optional but cool program you can use is [Jupyter Notebook](http://jupyter.org/) (comes pre-installed with Anaconda). These notebooks support Python code, as well as Markdown and LaTeX, so you can keep all of your code organized and easily testable (hint: use for easier data analysis!). You should be able to open this through Anaconda Navigator.

These are just recommendations - there are plenty of other GUI-based applications for Python development out there (like Enthought Canopy). If you have time, do some playing around and see what you like!

### Command line (Mac/Linux)

If you want to change other parameters of the simulation, you can use:

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


## Understanding the Simulation

There are three important python files in this simulation: `main.py`, `ising.py`, and `annealing.py`

`main.py` is the file you have to run for the simulation. The code in this file takes in the input parameters, runs the Ising model for each temperature step, gets the relevant data, saves it, and gives you a set of nice plots at the end. This is a lot, so we've broken this down into different functions to make it easier to understand/change. Here are the two most important ones:

* `run_simulation`: takes in all the input variables and runs the simulation.

* `calculate_and_save_values`: takes in the energy, magnetization, and spin values from the Ising code, calculates the appropriate statistical values, and saves them to a CSV file. **This is where you should implement code to calculate the other values you are interested in**.

`ising.py` calculates the Ising model at a certain temperature

## Multi-processing

If you feel adventurous and want to use multiprocessing (running the simulation on multiple cores), we got you covered! Use the `main-multiprocessing.py` file to get started. This is almost identical to the `main.py` file.

The code will automatically use all the cpu cores you have (you can change this in the code)

There isn't a clean way to add progress bars to this, so instead, the program will output the temperatures it is currently computing, as well as the temperatures it is finished with.

**NOTE:** Because of the asynchronous conditions of this program, it may not write your data to the CSV file in the right order. Make sure you sort by temperature before analyzing the data.

## Acknowledgements

Thank you to Jed Thompson for writing the original Matlab code, and to Nir Navon, Reina Maruyama, and Steve Lamoreaux for their help.