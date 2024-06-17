# Numerical Laplace Equation Solver

### Laplace Equation Solver

The LaplaceEquationSolver.py is the Python script that numerically solves Laplace's equation using the finite difference method on the 
PDE then a 4-point stencil method and relaxation method to iteratively solve for the electric potential. When the program runs, the user
can specify the grid size and input a limited range of boundary conditions. The program then plots the scalar potential as a contour map,
calculates the corresponding electric field, and plots it as a 2D vector arrow field and then as a contour map by taking its magnitude.
The user also has the option to retrieve the value of the potential at a particular point. 

This is the main code and should be what the user runs. The other files are just for testing.


### Analysis

Analysis.py is the Python code that was used to run unit tests on the result by calculating certain fractional errors.


### Particle Trajectory and Simulation

ParticleTrajectory.py is the Python file that will calculate the motion of a proton when placed inside the electric field produced
by the potential and then plot it.
ParticleSimulation.py is the corresponding Python script to the above and is the only one that should be tailored. One can
specify the starting location and velocity of the proton. They can also choose the time step and the number of time steps
the simulation runs for as well as the desired numerical integration technique used to solve for the proton's trajectory.


### V_field

This is an old Python code that hasn't been updated in a while. It is a rudimentary, initial version of LaplaceEquationSolver.py.
It was abandoned and is safe to ignore but is kept in to show the full history of the project.


## Credits

PHYS389 Modelling Project, Lancaster University, Physics Department. \
The Python code was written using Spyder.\
Author: Dharun Venkateswaran \
Date: 04/04/2024 \
Acknowledgements: Ian Bailey and Chris Arridge 


## Contact

Email: dharun.v2003@gmail.com
