"""
@author: Dharun Venkat
"""

import numpy as np
from LaplaceEquationSolver import LaplaceEquationSolver
from ParticleTrajectory import ParticleTrajectory



# running the particle trajectory simulation


ProtonMass = 1.67e-27
ProtonCharge = 1.6e-19
Proton = ParticleTrajectory(
    acceleration = np.array([0, 0]),
    name = "Proton",
    mass = ProtonMass,
    charge = ProtonCharge
)

LaplaceEquationSolver = LaplaceEquationSolver()
V, x_max, y_max = LaplaceEquationSolver.solve_and_plot()
Ex, Ey = LaplaceEquationSolver.calculate_electric_field(V)
#Ex, Ey = LaplaceEquationSolver.interpolation_func(V) - didn't get round to solving this



# the following lines are the inputs that can be changed

deltaT = 1e-5 # time step interval
Nt = int(500) # number of iterations to simulate 
my_method = 'Euler' # choose a numerical approximation technique to use here 
initial_position = np.array([50, 50]) # specify the starting location and velocity of the particle
initial_velocity = np.array([0, 0])

ParticleTrajectory = ParticleTrajectory()
positions = ParticleTrajectory.update(Ex, Ey, x_max, y_max, initial_position, initial_velocity, deltaT, Nt, my_method)
ParticleTrajectory.plot_trajectory(positions)
