"""
@author: Dharun Venkat
"""

import numpy as np
import matplotlib.pyplot as plt
class ParticleTrajectory:
    
    
    def __init__( # initialising all the values that are going to be used
        self,
        acceleration = np.array([0, 0], dtype=float),
        name ='Particle',
        mass = 1.0,
        charge = 1.0,
        e = 1.6e-19,
        
    ):
        self.mass = 1.67e-27
        self.name = name
        self.charge = 1.6e-19
        self.e = e
        
        self.acceleration = np.array(acceleration, dtype = float)
        
    
    
    def __str__(self): # represent the class objects as strings to make the code more readable and easier to debug
        return "Particle: {0}, Mass: {1}, Charge: {2}, Position: {3}, Velocity: {4}, Acceleration: {5}".format(
            self.name, self.mass, self.charge, self.position, self.velocity, self.acceleration
    )
        
        
    def update(self, Ex, Ey, x_max, y_max, initial_position, initial_velocity, deltaT, Nt, my_method):
        """
        This method uses the numerical approximation chosen by the user to calculate
        and update the position and velocities of the particle in the electric field and
        calculates them as 1x2 vector NumPy arrays. The positions are then returned and used
        to be plotted.

        Parameters
        ----------
        Ex/Ey : x_max x y_max NumPy array
            The electric field components calculated by the LaplaceEquationSolver class.
        x_max/y_max: integers
            The x and y size of the coordinate grid specified.
        initial_position: 1x2 NumPy array
            The initial position of the proton in the form (x,y).
        initial_velocity: 1x2 NumPy array
            The initial velocity of the proton in the form (v_x,v_y).
        deltaT : float 
            The time step interval.
        Nt : integer
            The number of time steps the simulation is run.
        my_method : string
            The desired numerical method to be used.

        Raises
        ------
        TypeError
            Raises an error if an incorrect method is chosen and lets the user know.

        Returns
        -------
        positions : list
            A list of the positions the particle has been at to then plot.

        """
        positions = [initial_position] # store the positions and velocities as a list of 1x2 NumPy arrays
        velocities = [initial_velocity] # positions stored as a list because that's easy to plot
        
        for i in range(Nt):
            
            #  the relevant electric field value is where the particle currently is
            
            coordinate = positions[-1] # store the current coordinate as can't call positions into electric field otherwise
            
            if coordinate[0] < 0 or coordinate[0] > x_max or coordinate[1] < 0 or coordinate[1] > y_max: # stops simulation if boundary is reached
                break
            
            Ex_x = Ex[int(coordinate[0]), int(coordinate[1])]
            Ey_y = Ey[int(coordinate[0]), int(coordinate[1])] 
            
            acceleration = (self.charge / self.mass) * np.array([Ex_x, Ey_y])
        
            if my_method == 'Euler':
                
                position = positions[-1] + deltaT * velocities[-1] # use the most recent position and velocity
                velocity = velocities[-1] + deltaT * acceleration
                
            elif my_method == 'Euler_Cromer':
                velocity = velocities[-1] + deltaT * acceleration
                position = positions[-1] + deltaT * velocities[-1]
                
            elif my_method == 'Verlet': # symplectic algorithm
                old_acceleration = acceleration
                position = positions[-1] + deltaT * velocities[-1] + 1/2 * acceleration * (deltaT ** 2)
                
                # update the acceleration 
                Ex_x = Ex[int(position[0]), int(position[1])]
                Ey_y = Ey[int(position[0]), int(position[1])] 
                
                acceleration = (self.charge / self.mass) * np.array([Ex_x, Ey_y])
                
                velocity = velocities[-1] + 1/2 * (old_acceleration + acceleration) * deltaT
                
            elif my_method == 'RK-2':
                position_half = positions[-1] + 1/2 * deltaT * velocities[-1]
                velocity_half = velocities[-1] + 1/2 * deltaT * acceleration
                position = positions[-1] + velocity_half * deltaT
                
                # update the acceleration 
                Ex_x = Ex[int(position_half[0]), int(position_half[1])]
                Ey_y = Ey[int(position_half[0]), int(position_half[1])]                
                acceleration_half = (self.charge / self.mass) * np.array([Ex_x, Ey_y])
                
                velocity = velocities[-1] + acceleration_half * deltaT
            
           
            else:
                raise TypeError("You have not chosen a valid numerical approximation method.")
                
            positions.append(position) # updating the position and velocity
            velocities.append(velocity)
            
        return positions
        
    
    def plot_trajectory(self, positions):
        """
        Plots the trajectory of the proton given a series of positions.

        Parameters
        ----------
        positions : list 
            A list containing the positions of the proton 
            at different time steps over the course of the simulation. 

        Returns
        -------
        None.
        
        """
        positions = np.array(positions)
        
        plt.plot(positions[:, 0], positions[:, 1], label = 'Trajectory')
        plt.scatter(positions[0, 0], positions[0, 1], color='red', label = 'Initial Position', marker = 'o')
        plt.scatter(positions[-1, 0], positions[-1, 1], color = 'green', label = 'Final Position', marker = 'o')
        plt.title('Proton Particle Trajectory')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.axis('equal') # set scaling to be equal
        plt.legend()
        plt.show()
