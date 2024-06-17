"""
@author: Dharun Venkat
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline


class LaplaceEquationSolver:
    
    def __init__(self):
        
        self.x_min, self.x_max, self.y_min, self.y_max, self.dx, self.x, self.y, self.xv, self.yv = self.grid_spacing()
        self.boundary_conditions = self.set_boundary_conditions()
        

    def grid_spacing(self):
        """
        Returns the list of x and y coordinates that make up the discrete grid as well as the coordinate matrices
        xv and yv to be used in other functions that is specified by the user via inputs. Raises errors and restarts the
        programme if a nonvalid input is attempted.

        Parameters
        ----------
        None.
        
        
        Raises
        ---------
        ValueError
            If the user does not not input an integer for the end points of the grid.
        
        
        Returns
        ---------
        x_min : integer
            The starting x coordinate.
        x_max : integer
            The end x coordinate.
        y_min : integer
            The starting y coordinate.
        y_max : integer
            The end y coordinate.
        dx = float
            The spacing between x and y coordinates.
        x : list
            The list of x coordinates
        y : list
            The list of y coordinates
        xv, yv: NumPy meshrgrid
            Matrix of x and y coordinates
        """
        
        self.x_min = 0
        
        while True: # loop until a valid input is given
            try:
                self.x_max = int(input("The start x coordinate is 0. Enter the ending x coordinate as an integer: "))
                break
            except ValueError:
                print("You have not typed in a valid ending x coordinate. It must be an integer.\n")
            
        
        self.y_min = 0
        
        while True: 
            try:
                self.y_max = int(input("The start y coordinate is 0. Enter the ending y coordinate as an integer: "))
                break
            except ValueError:
                print("You have not typed in a valid ending y coordinate. It must be an integer.\n")
            
            
        self.dx = 1 # standard integer spacing of 1
        
            
        self.x = np.arange(self.x_min, self.x_max + 1.0, self.dx) # returning evenly spaced values: a coordinate list
        self.y = np.arange(self.y_min, self.y_max + 1.0, self.dx)
        
        self.xv, self.yv = np.meshgrid(self.x, self.y) # returns coordinates matrices from coordinate vectors
        
        return self.x_min, self.x_max, self.y_min, self.y_max, self.dx, self.x, self.y, self.xv, self.yv
               

    def set_boundary_conditions(self):
        """
        Creates a dictionary list of the boundary conditions specified by the user.
        
        Parameters
        ----------
        None.
        
        
        Raises
        ---------
        TypeError
            If the user does not type in a valid type of boundary condition available.
        
        ValueError
            If the user types in a boundary condition outside of the specified grid.
        
        Returns
        ----------
        boundary_conditions: dictionary 
            Contains the boundary conditions specified by the user in the form {(i,j}: potential_value)}
        """
        
        
        boundary_conditions = {} # dictionary to store the boundary conditions specified
        
        print("Now that the coordinate grid has been created, you can now add in the boundary conditions that you wish to have. \n")
        
        while True:            
            while True:
                boundary_type = str(input("What type of boundary condition would you like to add? Type c for circle, l for line, p for point, or n to stop: "))
            
                if boundary_type in ["c", "l", "p", "n"]: # only boundary condition types accepted
                    break
                else:
                    print("You have not chosen a valid type of boundary condition. Make sure the letter is lowercase. \n")
            
            if boundary_type == 'n':
                break
            
            
            elif boundary_type == 'p':
                while True:
                    try:
                        point_x, point_y = map(int, input("Which point would you like the potential to be? Enter it like 'x y' where x and y are integers. \n").split())
                        
                        if point_x < self.x_min or point_x > self.x_max or point_y < self.y_min or point_y > self.y_max: # checking if the coordinates are within the grid boundaries
                            raise ValueError("Coordinates are outside the specified grid. Please try again. \n")
                        break  # breaks the loop if valid coordinates are provided
                    except ValueError:
                        print("Invalid coordinate inputs. Please make sure that they are both valid integers. \n")
                            
                while True:
                    try:
                        point_potential = float(input("Enter the value of the potential at this point: ")) # write a function to check if the potential is a float as this gets repeated
                        
                        break
                    except ValueError:
                        print("You have not entered a valid electric potential value. Make sure it is a float and try again. \n")
                
                extra_boundary_condition = {(int(point_x), int(point_y)): point_potential} 
                boundary_conditions.update(extra_boundary_condition)
                
            
            elif boundary_type == 'l':
                while True:
                    line_type = str(input("What type of line would you like? Enter d for diagonal, h for horizontal, or v for vertical: "))
                
                    if line_type in ["d", "h", "v"]:
                        break
                    else:
                        print("You have not chosen a valid type of line. Make sure the letter is lowercase. \n")
                
                
                if line_type == 'h':
                    while True:
                        try:
                            h_line = int(input("What value of y would you like to horizontal line to be? It must be an integer: "))                   
                            
                            if h_line < self.y_min or h_line > self.y_max:
                                raise ValueError("The y coordinate is outside of the grid. Please try again. \n")
                            break                         
                        except ValueError:
                            print("Invalid y input. Please make sure that it's an integer. \n")
                    
                    while True:
                        try:
                            h_line_potential = float(input("Enter the value of the potential along this line: "))
                            
                            break
                        except ValueError:
                            print("You have not entered a valid electric potential value. Make sure it is a float and try again. \n")
                            
                    extra_boundary_condition = {(int(i), int(h_line)): h_line_potential for i in range(len(self.y))}
                    boundary_conditions.update(extra_boundary_condition)
               
                    
                elif line_type == 'v':
                    while True:
                        try:
                            v_line = int(input("What value of x would you like to vertical line to be? It must be an integer: "))
                            
                            if v_line < self.x_min or v_line > self.x_max:
                                raise ValueError("The x coordinate is outside of the grid. Please try again. \n")
                            break
                        except ValueError:
                            print("Invalid x input. Please make sure that it's an integer. \n")
                            
                    while True:
                        try:
                            v_line_potential = float(input("Enter the value of the potential along this line: "))
                    
                            break
                        except ValueError:
                            print("You have not entered a valid electric potential value. Make sure it is a float and try again. \n")
                    
                    extra_boundary_condition = {(int(v_line), int(i)): v_line_potential for i in range(len(self.x))}
                    boundary_conditions.update(extra_boundary_condition)
                
                
                elif line_type == 'd':
                    while True:
                        try:
                            d_line_xstart, d_line_ystart = map(int, input("Which point would the potential line start at? Enter it like 'x y' where x and y are integers. \n").split())
                            
                            if d_line_xstart < self.x_min or d_line_xstart > self.x_max or d_line_ystart < self.y_min or d_line_ystart > self.y_max:
                                raise ValueError("Coordinates are outside of the specified grid. Please try again. \n")
                            break
                        except ValueError:
                            print("Invalid coordinate inputs. Please make sure that they are both valid integers. \n")
                    
                    while True:
                        try:
                            d_line_xend, d_line_yend = map(int, input("Which point would the potential line start at? Enter it like 'x y' where x and y are integers. \n").split())
                            if d_line_xend < self.x_min or d_line_xend > self.x_max or d_line_yend < self.y_min or d_line_yend > self.y_max:
                                raise ValueError("Coordinates are outside of the specified grid. Please try again. \n")
                            break
                        except ValueError:
                            print("Invalid coordinate inputs. Please make sure that they are both valid integers. \n")
                    
                    while True:
                        try:
                            d_line_potential = float(input("Enter the value of the potential along this line: "))
                            
                            break
                        except ValueError:
                            print("You have not entered a valid electric potential value. Make sure it is a float and try again. \n")
                    
                    
                    m = (d_line_yend - d_line_ystart) / (d_line_xend - d_line_xstart) # calculating the straight line gradient
                    c = d_line_ystart - (m * d_line_xstart) # calculating the straight line y-intercept
                    
                    for i in range(0, (d_line_xend - d_line_xstart + 1)):
                        
                        extra_boundary_condition = {(int(d_line_xstart + i), int(m * (d_line_xstart + i) + c)): d_line_potential}
                        boundary_conditions.update(extra_boundary_condition)
                    
                    
            elif boundary_type == 'c':
                while True:
                    try:
                        x_centre, y_centre = map(int, input("Which point would you like the circle to be cenetred? Enter it like 'x y' where x and y are integers. \n").split())
                        
                        if x_centre < self.x_min or x_centre > self.x_max or y_centre < self.y_min or y_centre > self.y_max:
                            raise ValueError("Coordinates are outside of the specified grid. Please try again. \n")
                        break
                    except ValueError:
                        print("Invalid coordinate inputs. Please make sure that they are both valid integers. \n")
                 
                while True:
                    try:
                        radius = float(input("What would the radius of the circle be? "))
                        
                        break
                    except ValueError:
                        print("You have not typed in a valid radius. It must be a float. \n")
                             
                while True:
                    try:
                        circle_potential = float(input("Enter the value of the potential along the circle: "))
                        
                        break
                    except ValueError:
                        print("You have not entered a valid electric potential value. Make sure it is a float and try again. \n")
                
                 
                
                theta = np.linspace(0, 2 * np.pi, 50) # arbitrary spacing chosen - can vary it
                x_theta = x_centre + radius * np.cos(theta)
                y_theta = y_centre + radius * np.sin(theta)
                 
                for i in range(0, len(x_theta)):
                     
                    extra_boundary_condition = {(int(x_theta[i]), int(y_theta[i])): circle_potential}
                    boundary_conditions.update(extra_boundary_condition)
                    
        
        return boundary_conditions


    def finite_difference(self, tolerance = 1e-5, max_iterations = 5000):
        """
        Returns the scalar potential field using the finite difference method of averaging the potential from
        the 4 nearby adjacent points - 4-point average stencil method. It iterates using the method of relaxation
        until the true solution is found.
        
        Parameters
        ----------
        tolerance : float
            Determines at what tolerance does the simulation stop updating the potential using the stencil method
        max_iterations : integer
            Determines the max amount of iterations the simulation will go through if the tolerance if never reached
            
        
        Returns
        ---------
        V : NumPy array
            The electric potential at every point on the discrete grid.
        """
        
        
        V = np.zeros_like(self.xv, dtype = np.float64) #  creating an array of 0s with same size as coordinate matrix to initialise
        
        for (i, j), potential_value in self.boundary_conditions.items(): # apply the boundary conditions to V
            V[j, i] += potential_value # this is swapped because of how NumPy arrays work: V[row = y_value, column = x_value]

        
        for n in range(0, max_iterations):
            
            V_prev = V.copy() # to store and compare to the next value whether the potential has significantly changed or not
            
            for i in range(1, len(self.x) - 1): #  running through the interior points only
                
                for j in range(1, len(self.y) - 1):

                    if (i, j) not in self.boundary_conditions: # make sure the boundary conditions don't get replaced
                        
                        V[j, i] = (V[j, i + 1] + V[j, i - 1] + V[j + 1, i] + V[j - 1, i]) / 4 # average of the 4 points around it - stencil method
            
                    
            if np.max(np.abs(V - V_prev)) < tolerance: # stops running if potential no longer changes after an iteration
                break
        
        
        return V


    def plot_potential_field(self, V):
        """
        Plots a contour map of the scalar electric potential field.

        Parameters
        ----------
        V : NumPy array
            The scalar potential field at all points in the grid specified.

        Returns
        -------
        None.

        """
        
        
        plt.contourf(self.xv, self.yv, V, 50) # can change number of contour lines if needs be
        plt.axis('scaled')
        plt.colorbar(label = 'Electric Potential (V)')
        plt.title('Electrostatic Potential Distribution')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.show()


    def get_potential_field_value(self, V):
        """
        Allows the user to retrieve the electric potential value at a particular
        point.

        Parameters
        ----------
        V : NumPy array
            The electric potential value at all points on the grid.

        Returns
        -------
        None.

        """
        
        
        print("Now that the potential field has been found, you can request to know the exact value of it at a particular coordinate. \n")
        
        while True:
            while True:
                ask = str(input("Would you like to know the potential at a particular point? Type y for yes, or n for no: "))
                
                if ask in ['y', 'n']:
                    break
                else:
                    print("You have not typed in either valid of the valid inputs 'y' or 'n'. Please try again. \n")
            
            if ask == 'n':
                break
            
            elif ask == 'y':
                while True: 
                    try:
                        x_coordinate, y_coordinate = map(int, input("Which point would you like to know the potential at? Enter it like 'x y' where x and y are integers. \n").split())
                        print("The potential at ({0},{1}) is {2}.".format(x_coordinate, y_coordinate, V[y_coordinate, x_coordinate]))
                        break
                                       
                    except ValueError:
                        print("Those are not valid coordinates for this grid. Please try again.")
                            

    #def interpolation_func(self, V): #  never got round to finishing this off
        

    #    interpolate_func = RectBivariateSpline(self.x, self.y, V)
        
    #    Ex = interpolate_func(self.x, self.y, dx = 1, dy = 0)
    #    Ey = interpolate_func(self.x, self.y, dx = 0, dy = 1)
        
    #    return Ex, Ey



    def calculate_electric_field(self, V):
        """
        Calculates the electric vector field from the scalar potential field using the NumPy grad function.
        
        Parameters
        ----------
        V : NumPy array
            The scalar potential field calculated previously using finiteDifference.
            
        Returns
        ---------
        Ex, Ey : NumPy array
            The x and y electric field obatined using the scalar potential field through E = -grad(V).
        """
        
        
        Ey, Ex = np.gradient(-V)
        
        return Ex, Ey
        

    def plot_vector_electric_field(self, Ex, Ey):
        """
        Plots the electric field as a vector field.
        
        Parameters
        ----------
        Ex, Ey : NumPy array
            The x and y components of the electric field.

        Returns
        -------
        None.

        """     
        
        
        xE, yE = np.meshgrid(np.arange(0, Ex.shape[0]), np.arange(0, Ex.shape[1])) # quiver needs to take in meshgrid of the electric field components
        plt.quiver(xE, yE, Ex, Ey, color = 'blue', units = 'width') # plt.quiver plots a 2D arrow field
        plt.title('Vector Electric Field')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.show()


    def plot_mag_electric_field(self, Ex, Ey):
        """
        Plots the magnitude of the electric field as a heat map.

        Parameters
        ----------
        Ex, Ey : NumPy array
            The x and y components of the electric field.

        Returns
        -------
        None.

        """
        

        E_mag = np.sqrt(Ex **2 + Ey **2)
        plt.contourf(self.xv, self.yv, E_mag, 50) # keep same number of contour lines as V field    
        plt.axis('scaled')
        plt.colorbar(label = 'Electric Field Magnitude (V/m)')
        plt.title('Electric Field Magnitude')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.show()
        
        
    def solve_and_plot(self):
        """
        The function that calls out the calculate and plot functions all in one function.
        
        Parameters
        ----------
        None.

        Returns
        -------
        V : NumPy array
            The electric potential at all points on the discrete grid.
        x_max/y_max : integers
            The end x and y coordinates of the grid.

        """
        
        V = self.finite_difference()
        self.plot_potential_field(V)
        self.get_potential_field_value(V)
        Ex, Ey = self.calculate_electric_field(V)
        self.plot_vector_electric_field(Ex, Ey)
        self.plot_mag_electric_field(Ex, Ey)
        #Ex, Ey = self.interpolation_func(V) - didn't get to finishing this off
        
        return V, self.x_max, self.y_max


# activating the class


if __name__ == "__main__":
    solver = LaplaceEquationSolver()
    solver.solve_and_plot()
