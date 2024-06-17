"""
@author: Dharun Venkat
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from LaplaceEquationSolver import LaplaceEquationSolver

LaplaceEquationSolver = LaplaceEquationSolver()
V = LaplaceEquationSolver.solve_and_plot() # or change to numerical_solution

difference_list = []

for i in range(0, 51): # finding the fractional difference for point and lines
    
    difference = np.abs(V[50 + i, 50] - V[50 - i, 50])
    fractional_difference = difference / V[50 + i, 50]
    difference_list.append(fractional_difference)


#for i in range(5, 46): # finding the fractional difference radially for a circle
#    
#    right_pot_value = V[50, 50 + i]
#    left_pot_value = V[50, 50 - i]
#    top_pot_value = V[50 + i, 50]
#    bottom_pot_value = V[50 - i, 50]
#    
#    max_value = max(right_pot_value, left_pot_value, top_pot_value, bottom_pot_value)
#    min_value = min(right_pot_value, left_pot_value, top_pot_value, bottom_pot_value)
#    
#    difference = abs(max_value - min_value)
#    fractional_difference = difference / max_value
#    difference_list.append(fractional_difference)


# plotting fractional difference plot

x_values = np.arange(51) # change to 41 for circle

plt.plot(x_values, difference_list, marker='o', linestyle='-')
plt.title('Fractional Difference Against Distance From Potential')
plt.xlabel('Distance Away From Potential')
plt.ylabel('Fractional Difference')
plt.grid(True)
plt.show()


#def V_anal1(x, y, L, H, V_0, n_max): # the retangular plate analytical solution
#    sum_result = 0
#    for n in range(1, n_max + 1, 2):  # summing over odd n only
#        term = (4 / np.sinh(n * np.pi * H / L)) * (1 / (n * np.pi)) * np.sin(n * np.pi * x / L) * np.sinh(n * np.pi * (H - y) / L)
#        sum_result += term
#    return  sum_result


#def V_anal2(x, y, a, b, V_0, n_max): # the rectangular pipe analytical solution
#    sum_result = 0
#    for n in range(1, n_max + 1, 2):  # summing over odd n only
#        term = (1/n) * (np.cosh(n * np.pi * (x - b) / a) / np.cosh(n * np.pi * b / a)) * np.sin(n * np.pi * y / a)
#        sum_result += term
#    return (4 * V_0 / np.pi) * sum_result


# parameters of above

L = 100
H = 100
a = 100
b = 50
V_0 = 1
n_max = 150 # choosing an appropriately large n value

# generate 101 x 101 grid points of even spacing - 101 because 0 points are counted
x_vals = np.linspace(0, L, 101)
y_vals = np.linspace(0, H, 101)
X, Y = np.meshgrid(x_vals, y_vals)

# evaluate the function for each point in the grid

#analytical_solution = np.zeros((101, 101)) # initalise 0 matrix
#for i in range(101):
#    for j in range(101):
#        analytical_solution[i, j] = V_anal1(X[i, j], Y[i, j], L, H, V_0, n_max)


# Calculate fractional difference, skipping over points where analytical solution is zero

#fractional_difference = np.zeros((99, 99))
#for i in range(99):
#    for j in range(99):
#        if analytical_solution[i + 1, j + 1] != 0:
#            fractional_difference[i, j] = np.abs((numerical_solution[i + 1, j + 1] - analytical_solution[i + 1, j + 1]) / analytical_solution[i + 1, j + 1])
#        else:
#            fractional_difference[i, j] = np.nan

# Plotting
#plt.figure()
#plt.imshow(fractional_difference, cmap='jet', origin='lower')
#plt.colorbar(label='Fractional Difference')
#plt.title('Fractional Difference Between Numerical and Analytical Solutions')
#plt.xlabel('X')
#plt.ylabel('Y')
#plt.show()

#plt.figure(figsize=(8, 6))
#plt.contourf(X, Y, analytical_solution, levels=50, cmap='viridis')
#plt.colorbar(label='Voltage (V)')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Contour Plot of V_anal1')
#plt.grid(True)
#plt.show()




        
