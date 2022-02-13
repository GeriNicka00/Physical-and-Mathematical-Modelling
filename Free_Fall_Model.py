# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:30:55 2019

@author: BlackOp
"""

# Simulation of a free falling skydiver. #

# Importing required libraries. #

import numpy as np
import matplotlib.pyplot as plt

# Defining height and vertical velocity functions with the analytical method. #

def Height(time):
    y_an = 1000 - (m/k)*np.log(np.cosh(np.sqrt((k*g)/m)*time))
    return y_an 

def VerticalVelocity(time):
    Vy_an = -np.sqrt((m*g)/k)*np.tanh((np.sqrt((g*k)/m))*time)
    return Vy_an

# Creating the menu to allow user to control the problem parameters. #

print('YOU CAN CONTROL THE PARAMETERS OF THE FREE FALL.')

Input_m = float(input('Enter a positive, non zero value for mass(kg):'))
m = float(Input_m)

Input_Cd = float(input('Enter a positive, non zero value for the drag coefficient(no units):'))
Cd = float(Input_Cd)

Input_A = float(input('Enter a positive non zero value for the cross sectional area of the object(m^2):'))
A = float(Input_A)
         
Input_tf = float(input('Enter a positive value for time of free fall(s):'))
print('\n')

ρ0 = 1.2               # Air density (kg/m^3). #
tf = float(Input_tf)   # Duration of free fall. #
k = (Cd*A*ρ0)/2        # Constant k (kg/m). #
g  = 9.81              # Gravitational acceleration (m/s^2). #
h = 7640               # Scale height (m). #

# Using the two Euler approximations to plot height and vertical velocity. #
    
print('YOU CAN CONTROL THE PARAMETERS OF THE EULER APPROXIMATIONS.')
                                         
Input_Δt = input('Chose length of time increment steps:')
Δt = float(Input_Δt)
    
tsteps = int(tf/Δt)  

# Initializing the required arrays used in the following loops. #

y_eu = np.zeros([tsteps+1])          # Position values for Euler method stored. #
v_eu = np.zeros([tsteps+1])          # Velocity values for Eulear method stored. #

v_mid = np.zeros([tsteps+1])         # Velocity at middle point of each time step stored. #
y_meu = np.zeros([tsteps+1])         # Position values for Modified Euler method stored. #
v_meu = np.zeros([tsteps+1])         # Velocity values for Modified Euler method stored. #

k_y_eu = np.zeros([tsteps+1])        # Constant k w.r.t position used with Euler method. #
vel_den_eu = np.zeros([tsteps+1])    # Velocity values with varying density using Euler method stored. #
y_den_eu = np.zeros([tsteps+1])      # Position values with varying density using Euler method stored. #

k_y_meu = np.zeros([tsteps+1])       # Constant k w.r.t position used with Modified Euler method. # 
vel_den_meu = np.zeros([tsteps+1])   # Velocity values with varying density using Modified Euler method stored. #
y_den_meu = np.zeros([tsteps+1])     # Position values with varying density using Modified Euler method stored. #
                                 
time = np.linspace(0, tf, tsteps+1)  # Array of time values combined with arrays above which results in the desired graphs. #


y_eu[0] = 1000         # Initial height with Euler method. #
v_eu[0] = 0            

v_mid[0] = 0
y_meu[0] = 1000        # Initial height with Modified Euler method. #
v_meu[0] = 0  

vel_den_eu[0] = 0
y_den_eu[0] = 39045    # Initial height for varying density with Euler method. #

vel_den_meu[0] = 0
y_den_meu[0] = 39045   # Initial height for varying density with Modified Euler method. #

for step in range(tsteps): # Loop allows approximation of position and velocity with Euler and Modified Euler methods. #
    # Euler. #
    v_eu[step+1] = v_eu[step] - Δt*(g + (k/m)*np.absolute(v_eu[step])*v_eu[step])
    y_eu[step+1] = y_eu[step] + Δt*v_eu[step]
    # Terminates loop when ground is reached. #
    if y_eu[step+1] <= 0 :
        break
    # Modified Euler. #
    v_mid[step+1] = v_meu[step] - (Δt/2)*(g + (k/m) *np.absolute(v_meu[step])*v_meu[step])
    v_meu[step+1] = v_meu[step] - (Δt)*(g + (k/m) *np.absolute(v_mid[step+1])*v_mid[step+1])
    y_meu[step+1] = y_meu[step] + Δt * v_mid[step+1]
    # Terminates loop when ground is reached. #
    if y_meu[step+1] <= 0 :
        break
    
v1 = np.zeros(tsteps+1)     # Array stores position values obtained from analytical method. #
y1 = np.zeros(tsteps+1)     # Array stores velocity values obtained from analytical method. #

# Loop assigns position and velocity values to each entry in the arrays so that the loop eventually identifies height zero. #  
for j in np.arange(tsteps):
    v1[j] = VerticalVelocity(time[j])
    y1[j] = Height(time[j])
    # Terminates loop when ground is reached for Analytical method. #
    if y1[j] < 0 :
        break

# Introducing constant k dependent on height with Euler method. #
# Loop allows approximation of position and velocity for varying density with Euler Method. #
        
for step in np.arange(tsteps):
    k_y_eu[step+1] = (Cd*ρ0*(np.exp(-y_den_eu[step]/h))*A)/2
    vel_den_eu[step+1] = vel_den_eu[step] - Δt*(g + ((k_y_eu[step+1])/m)*np.absolute(vel_den_eu[step])*vel_den_eu[step]) 
    y_den_eu[step+1] = y_den_eu[step] + Δt*vel_den_eu[step]
    # Terminates loop when ground is reached. #
    if y_den_eu[step+1] <= 0 :
        break

# Introducing constant k dependent on height with Modified Euler method. #
# Loop allows approximation of position and velocity for varying density with Modified Euler Method. #  
        
for step in np.arange(tsteps):    
    k_y_meu[step+1] = (Cd*ρ0*(np.exp(-y_den_meu[step]/h))*A)/2
    vel_den_meu[step+1] = vel_den_meu[step] - Δt*(g + ((k_y_meu[step+1])/m)*np.absolute(vel_den_meu[step])*vel_den_meu[step]) 
    y_den_meu[step+1] = y_den_meu[step] + Δt*vel_den_meu[step]
    # Terminates loop when ground is reached. #
    if y_den_meu[step+1] <= 0 :
        break
    
# Subplot of plots of position versus time for all methods. #
    
plt.subplot(2,2,1)
plt.title('1. Position-Time plots for all methods.')
plt.plot(time, y1, 'r', label="Analytic")
plt.plot(time, y_eu, 'c', label="Euler")
plt.plot(time, y_meu, 'm', label="Modified Euler")
plt.xlabel('Time(s)')
plt.ylabel('Height(m)')
plt.legend(loc="best")
plt.show()

# Subplot of plots of velocity versus time for all methods. #

plt.subplot(2,2,2)
plt.title('2. Velocity-Time plots for all methods.')
plt.plot(time, v1, 'r', label="Analytic")
plt.plot(time, v_eu, 'c', label="Euler")
plt.plot(time, v_meu, 'm', label="Modified Euler")
plt.xlabel('Time(s)')
plt.ylabel('Vertical Velocity(m/s)')
plt.legend(loc="best")
plt.show()

# Subplot of plots of position versus time with varying density for both Euler methods. # 
 
plt.subplot(2,2,3)
plt.title('3. Position-Time with varying density for both Euler methods.')
plt.plot(time, y_den_eu, 'y' , label="Varying Density with Euler")
plt.plot(time, y_den_meu, 'k',linestyle='dotted', label="Varying Density with Modified Euler")
plt.xlabel('Time(s)')
plt.ylabel('Height(m)')
plt.legend(loc="best")
plt.show()

# Subplot of plots of velocity versus time with varying density for both Euler methods. #

plt.subplot(2,2,4)
plt.title('4. Velocity-Time with varying density for both Euler methods.')
plt.plot(time, vel_den_eu, 'y', label="Varying Density with Euler")
plt.plot(time, vel_den_meu, 'k',linestyle='dotted', label="Varying Density with Modified Euler")
plt.xlabel('Time(s)')
plt.ylabel('Vertical Velocity(m/s)')
plt.legend(loc="best")
plt.show()

# Variables below store terminal velocities predicted with each approximation method. #

term_eu = min(v_eu)
term_meu = min(v_meu)
term_den_eu = min(vel_den_eu)
term_den_meu = min(vel_den_meu)

# Printing terminal velocities predicted with each approximation method. #

print("Terminal velocity predicted with Euler method:", term_eu)
print("Terminal velocity predicted with Modified Euler method:", term_meu)
print("Terminal velocity predicted with Euler method in vaarying density:", term_den_eu)
print("Terminal velocity predicted with Modified Euler method in varying density:", term_den_meu)


# END OF SIMULATION CODE #





