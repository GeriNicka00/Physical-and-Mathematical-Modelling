# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 11:26:05 2020

@author: BlackOp
"""

# -*- coding: utf-8 -*-
"""
FRESNEL DIFFRACTION SIMULATION
"""
# Import required libraries #

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# IMPORTANT NOTES FOR THE USER #
print('Please input (SI) units only.\n')
print('TO SEE NEAR-FIELD DIFFRACTION THE APERTURE AND SCREEN SIZE, THE DISTANCE BETWEEN THEM AND N MUST BE ADJUSTED ACCORDINGLY THROUGH MANY TRIALS\n')
print('NOTE THAT IF YOU CHOOSE N LARGER 100 FOR PART C TOO MUCH TIME IS NEEDED TO PLOT RESULTS. DONT TRY PLOTTING PARTS A,B,C TOGETHER WITH TOO LARGE N\n')
print('FOR PART A AND B LARGE N VALUES CAN BE USED TO ACHIEVE FRESNEL DIFFRACTION RAPIDLY\n')
print('INVESTIGATE PART C INDEPEDENTLY AND MAKE SURE CODE AND APERTURE SHAPE FILES ARE IN THE SAME DESTINATION')

# Define variables and allow user to control  number of them to investigate diffraction patterns. #
E0 = 1
λ  = 1e-6
k  = 2*(np.pi)/λ
# Choosing aperture width #
Input_x1_p = input('Enter a negative value for the origin of the square x1_p and the aperture width will be double its absolute value:') 
x1_p = float(Input_x1_p)
x2_p = -float(Input_x1_p)
y1_p = float(Input_x1_p)
y2_p = -float(Input_x1_p)
# Choosing screen width #
Input_x1 = input('Enter a negative value for the origin of the square x1 and the screen width will be double its absolute value:') 
x1 = float(Input_x1)
x2 = -float(Input_x1)
y1 = float(Input_x1)
y2 = -float(Input_x1)

Input_z = input('Enter a value for the distance between the aperture and the screen in meters (values as in the script):')
z = float(Input_z)

Input_N = input('Enter a value for the number of subintervals N to perform Simpson rule (EVEN INTEGER):')
N = float(Input_N)

x_p = np.linspace(x1_p,x2_p,N+1)
y_p = np.linspace(y1_p,y2_p,N+1)

print("Type 'Slit.png' or 'TwoSlit.png' or 'TwoCircles.png' or '3Circles.png' or 'Triangle.png' for your aperture shape in part C.")
Input_aperture = input('Chose your aperture shape for part C :')
aperture = Input_aperture

# Define functions #

def Simpson(f,xi,x1_p,x2_p,N):
    '''
    Approximate the integral of f(x) using Simpson's Rule.
    
    Parameters:
    -----------
    -> f: Vectorized function of a single variable x.
    -> x1,x2: Limits of integration in aperture coordinates.
    -> N: Number of subintervals within the limits of integration [x1_p,x2_p] 
    
    Calculates:
    -----------
    -> Approximation of the integral of f(x) in 1-D from x1 to x2 using
    Simpson's rule with N subintervals of equal length.
    
    '''
    dx_p = (x2_p-x1_p)/N
    x_p = np.linspace(x1_p,x2_p,N+1)
    y = f(xi, x_p)
    S = dx_p/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2]) 
    return S

def Integration2D_Simpson(f, xi, yj, x1_p, x2_p, y1_p, y2_p, N, anyshape):
    """
    Approximates the integral of f(x,y) using Simpson's Rule
    twice in the x and y coordinates for 2-D.
    
    """
    dx_p =(x2_p-x1_p)/N 
    dy_p =(y2_p-y1_p)/N
    
    E_x_p = np.zeros((int(N+1),int(N+1)), dtype=complex)
    
    for i, x_pi in enumerate(x_p):                                              
        f_yp = f(xi, x_pi, yj, y_p, anyshape[i])                                # Loop assigns values in the x and y axes #
        E_x_p[i] = (dy_p/3)*np.sum(f_yp[0:-1:2]+f_yp[2::2]+4*f_yp[1::2])        # using Simpson's rule twice and  creates # 
                                                                                # sets of (x,y) which produce a 2-D image #
    Eintegral = (dx_p/3)*np.sum(E_x_p[0:-1:2]+E_x_p[2::2]+4*E_x_p[1::2])
        
    return Eintegral

def E2_D_AnyShape(xi, x_p, yj, y_p, anyshape):                                              # Allows the investigation of diffraction #
    return anyshape*(k*E0/(2*np.pi*z))*np.exp(((0+1j)*k/(2*z))*((xi-x_p)**2+(yj-y_p)**2))   # patterns in 2D using any aperture shape #
                                                                                             
def E1_D(xi,x_p):                                                               
    Ex = ((k*E0)/(2*(np.pi)*z))*np.exp(((0+1j)*k/(2*z))*(xi-x_p)**2) 
    return Ex

def Intensity(xi, x1_p, x2_p, N):                                                
    E1D = Simpson(E1_D, xi, x1_p, x2_p, N)                                      # Evaluates intensity in one dimension #               
    Intensity = np.power((np.absolute(E1D)),2)                                  # recalling the defined Simpson method #
    return Intensity

def IntensityPlot():
    x = np.linspace(-0.005,0.005,1000)
    Intensity_vals = np.zeros(len(x), dtype=complex)                                
    for i in np.arange(len(x)):                                                     
        Intensity_vals[i] = Intensity(x[i], x1_p, x2_p, N)                          
    plt.title(" Diffraction in 1D for square aperture. ")                           # Arrays of distance and  intensity values  are # 
    plt.xlabel(" Screen coordinate (m) ")                                           # initiated and plotted  against each other  to #
    plt.ylabel(" Relative Intensity ")                                              # visualise one dimensional Fresnel diffraction #
    plt.plot(x, Intensity_vals)
    plt.show()
    return IntensityPlot 

def Intensity2_DPlot():
    x = np.linspace(x1,x2,1001)
    a = np.zeros(len(x), dtype=complex)                                             
    y = np.linspace(y1,y2,1001)                                           
    b = np.zeros(len(x), dtype=complex)                                             # Two dimensional array with dimensions equal to the #
    Intensity2_D_Vals = np.zeros([len(x),len(y)])                                   # length of the x & y arrays stores intensity values #
    for i in np.arange(len(x)):                                                     # needed to plot Fresnel diffraction in 2 dimensions #
        a[i] = Intensity(x[i], x1_p, x2_p, N)
        b[i] = a[i]
    for i in range(1001):
        for j in range(1001):
            Intensity2_D_Vals[i,j] = a[j]*b[i]
    plt.title("Diffraction in 2D for square aperture.")        
    plt.imshow(Intensity2_D_Vals, cmap=plt.get_cmap('jet_r'), origin='lower')
    plt.xlabel('Horizontal Pixel Number')
    plt.ylabel('Vertical Pixel Number')
    plt.show()
    return Intensity2_DPlot

# Menu allows user to invetigate Fresnel Diffraction in 1-D and 2-D #

print('Enter A to investigate Fresnel Diffraction in 1D. \n')
print('Enter B to investigate Fresnel Diffraction in 2D.\n')
print('Enter C to investigate Fresnel Diffraction in 2D for any aperture shape.\n')
print('Enter E to quit.')

UserInput = '0'
while UserInput != 'E': 
    UserInput = input('Chose an operation to perform, "A", "B", "C" or "E":')   
    if UserInput == "A":                                                        # Choosing A plots intensity in 1D #
        IntensityPlot() 
                                                      
    elif UserInput =='B':                                                       # Choosing B plots intensity in 2D #
        Intensity2_DPlot()
        
    elif UserInput =='C':
        anyshape = np.array(Image.open(aperture).convert('L').resize((int(N+1),int(N+1))))          # Plots Fresnel Diffraction pattern in 2D for #
        xi = np.linspace(x1,x2,51)                                                                  # any arbitrary aperture shape using the same #
        yj = np.linspace(x1,x2,51)                                                                  # concept implemented for the square aperture #
        Eintegral = np.zeros((len(xi),len(yj)))
        for j in np.arange(len(yj)):
            for i in np.arange(len(xi)):
                Eintegral[i][j] = np.abs(Integration2D_Simpson(E2_D_AnyShape, xi[i], yj[j], x1_p, x2_p, y1_p, y2_p, N, anyshape))**2
        plt.title(' Diffraction pattern of ' +str(aperture)+ ' in 2D ') 
        plt.imshow(Eintegral, cmap=plt.get_cmap('jet_r'), origin='lower')
        plt.xlabel('Horizontal Pixel Number')
        plt.ylabel('Vertical Pixel Number')
        plt.show()    
        
    elif UserInput != 'E':                       # Choices other than A,B,C or E result in error #
        print('This is not a valid choice')

print(' You have chose to quit. --Goodbye-- ')   # Choosing option E terminates the program #





    
    







