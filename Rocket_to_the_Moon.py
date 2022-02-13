# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:18:18 2020

@author: BlackOp
"""
# Import libraries #

import matplotlib.pyplot as plt
import numpy as np

# Input functions #

def time_choice():             # Gives user option to choose orbital time period. #
    print("\nSimulation time")
    user_input= input("Choose to: \n(A) Orbit for 19,000 seconds (default). \n(B) Orbit with custom time period. \nCHOICE: ")
    while user_input!= 'A' or user_input!= 'B' or user_input!= 'C':
        if user_input == "A":
            tmax= 19000
            break
        elif user_input == "B":
            user_inputt = input("Enter the simulating time (s): ")
            tmax= int(user_inputt)
            break
    return tmax

def xy_choice():               # Gives user option to choose x and y coordinates. #
    print("\n[X AND Y COORDINATES]")
    user_input= input("Choose to: \n(A) Keep the default x and y coordinates (15000 km) and (0 km). \n(B) Change them. \nCHOICE: ")
    while user_input!= 'A' or user_input!= 'B' or user_input!= 'C':
        if user_input == "A":
            print("\nStarting coordinates chosen as 15000 km for x and 0 km for y.")
            x1= 15000000
            y1= 0
            break
        elif user_input == "B":
            print("\nCustom starting coordinates chosen.")
            user_inputx = input("Starting x coordinate (km): ")
            user_inputy = input("Starting y coordinate (km): ")
            x1= float(user_inputx)*1000
            y1= float(user_inputy)*1000
            break
    return x1,y1

def v_choice():                # Gives user option to choose Vx and Vy components. #
    print("\n[Vx AND Vy VALUES]")
    user_input= input("Choose to: \n(A) Keep default Vx and Vy values (0 m/s) and (5153 m/s). \n(B) Change them. \nCHOICE: ")
    while user_input!= 'A' or user_input!= 'B' or user_input!= 'C':
        if user_input == "A":
            print("\nStarting coordinates chosen to be 3,600 km for x and y.")
            Vx= 0
            Vy= 5153
            break
        elif user_input == "B":
            print("\nCustom starting velocities chosen.")
            print('Default: 0m (Vx) and 5153.2m (Vy)')
            user_inputVx = input("Starting Vx value (m/s): ")
            user_inputVy = input("Starting Vy value (m/s): ")
            Vx= float(user_inputVx)
            Vy= - float(user_inputVy)
            break
    return Vx,Vy

# Definining equation to be solved with Runge-Kutta. #

def dx_dt(Vx):
    dx_dt= Vx
    return dx_dt


def dy_dt(Vy):
    dy_dt= Vy
    return dy_dt


def dVx_dt(x,y):
    r=(float(x)**2 + float(y)**2)**(1/2)
    rm=((Rm-x)**2 + y**2)**(1/2)
    dVx_dt= ((-G * M * x)/ abs(r)**3)+ ((G * M_moon * (Rm-x))/ abs(rm)**3)
    return dVx_dt


def dVy_dt(x,y):
    r= (float(x)**2 + float(y)**2)**(1/2)
    rm= ((Rm-x)**2 + y**2)**(1/2)
    dVy_dt= ((-G * M * y) / abs(r)**3)- ((G * M_moon * (y))/ abs(rm)**3)
    return dVy_dt  

def Runge_Kutta(t,x,y,Vx,Vy):
    xvals=[]
    yvals=[]
    Vxvals=[]                       # Creating empty arrays to store values. #
    Vyvals=[]
    tvals=[]
    
    xvals_set=[]
    yvals_set=[]
    Vxvals_set=[]
    Vyvals_set=[]
    
    KE = []
    PE = []
    TE = []
    for i in range(0,tmax, h):      # Defining Runge-Kutta coefficients. #
    
        k1x= dx_dt(Vx)
        k1y= dy_dt(Vy)
        k1Vx=dVx_dt(x,y)
        k1Vy=dVy_dt(x,y)
        
        k2x= dx_dt(Vx + ((h*k1Vx)/2))
        k2y= dy_dt((Vy + ((h* k1Vy)/2)))
        k2Vx= dVx_dt((x + ((h*k1x)/2)),(y + ((h*k1y)/2)))
        k2Vy= dVy_dt((x + ((h*k1x)/2)),( y+ ((h*k1y)/2)))
        
        k3x= dx_dt((Vx + ((h*k2Vx)/2)))
        k3y=dy_dt((Vy + ((h*k2Vy)/2)))
        k3Vx=dVx_dt((x + ((h*k2x)/2)),(y + ((h*k2y)/2)))
        k3Vy=dVy_dt((x + ((h*k2x)/2)),(y + ((h*k2y)/2)))
        
        k4x= dx_dt((Vx + h* k3Vx))
        k4y=dy_dt((Vy + h*k3Vy))
        k4Vx=dVx_dt((x + h* k3x),(y + h* k3y))
        k4Vy=dVy_dt((x + h* k3x),(y + h* k3y))
            
        x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)     # Calculating position coordinates. #
        xvals.append(x)
        xt= round(x)
        xvals_set.append(xt)
        
        y= y + (h/6)*(k1y + 2*k2y + 2*k3y + k4y)
        yvals.append(y)
        yt= round(y)
        yvals_set.append(yt)
        
        Vx= Vx + (h/6)*(k1Vx + 2*k2Vx + 2*k3Vx + k4Vx)  # Calculating velocity coordinates. #
        Vxvals.append(Vx)
        Vxt= round(Vx)
        Vxvals_set.append(Vxt)
        
        Vy= Vy + (h/6)*(k1Vy + 2*k2Vy + 2*k3Vy + k4Vy)
        Vyvals.append(Vy)
        Vyt = round(Vy)
        Vyvals_set.append(Vyt)
        
        ke = 0.5 * 1 * (Vx**2 + Vy**2)               # Defining energies of the orbit. #
        KE.append(ke)
        pe = -(6.67e-11*5.97e24)/(x**2 + y**2)**0.5
        PE.append(pe)
        te = ke + pe
        TE.append(te)
        t = t + h
        tvals.append(t)

    return  tvals, xvals, yvals, Vxvals, Vyvals, KE, PE, TE, xvals_set, yvals_set, Vxvals_set, Vyvals_set


user_input= "0"
while user_input != "Q":

    user_input= input("MAIN MENU \nChoose to: \n(A) Orbit the earth. \n(B) Launch satellite to the moon. \n(Q) Quit. \nChoice: ")
    
# Option A - To orbit the earth. #
             
    if user_input == 'A': 
        
        G= 6.67e-11
        M=5.97e24
        M_moon= 7.347e22 
        h = 300
        x,y = xy_choice()
        Vx, Vy = v_choice() # Velocities required to escape the earth's gravitational pull. #
        t = 0
        tmax = time_choice()
        
        R=0
        Rm=384400000
        rm=0
        
        plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[1], Runge_Kutta(t,x,y,Vx,Vy)[2]) # Plot trajectory of rocket orbit around Earth. #
        plt.scatter(0,0)
        circle1=plt.Circle((0,0),6300e3, color = 'k')
        plt.gcf().gca().add_artist(circle1)
        plt.xlabel("x Coordinates (m)")
        plt.ylabel("y Coordinates (m)")
        plt.axis('equal')
        plt.show()
        
        r=np.zeros(len(Runge_Kutta(t,x,y,Vx,Vy)[1]))             # Calculating ellipticity. #
        for i in np.arange(len(Runge_Kutta(t,x,y,Vx,Vy)[1])):
            r[i] = np.sqrt(Runge_Kutta(t,x,y,Vx,Vy)[1][i]**2+Runge_Kutta(t,x,y,Vx,Vy)[2][i]**2)
        r_a=np.max(r)   # Apoapsis # 
        r_p=np.min(r)   # Periapsis #
        e = (r_a-r_p)/(r_a+r_p)
        print('Ellipticity is:', e)
                
        input0=0
        input0= input("Would you like to see the change in energy of this orbit? \n(A) No \n(B) Yes \nChoice: ")
        while input0 != 'A' or input0 != 'B':
            if input0 == 'A':
                break
            elif input0 == 'B':  # Plot energies for rocket's orbit around the Earth. #
                plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[5], label=("Kinetic Energy"))
                plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[6], label=("Potential Energy"))
                plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[7], label=("Total Energy"))
                plt.xlabel("Time (s)")
                plt.ylabel("Energy (J)")
                plt.legend()
                plt.show()
                
                break
            else:
                print('Please select "A" or "B": ' )                    

# Option B - Slingshot #
        
    elif user_input == 'B':
        print('You have chosen part (B)')
        print("")
        
        G= 6.67e-11
        M=5.972e24
        M_moon= 7.35e22    # Moon Mass #
        h = 50

        t = 0
        tmax = 900000

        Rm=384400000
        
        user_input_option=0
        user_input_option= input("Chose to: \n(A) Plot one graph \n(B) Compare 4 plots of different positions. \n(C) Compare 4 plots of different velocities \nChoice: ")
        while user_input_option != 0:
            if user_input_option == 'A':

                x,y = -6500e3, 0
                Vx, Vy = 0, 10964     
                
                user_input_energies = 0
                user_input_energies = input('Display the energy changes due to orbit? \n(A) No \n(B) Yes \nChoice: ')
                while user_input_energies != 'A' or user_input_energies != 'B':
                    
                    if user_input_energies == 'A':
                        print('Will not display energy changes.')
                        break
                    
                    elif user_input_energies == 'B':  # Plot energies for slingshot. #
                        plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[5], label=("Kinetic Energy"))
                        plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[6], label=("Potential Energy"))
                        plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[0], Runge_Kutta(t,x,y,Vx,Vy)[7], label=("Total Energy"))
                        plt.xlabel("Time (s)")
                        plt.ylabel("Energy (J)")
                        plt.legend()
                        plt.show()
                        break

                    else:
                        print('Enter a valid choice: ')
                        break
                    
                plt.plot(Runge_Kutta(t,x,y,Vx,Vy)[1], Runge_Kutta(t,x,y,Vx,Vy)[2])  # Plot slingshot trajectory. #
                plt.scatter(0,0, s=90)
                plt.scatter(384400000,0)
                plt.xlabel("x Coordinates (m)")
                plt.ylabel("y Coordinates (m)")
                plt.axis('equal')
                plt.ylim(-1e8,1e8)
                plt.show()  
                
                break
            
            
            elif user_input_option == 'B':    # Plot 4 trajectories of different starting positions. # 
                x1 = float(input("1st x position (m): "))
                y1 = float(input("1st y position (m): "))
                x2 = float(input("2nd x position (m): "))
                y2 = float(input("2nd y position (m): "))
                x3 = float(input("3rd x position (m): "))
                y3 = float(input("3rd y position (m): "))
                x4 = float(input("4th x position (m): "))
                y4 = float(input("4th y position (m): "))                
                
                Vx, Vy = 0, 10964
                plt.plot(Runge_Kutta(t,x1,y1,Vx,Vy)[1], Runge_Kutta(t,x1,y1,Vx,Vy)[2], label=(x1, "m"))
                plt.plot(Runge_Kutta(t,x2,y2,Vx,Vy)[1], Runge_Kutta(t,x2,y2,Vx,Vy)[2], label=(x2, "m"))
                plt.plot(Runge_Kutta(t,x3,y3,Vx,Vy)[1], Runge_Kutta(t,x2,y2,Vx,Vy)[2], label=(x3, "m"))
                plt.plot(Runge_Kutta(t,x4,y4,Vx,Vy)[1], Runge_Kutta(t,x2,y2,Vx,Vy)[2], label=(x4, "m"))
                
         
                plt.scatter(0,0, s=90)
                plt.scatter(384400000,0)
                plt.xlabel("x Coordinates (m)")
                plt.ylabel("y Coordinates (m)")
                plt.axis('equal')
                plt.ylim(-1e8,1e8)
                plt.show()  

                break
            elif user_input_option == 'C':    # Plot 4 trajectories of different starting velocities. #
                Vx1 = float(input("1st Vx value(m/s): "))
                Vy1 = float(input("1st Vy value(m/s): "))
                Vx2 = float(input("2nd Vx value(m/s): "))
                Vy2 = float(input("2nd Vy value(m/s): "))
                Vx3 = float(input("3rd Vx value(m/s): "))
                Vy3 = float(input("3rd Vy value(m/s): "))
                Vx4 = float(input("4th Vx value(m/s): "))
                Vy4 = float(input("4th Vy value(m/s): ")) 
                
                x,y = -6500e3, 0
                
                plt.plot(Runge_Kutta(t,x,y,Vx1,Vy1)[1], Runge_Kutta(t,x,y,Vx1,Vy1)[2], label=(Vy1, "m/s"))
                plt.plot(Runge_Kutta(t,x,y,Vx2,Vy2)[1], Runge_Kutta(t,x,y,Vx2,Vy2)[2], label=(Vy2, "m/s"))
                plt.plot(Runge_Kutta(t,x,y,Vx3,Vy3)[1], Runge_Kutta(t,x,y,Vx3,Vy3)[2], label=(Vy3, "m/s"))
                plt.plot(Runge_Kutta(t,x,y,Vx4,Vy4)[1], Runge_Kutta(t,x,y,Vx4,Vy4)[2], label=(Vy4, "m/s"))
                
                
                plt.scatter(0,0, s=90)
                plt.scatter(384400000,0)
                plt.xlabel("x Coordinates (m)")
                plt.ylabel("y Coordinates (m)")
                plt.axis('equal')
                plt.ylim(-1e8,1e8)
                plt.legend()
                plt.show()  
                break
            
        else:
            print('Please enter a valid choice: ')


    
    elif user_input == 'Q': 
        print('You chose to quit')
        
    else:
        print("Choice invalid.")