# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:45:58 2019

@author: BlackOp
"""

# This the beginning of my program #

import math as mt                     # Import the necessary libraries #

def ArcTanTaylor(x,N):
    """ Function ArcTanTaylor(x,N) will evaluate the inverse tangent function
              of x by testing the convergence of Taylor Series about x.
    """
    tex_sum = 0
    if abs(x) <= 1:                                               # The set of conditionals checks #
        for n in range (0,N):                                     # the ranges in which x belongs. #
            tex_sum += ((-1)**n)*(x**(2*n + 1))/(2*n + 1)                                          
        return tex_sum                                            # For each range , different functions # 
    elif x > 1:                                                   # are set to  approximate the  inverse # 
        for n in range (0,N):                                     # tangent function of a given value x. #
            tex_sum += - ((-1)**n)*(1/x)**(2*n + 1)/(2*n + 1)     
        return tex_sum + (mt.pi)/2                                # The Taylor series expansion of the inverse #
    elif x < -1 :                                                 # tangent function is set to compute as many #
        for n in range (0,N):                                     # terms the user instructs it to by  chosing #
            tex_sum += - ((-1)**n)*(1/x)**(2*n + 1)/(2*n + 1)     # an integer number for the value or N terms.#
        return tex_sum - (mt.pi)/2                                
    
print('Enter A to approximate the inverse tangent of x using the Taylor series expansion.\n')
print('Enter B to compare the approximation with that of the built in inverse tangent function.\n')
print('Enter Q to quit.')

UserInput = '0'
while UserInput != 'Q':                                                         
    UserInput = input('Chose an operation to perform, "A", "B", or "Q":')       # Allows user to choose operation. #
    print('Operation ' + UserInput + ' was chosen.\n')

    if UserInput =='A':
        print('Taylor Series expansion for a chosen value x will follow.')
        Input_x = input('Enter a value for x (floating point number):')          # Allows user to choose any value for x. #
        x = float(Input_x)
        Input_N = input('Enter a value for N (positive integer number):')       # Allows user to choose  positive integer values for N. #
        N = int(Input_N)
        print('The approximation of the arc-tangent of ' + str(x) + ' is :' ,ArcTanTaylor(x,N))
        print('Below, the precision of the expansion for known points will be tested.\n') 
        print('The arc-tan of 1/sqrt(3) is: ' ,ArcTanTaylor(1/mt.sqrt(3),10) , ' which is very close to Pi/6 =' ,mt.atan(1/mt.sqrt(3)))
        print('Difference from real value is: ' ,abs(ArcTanTaylor(1/mt.sqrt(3),10) - mt.atan(1/mt.sqrt(3))),'.\n')
        print('The arc-tan of 1 is: ' ,ArcTanTaylor(1,10) , ' which is close to Pi/4 =' ,mt.atan(1))
        print('Difference from real value is: ' ,abs(ArcTanTaylor(1,10) - mt.atan(1) ),'.\n')
        print('The difference for x=1 is relatively bigger than any other x value as the series converges to that point.')
        
        
    elif UserInput =='B':
        print('The approximation of the arc-tangent of x will be compared to that of the built in arc-tangent function.')
        print('You must chose a value for x in the range [-2,2].')              
        x = 0
        Input_x = input('Enter a value for x (floating point number):' )
        x = float(Input_x)
        if abs(x) <= 2:                                                                    # For x in [-2,2] the approximation #
            Input_N = input('Enter a value for N (positive integer number):')              # of arc-tan with an expansion of N #
            N = int(Input_N)                                                               # terms is compared to the result a #
            print('The computed arc-tangent of ' + str(x) + ' is :' ,ArcTanTaylor(x,N))    # built in arc-tan function  gives, #
            print('The result from the built in function is :' , mt.atan(x))               #   and the difference is printed.  #
            difference = ArcTanTaylor(x,N) - mt.atan(x)
            print('The difference between them is:' , abs(difference))
            data = [['Input x', 'Approximation', 'Computed result', 'Difference'],
                    [-2.0, ArcTanTaylor(-2.0,N), mt.atan(-2.0), abs(ArcTanTaylor(-2.0,N) - mt.atan(-2.0))],
                    [-1.8, ArcTanTaylor(-1.8,N), mt.atan(-1.8), abs(ArcTanTaylor(-1.8,N) - mt.atan(-1.8))],
                    [-1.6, ArcTanTaylor(-1.6,N), mt.atan(-1.6), abs(ArcTanTaylor(-1.6,N) - mt.atan(-1.6))],
                    [-1.4, ArcTanTaylor(-1.4,N), mt.atan(-1.4), abs(ArcTanTaylor(-1.4,N) - mt.atan(-1.4))],
                    [-1.2, ArcTanTaylor(-1.2,N), mt.atan(-1.2), abs(ArcTanTaylor(-1.2,N) - mt.atan(-1.2))],
                    [-1.0, ArcTanTaylor(-1.0,N), mt.atan(-1.0), abs(ArcTanTaylor(-1.0,N) - mt.atan(-1.0))],
                    [-0.8, ArcTanTaylor(-0.8,N), mt.atan(-0.8), abs(ArcTanTaylor(-0.8,N) - mt.atan(-0.8))],  
                    [-0.6, ArcTanTaylor(-0.6,N), mt.atan(-0.6), abs(ArcTanTaylor(-0.6,N) - mt.atan(-0.6))],  
                    [-0.4, ArcTanTaylor(-0.4,N), mt.atan(-0.4), abs(ArcTanTaylor(-0.4,N) - mt.atan(-0.4))],  # This list contains a set of four columns. #
                    [-0.2, ArcTanTaylor(-0.2,N), mt.atan(-0.2), abs(ArcTanTaylor(-0.2,N) - mt.atan(-0.2))],  # Each column displays the x values between #
                    [+0.0, ArcTanTaylor(+0.0,N), mt.atan(+0.0), abs(ArcTanTaylor(+0.0,N) - mt.atan(+0.0))],  # [-2,2] in steps of 0.2,the approximation, #
                    [+0.2, ArcTanTaylor(+0.2,N), mt.atan(+0.2), abs(ArcTanTaylor(+0.2,N) - mt.atan(+0.2))],  # the built  in  function result  and  the  #
                    [+0.4, ArcTanTaylor(+0.4,N), mt.atan(+0.4), abs(ArcTanTaylor(+0.4,N) - mt.atan(+0.4))],  #     absolute value of their difference.   #
                    [+0.6, ArcTanTaylor(+0.6,N), mt.atan(+0.6), abs(ArcTanTaylor(+0.6,N) - mt.atan(+0.6))],  
                    [+0.8, ArcTanTaylor(+0.8,N), mt.atan(+0.8), abs(ArcTanTaylor(+0.8,N) - mt.atan(+0.8))],
                    [+1.0, ArcTanTaylor(+1.0,N), mt.atan(+1.0), abs(ArcTanTaylor(+1.0,N) - mt.atan(+1.0))],
                    [+1.2, ArcTanTaylor(+1.2,N), mt.atan(+1.2), abs(ArcTanTaylor(+1.2,N) - mt.atan(+1.2))],
                    [+1.4, ArcTanTaylor(+1.4,N), mt.atan(+1.4), abs(ArcTanTaylor(+1.4,N) - mt.atan(+1.4))],
                    [+1.6, ArcTanTaylor(+1.6,N), mt.atan(+1.6), abs(ArcTanTaylor(+1.6,N) - mt.atan(+1.6))],
                    [+1.8, ArcTanTaylor(+1.8,N), mt.atan(+1.8), abs(ArcTanTaylor(+1.8,N) - mt.atan(+1.8))],
                    [+2.0, ArcTanTaylor(+2.0,N), mt.atan(+2.0), abs(ArcTanTaylor(+2.0,N) - mt.atan(+2.0))]]

            dash = '-'*80

            for i in range(len(data)):                # Setting a loop which allows the table to print by recalling the list made previously. #
                if i == 0:
                    print(dash)
                    print('{:<8s}{:>20s}{:>20s}{:>20s}'.format(data[i][0],data[i][1],data[i][2],data[i][3]))  # Headings of the table are printed. #
                    print(dash)
                else:
                    print('{:<8f}{:>20f}{:^20f}{:>20f}'.format(data[i][0],data[i][1],data[i][2],data[i][3]))  # All table contents are printed. #
            print(dash)        
            print('Table above shows the difference in the approximation for ' + str(N) + ' terms of the Taylor expansion.')
        
        
        elif abs(x) > 2:                                                        # If the user input x is not #
            print('Chosen value of x was not in the interval [-2,2].')          # in [-2,2] error will arise.#
            
            
            
                
    elif UserInput != 'Q':                      # Any value other than A,B,Q results in error. #
        print('This is not a valid choice')

print(' You have chose to quit this program. --Goodbye-- ')