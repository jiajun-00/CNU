# Define all functions for Task 3 here.
import numpy as np


def counter(original_f):
    '''
    Decorator to count the number of times original_f is evaluated.
    '''
    # Define our new, decorated function, with added features
    def decorated_f(x):
        # Increment the number of .evals (depending on whether x is a number or an array)
        try:
            l = len(x)
        except:
            l = 1
        decorated_f.evals += l
        
        # Still return the result of the original function
        return original_f(x)
    
    # Initialise the number of evaluations, store it in a new .evals attribute
    decorated_f.evals = 0
    
    # Return the new, decorated function, which now has a .evals value
    return decorated_f


def define_Ka(a):
    '''
    Function which defines and returns a function Ka for a given parameter a.
    The returned function Ka is decorated with counter().
    '''
    @counter
    def Ka(x):
        return 1 / np.sqrt(a - np.cos(x))
    
    # Returns the function itself as an object
    return Ka

def quadrature_trapz(f, xk, wk, a, b,last):
    '''
    Approximates the integral of f over [a, b],
    using the quadrature rule with weights wk
    and nodes xk and record the value calculated in the function.
    
    Input:
    f (function): function to integrate (as a Python function object)
    xk (Numpy array): vector containing all nodes
    wk (Numpy array): vector containing all weights
    a (float): left boundary of the interval
    b (float): right boundary of the interval
    last(float): the value of last used function
    
    Returns:
    return_List (list): the approximate value of the integral
        of f over [a, b], using the quadrature rule and the record function value
    '''
    # Define the shifted and scaled nodes
    yk = (b-a)*(xk+1)/2 + a
    
    # record the new calculated value to use less funtion 
    record_new = f(yk[1])
    # Compute the weighted sum
    I_approx = (b - a)/2 * (wk[0] * last + wk[1] * record_new)
    
    return_list = [I_approx, record_new]
    return return_list

def composite_trapz(f, a, b, M):
    '''
    Returns the approximation of the integral of f
    over [a, b], using the composite trapezoid rule
    with M-1 equal-width partitions.
    '''
    # Find each sub-interval
    bounds = np.linspace(a, b, M)
    
    # Define weights and nodes for trapezoid rule
    xk = np.array([-1., 1.])
    wk = np.array([1., 1.])
    
    # Loop to compute each small integral
    I_approx = 0
    # use last to record the last calculated value to use less function as possible
    last = f(0)
    for i in range(M-1):
        # count each interval and add them up
        current = quadrature_trapz(f, xk, wk, bounds[i], bounds[i+1],last)
        I_approx += current[0]
        # update the new calculated fuction value
        last = current[1]
    
    return I_approx

def Kintegral_trapz(Ka, n):
    '''
    Get the approximate value of kintegral integration from -pi to pi
    param: Ka the input kintegral function
    param: n the number of using Kintegral fuction 
    return: the approximate inetragtion of kintegral function from -pi to pi
    '''
    # for special case n=1, we can only count the midpoint
    if n == 1:
        return Ka(np.pi/2)*2*np.pi
    # else we use trapz and double the interval from 0 to pi, then we will get the whole integration by symmetric
    else:
        return composite_trapz(Ka, 0, np.pi, n)*2



# Other method of Kintegral
def composite_plus(f, a, b, M):
    '''
    Returns the approximation of the integral of f
    over [a, b], using the composite trapezoid rule
    with M-1 equal-width partitions.
    '''
    # Find each sub-interval
    bounds = np.linspace(a, b, M)
    
    # Define weights and nodes for trapezoid rule
    xk_trapz = np.array([-1., 1.])
    wk_trapz = np.array([1., 1.])

    # Define weights and nodes for simpson rule
    xk_simpson = np.array([-1., 0., 1.])
    wk_simpson = np.array([1/3, 4/3, 1/3])

    # Loop to compute each small integral
    I_approx = 0
    record = [0,f(a)]
    start = 0
    # if the number of equal-width partitions is even, we should use trapz once in the middle partition, 
    # since simposn can only be used in odd number of equal-width partitions
    while start != M - 1:
        if M % 2 == 0 and start == (M//4)*2:
            # use trapz rule for middle partition
            current = quadrature_plus(f, xk_trapz, wk_trapz, bounds[start], bounds[start+1],record[0],record[1])
            record = current[1:]
            I_approx += current[0]
            start = start + 1
        else:
            # use simpson for all other partitions
            current = quadrature_plus(f, xk_simpson, wk_simpson, bounds[start], bounds[start + 2],record[0],record[1]) 
            record = current[1:]
            I_approx += current[0]
            start = start + 2   
    return I_approx

def quadrature_plus(f, xk, wk, a, b, record_1, record_2):
    '''
    Approximates the integral of f over [a, b],
    using the quadrature rule with weights wk
    and nodes xk.
    
    Input:
    f (function): function to integrate (as a Python function object)
    xk (Numpy array): vector containing all nodes
    wk (Numpy array): vector containing all weights
    a (float): left boundary of the interval
    b (float): right boundary of the interval
    record_1: record the first value calculated before
    record_2: record the second value caliculated before
    
    Returns:
    return_List (list): the approximate value of the integral
        of f over [a, b], using the quadrature rule and the record function values
    '''
    # Define the shifted and scaled nodes
    yk = (b-a)*(xk+1)/2 + a
    record_newlist = [0,0]
    # Compute the weighted sum
    # when wk is 2 showing use the trapz rule
    if len(wk) == 2:
        # trapz only use two points
        record_newlist = [0 , f(yk[1])]
        I_approx = ((b - a)/2) * (record_2 * wk[0] + record_newlist[1] * wk[1])
    # when wk is 3 showing use the simpson rule    
    elif len(wk) == 3:
        # simpson use three points
        record_newlist = [f(yk[1]) , f(yk[2])]
        I_approx = ((b - a)/2) * (record_2 * wk[0] + record_newlist[0] * wk[1] + record_newlist[1] * wk[2]) 
    # add the approximate integration and record list into one list      
    return_list = [I_approx] + record_newlist
    
    return return_list 

def Kintegral_simpson_trapz(Ka, n):
    '''
    Get the approximate value of kintegral integration from -pi to pi
    param: Ka the input kintegral function
    param: n the number of using Kintegral fuction 
    return: the approximate inetragtion of kintegral function from -pi to pi
    '''
    # for special case n=1, we can only count the midpoint
    if n == 1:
        return Ka(np.pi/2)*2*np.pi
    # else we use trapz and double the interval from 0 to pi, then we will get the whole integration by symmetric
    else:
        return composite_plus(Ka, 0, np.pi, n)*2
    
def Kintegral(Ka,n):
    '''
    Get the approximate value of kintegral integration using midpoint rule from -pi to pi
    param: Ka the input kintegral function
    param: n the number of using Kintegral fuction 
    return: the approximate inetragtion of kintegral function from -pi to pi
    '''
    #find all the midpoint in the interval of [0,pi] and the end point should not be seem as midpoint
    xk = np.linspace((np.pi /(2*n)),np.pi+ (np.pi /(2*n)), n+1)
    # calculate each value with Ka function(exclude the last point, since it is not in the interval[0,pi])
    # times two to get the whole inetragtion from -pi to pi
    yk = (np.pi)/n * Ka(xk[:-1]) * 2 
    # sum them up to get the the approximate integration
    return np.sum(yk)

def quadrature_gussian(f, xk, wk, a, b):
    '''
    Approximates the integral of f over [a, b],
    using the quadrature rule with weights wk
    and nodes xk.
    
    Input:
    f (function): function to integrate (as a Python function object)
    xk (Numpy array): vector containing all nodes
    wk (Numpy array): vector containing all weights
    a (float): left boundary of the interval
    b (float): right boundary of the interval
    
    Returns:
    I_approx (float): the approximate value of the integral
        of f over [a, b], using the quadrature rule.
    '''
    # Define the shifted and scaled nodes
    yk = (b-a)*(xk+1)/2 + a
    
    # Compute the weighted sum
    I_approx = (b - a)/2 * np.sum(wk * f(yk))
    
    return I_approx

def Kintegral_gussian(Ka,n):
    '''
    Use Gussian Legendre rule to calculate the integration
    :param Ka: the given function
    :param n: number of using the function
    :return: the approxiamte integration
    '''
    return quadrature_gussian(Ka, np.polynomial.legendre.leggauss(n)[0], np.polynomial.legendre.leggauss(n)[1], 0, np.pi)*2