# Define your functions for Task 3 here.
import numpy as np

def trajectory(horizontal_speed, initial_height, xmax, nmax):
    """Calculate the trajectory of a projectile.
    
    Parameters
    ----------
    horizontal_speed : float
        The initial horizontal speed of the projectile.
    initial_height : float
        The initial height of the projectile.
    xmax : float
        The maximum x-coordinate of the trajectory.
    nmax : int
        The maximum number of points to calculate.
    
    Returns
    -------
    x : numpy.ndarray
        The x-coordinates of the trajectory.
    y : numpy.ndarray
        The y-coordinates of the trajectory.
    """
    a = 1
    b = 0.25
    g = 980
    kappa = 0.5
    # if the number of points is 1, return the initial value
    if nmax == 1:
        return np.array([0]), np.array([initial_height])
    dx = xmax / (nmax - 1)
    # set the x-coordinates with nmax points
    x = np.linspace(0, xmax, nmax)
    y = np.zeros(nmax)
    # set the initial value
    y[0] = initial_height
    y[1] = - dx ** 2 * g *(1 - np.exp(-(y[0] - a)/b))/(2* horizontal_speed ** 2) + y[0]
    # calculate the trajectory by given formula
    for i in range(2, nmax):
        y[i] = (- g *(1 - np.exp(-(y[i-1]-a)/b))+ \
                (horizontal_speed ** 2) *(2* y[i-1] - y[i-2])/(dx ** 2) +\
                      kappa * horizontal_speed * y[i -2]/(2 * dx))/\
                        ((horizontal_speed ** 2)/ (dx ** 2) + (kappa * horizontal_speed)/ (2 * dx))
    return x, y

def residual(horizontal_speed, initial_height, x_teacup, y_teacup, nmax):
    """
    Calculate the residual between the trajectory and the teacup in y axis with x axis of tea cup.

    Parameters
    ----------
    horizontal_speed : float
        The initial horizontal speed of the projectile.
    initial_height : float
        The initial height of the projectile.
    x_teacup : float
        The x-coordinate of the teacup.
    y_teacup : float
        The y-coordinate of the teacup.
    nmax : int
        The maximum number of points to calculate.

    Returns
    -------
    residual : float
        The residual between the trajectory and the teacup in y axis with x axis of tea cup.
    """
    x, y = trajectory(horizontal_speed, initial_height, x_teacup, nmax)
    # calculate the residual
    return y[-1] - y_teacup

def times_changing_direction(x, y):
    '''
    Count the times of changing direction of the trajectory.

    Parameters
    ----------
    x : numpy.ndarray
        The x-coordinates of the trajectory.
    y : numpy.ndarray
        The y-coordinates of the trajectory.

    Returns
    -------
    times : int
        The times of changing direction of the trajectory.
    '''
    # calculate the deta x
    dx = x[1] - x[0]
    # calculate the derivative of y
    dy = (y[2:] - y[:-2]) / (2 * dx)
    times = 0
    # count the times of changing direction
    for i in range(len(dy) - 1):
        # if the derivative of y changes sign then the ball changes direction
        if dy[i] * dy[i + 1] < 0:
            times += 1
    return times

def find_horizontal_speed(initial_height, x_teacup, y_teacup, nmax, bounce = False):
    '''
    Find the horizontal speed of the projectile.

    Parameters
    ----------
    initial_height : float
        The initial height of the projectile.
    x_teacup : float
        The x-coordinate of the teacup.
    y_teacup : float
        The y-coordinate of the teacup.
    nmax : int
        The maximum number of points to calculate.
    bounce : bool
        If the projectile bounces on the ground.

    Returns
    -------
    horizontal_speed : float
        The horizontal speed of the projectile.
    '''
    if bounce:
        # set the value of gravity
        g = 980
        # set the range of the root
        high = x_teacup * np.sqrt(g / (2 * (initial_height - y_teacup )))
        low = high /(1 + 2 * np.sqrt((initial_height + y_teacup)/ (initial_height )))
        # Initial x-intercept
        mid = (low + high) / 2
        x, y = trajectory(mid, initial_height, x_teacup, nmax)
        # count the times of changing direction 
        times = times_changing_direction(x, y)
        # Loop until the root is found using midpoint method    
        while np.abs(residual(mid, initial_height, x_teacup, y_teacup,nmax)) > 0.01 or times != 2:
            if times == 2:
                if residual(low, initial_height, x_teacup, y_teacup,nmax)*\
                    residual(mid, initial_height, x_teacup, y_teacup,nmax) < 0:
                    high = mid
                else:
                    low = mid
            # more than one bounce showing a higher speed        
            elif times > 2:
                low = mid
            # less than one bounce showing a lower speed     
            else:
                high = mid
            mid = (low + high) / 2
            x, y = trajectory(mid, initial_height, x_teacup, nmax)
            # get the new times of changing direction
            times = times_changing_direction(x, y)    
        return mid
    else:
        # set the value of gravity
        g = 980
        horizontal_speed = x_teacup * np.sqrt(g / (2 * (initial_height - y_teacup )))
        # set the range of the root
        low = horizontal_speed * 0.8
        high = horizontal_speed * 1.2
        # Initial x-intercept by using falsi method
        mid = ((low*residual(high, initial_height, x_teacup, y_teacup,nmax)) - \
                high*residual(low, initial_height, x_teacup, y_teacup,nmax))/ \
                (residual(high, initial_height, x_teacup, y_teacup,nmax) - \
                residual(low, initial_height, x_teacup, y_teacup,nmax))
        # Loop until the root is found using falsi method
        while np.abs(residual(mid, initial_height, x_teacup, y_teacup,nmax)) > 0.01:
            # Check if the root is in the left or right half of the interval
            if residual(low, initial_height, x_teacup, y_teacup,nmax) * \
                residual(mid, initial_height, x_teacup, y_teacup,nmax) < 0:
                high = mid
            else:
                low = mid
            # Calculate the new x-intercept by using falsi method    
            mid = ((low*residual(high, initial_height, x_teacup, y_teacup,nmax)) - \
                   high*residual(low, initial_height, x_teacup, y_teacup,nmax))/ \
                    (residual(high, initial_height, x_teacup, y_teacup,nmax) - \
                     residual(low, initial_height, x_teacup, y_teacup,nmax))
        return mid

