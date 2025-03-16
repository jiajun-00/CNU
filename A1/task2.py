# Define all functions for Task 2 here.
import numpy as np

def price(widget, r, y, b):
    '''
    Function calculate the price of the `widget`
    :param weidget: `widget`, a list (possibly nested) representing a given widget
    :param r y b: epresenting respectively the prices of the `'RED'`, `'YELLOW'`, and `'BLUE'` raw materials, given in £ (GBP) 
    :return: the price of the `widget` rounded to the nearest £0.01
    '''
    Q = widget.copy()
    # Initialize the total price to zero
    sumprice = 0  
    while Q != []:
        # Check if the first widget is red
        if Q[0] == 'RED':
            # Add the price of a red widget
            sumprice = sumprice + r  
        # Check if the first widget is yellow
        elif Q[0] == 'YELLOW':  
            # Add the price of a yellow widget
            sumprice = sumprice + y  
        # Check if the first widget is blue
        elif Q[0] == 'BLUE':  
            # Add the price of a blue widget
            sumprice = sumprice + b  
        # if the first weight is not these primary weights, expand it
        else: 
            for i in Q[0]: 
                Q.append(i) 
        Q.pop(0) 
    return round(sumprice, 2)  # Return the total price rounded to two decimal places

def constituents(widget):
    '''
    This function count the number of each primary widgets
    param weidget: widget need to be split to primary widget
    return: the number of each primary widget red, yellow and blue
    '''
    Q = widget.copy()  # Make a copy of the widget list
    r = 0  
    y = 0 
    b = 0 
    while Q != []:
        # Check if the first widget is priamry and count it to the list
        if Q[0] == 'RED':  
            r = r + 1  
        elif Q[0] == 'YELLOW':  
            y = y + 1 
        elif Q[0] == 'BLUE':  
            b = b + 1 
        # if the first weight is not these primary weights, expand it    
        else:  
            for i in Q[0]:  
                Q.append(i)  
        Q.pop(0) 
    return [r,y,b]  # Return the counts of red, yellow, and blue widgets

def reconstruct_prices(orders, totals):
    '''
    The fuction will infer the price of each primary widget
    param orders: The given order from customers
    param: The price paid by them
    return: the price of each primary widget
    '''
    constituents_list = []  # Initialize an empty list to record primary widgets
    n = len(orders)  # Get the number of orders
    price = []  # Initialize an empty list to hold prices
    solution = [0,0,0]  # Initialize the solution to [0, 0, 0]
    # Record the primary widget for each order and make it as a matrix
    for i in orders:
        constituents_list.append(constituents(i))  
    # Transpose the constituents matrix
    A = np.stack(constituents_list)
    # Use least squares to solve for the prices
    price = np.linalg.lstsq(A,totals,rcond=None)  
    for i in range(0,len(price[0])):
        solution[i] = np.round(price[0][i],2)  # Round the price to two decimal places
    return solution  # Return the solution
