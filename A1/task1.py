# Define all functions for Task 1 here.
# Any necessary import statements should go at the top of each module.
# For example, if you need Numpy in your Task 1 functions:
import numpy as np
import random

def check_floor_root(number, root, power = 2):
    """
    This function checks check whether r is the p-th floor root of n.
    :param number: The number need to be checked
    :param root: The floor root
    :param power: The power of the floor root
    :return: True if r is the p-th floor root of n.
    """
    if root ** (power) <= number and number < (root + 1) ** power:
        return True
    else:
        return False


def unsafe_floor_sqrt(x):
    return int((x + 0.5)**(1/2))

def random_number(digits):
    """
    This function returns a random number with the specified number of digits.
    :param digits: The number of digits for the random number.
    :return: A random number with the specified number of digits.
    """ 
    return random.randint(10**(digits-1),10**(digits)-1) 

def unsafe_failure_rate(number_sizes, samples = 500):
    """
    This function calculates the failure rate of the unsafe_floor_sqrt function for a list of digits.
    :param number_sizes: The list of digits to calculate the failure rate for each digit.
    :param samples: The number of samples to use for each number of digits.
    :return: A list of the failure rates for each number of digits.
    """
    frequencies = []
    for i in number_sizes:
        # count the failure numbers of each number
        timesFaliure = 0
        # Test this digit for samples times
        for j in range(0,samples):
            random = random_number(i)
            floor_sqrt = unsafe_floor_sqrt(random)
            # count the error times
            if not check_floor_root(random, floor_sqrt):
                timesFaliure = timesFaliure + 1
        # get the final frequency for this number of digits        
        frequencies.append(timesFaliure/samples)      
    return frequencies


def floor_square_root(number):
    """
    This function returns the floor of the square root of a number based on given algorithm.
    :param number: The number to find the square root of.
    :return: The floor of the square root of the number.
    """
    digits = 0 
    left_n = number
    # get the digits of number
    while left_n != 0:
        left_n = left_n//10
        digits = digits +1    
    largest_even = digits//2*2
    root = 0
    # get the root of number from the largest digit
    for i in range(largest_even, -1, -2):
        remain = number // (10**i)
        # enable the next digit to be calculated through known digits
        root = root * 10
        choice = 0
        # get the next digit by using the known digits
        for j in range(0,10):
            # find the floor square of remaining numbers with power 2            
            if remain >= (root + j)**2 and remain<(root + j + 1)**2:
                choice = j
                break   
        root = root + choice
    return root 

def floor_root(number, power=2):
    """
    This function returns the floor of the square root of a number based on given algorithm with given power.
    :param number: The number to find the square root of.
    :param power: The power to use for the nth root.
    :return: The floor of the square root of the number.
    """
    digits = 0
    left_n = number
    # get the digits of n
    while left_n != 0:
        left_n = left_n//10
        digits = digits +1   
    # get the times of repeating loop     
    largest_power = digits//power*power
    root = 0
    # get the digits by going from largest power to 0 by power steps
    for i in range(largest_power, -1, -power):
        remain = number // (10**i)
        root = root * 10
        choice = 0
        # test from 0 to 9 to find the suitable digit for current
        for j in range(0,10):
            # find the floor square of remaining numbers with power
            if remain >= (root + j)**power and remain<(root + j + 1)**power:
                choice = j
                break   
        root = root + choice
    return root 



    





