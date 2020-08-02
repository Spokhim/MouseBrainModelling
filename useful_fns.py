    """[summary]
    A bunch of useful functions.    
    """

import numpy as np
import pandas as pd
import pylab
import matplotlib.pyplot as plt
from scipy import stats
import inspect
import os
import csv 
import time
import sys

def par_extract(file):
    # Extracts the characters between square brackets.
    s_filter = ""
    x = False
    for i in file:
        if i == "[": x = True
        elif i == "]": x = False
        elif x: s_filter += i

    return s_filter

def sorter(X,Y):
    # Sort X based on Y
    Z = [x for _,x in sorted(zip(Y,X))]
    return Z

def current_calculator(V,G,SCM):
    """

    Parameters
    ----------
    V : [type]
        [description]
    """
    V_med = np.median(V)
    V_min = np.min(V)
    V_max = np.max(V)

    J_med = np.median(V)
    J_min = np.min(V)
    J_max = np.max(V)


    # This returns as a tuple!  Oooh! (not a list)
    return J_med,J_min,J_max