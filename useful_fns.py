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
import glob
import pandas as pd

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *

# Input Simulation Pipeline
from SimulationPipeline import *

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

def params_extractor(file):
    """[summary]

    Parameters
    ----------
    file : csv file
        File contianing the parameters and their values.

    Returns
    -------
    dict
        ParamsDict contains the parameters and their values in the form of a dictionary.
    """


    return ParamsDict

def current_calculator(V,G,SCM):
    """
    Calculates the external current entering each node.

    Parameters
    ----------
    V : Time Series Matrix
        
    G : Coupling constant

    SCM : Structural Connectivity Matrix

    Returns
    -------
    tuple
        J_med,J_min,J_max
    """
    
    V_med = list(map(np.median, V))
    V_min = list(map(np.min, V))
    V_max = list(map(np.max, V))

    J = V_med * G * SCM 
    J_med = list(map(sum, V_med * G * SCM))
    J_min = list(map(sum, V_min * G * SCM))
    J_max = list(map(sum, V_max * G * SCM))

    print(J_med)

    plt.plot(J_med)
    plt.plot(J_min)
    plt.plot(J_max)
    plt.show()

    '''
    J = []

    for i in np.arange(V.shape[1]):
        item = V[:500,i]
        J_i = list(map(sum,item * G * SCM)) 
        J.append(J_i)
        # print(J)
        plt.plot(J_i)
    plt.show()
    print(len(J[0]))

    plt.plot(J)
    plt.show()

    J_med = 1
    J_min = 2
    J_max = 3
    '''
    # This returns as a tuple!  Oooh! (not a list)
    return J_med,J_min,J_max

