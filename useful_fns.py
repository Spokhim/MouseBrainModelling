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

# Need a params extractor function.

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

'''
# Empty dict
ParamsDict = { }
ParamsDict["name"] = "MouseCortex"
ParamsDict["G"] = np.array([1.9]) 
ParamsDict["REMOVE"] = [7]
ParamsDict["BINARY"]=True

# Ye dunno why having the closing square bracket messes up glob glob. 
Sim_run_files = glob.glob("do-not-track/R_LCycle_G[0.4*_.csv")
print(Sim_run_files)

# Read file import data
#df = pd.read_csv(all_files[11],delimiter="\t",header=None)
# Genfromtxt gives us a np array. 
df = np.genfromtxt(Sim_run_files[-1],delimiter="\t")

bold_time = df[0]
bold_data = df[1:]

# Load the connectivity data from a zip file. 
con = connectivity.Connectivity.from_file(os.getcwd() +"/Connectomes/" + ParamsDict["name"] + ".zip")

# Remove the ith row and column in centres, tract_lengths and weights. i.e. the specified region(s)
con.centres = np.delete(con.centres,ParamsDict["REMOVE"])
con.weights = np.delete(con.weights,obj=ParamsDict["REMOVE"],axis=0)
con.weights = np.delete(con.weights,obj=ParamsDict["REMOVE"],axis=1)
con.tract_lengths = np.delete(con.tract_lengths,obj=ParamsDict["REMOVE"],axis=0)
con.tract_lengths = np.delete(con.tract_lengths,obj=ParamsDict["REMOVE"],axis=1)

if ParamsDict["BINARY"]==True:
    con.weights = con.weights!=0

current_calculator(bold_data,ParamsDict["G"],con.weights)
'''