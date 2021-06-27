# This Python File contains all the functions to generate the appropriate sub-figures in the paper.
# Run Reproducibility.ipynb to see all the figures and how they're generated (before they are edited with inkscape for use in publication)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pylab
import matplotlib.pyplot as plt
from scipy import stats
import inspect
import os
import csv 
import time

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *
from matplotlib.colors import ListedColormap
from turbo_colormap import *

from scipy.io import loadmat

###  Figure 1: 

# Fig 1. b) - SCM 

def Show_SCM():

    # Empty dict
    ParamsDict = { }
    ParamsDict["name"] = "MouseCortex"
    #ParamsDict["G"] = np.array([G_value]) 
    ParamsDict["REMOVE"] = [7]
    ParamsDict["BINARY"]=True

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
    SCM = con.weights

    # ListedColormap(turbo_colormap_data)
    #matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
    SCM = 1/SCM
    # Rearrange SCM
    SCM = SCM [index_fg] [:,index_fg]

    cs=plt.imshow(SCM, cmap="jet", aspect='equal', interpolation='none')
    #plt.title('Structural connectivity matrix', fontsize=20)
    #axcb=plt.colorbar(cs)
    #axcb.set_label('Weight', fontsize=20)
    #plt.yticks([0,1,2,3])
    #axcb.ax.tick_params(labelsize=16)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.xlabel("Region", fontsize=20)
    plt.ylabel("Region", fontsize=20)
    #plt.savefig("do-not-track\\SCM.pdf",bbox_inches='tight')
    plt.show()

    return


def Show_FCM():
    # Fig 1. g)- FCM 

    annots = loadmat('do-not-track\\timeSeriesData.mat') 
    DataMatrix = annots['timeSeriesData']

    # What needs to happen is to obtain each matrix (100 total for the no. of mice)

    # Create empty 3D Matrix to hold FC analysis
    FCMatrix = np.empty((38,38,100))

    for i in range(100):
        # For now we only take the cortical areas, so the first 38 rows.
        run = DataMatrix[0:38,:,i]

        # np.size(DataMatrix[:,:,0]) # Verify correct. 

        # Run Analysis - Pearson Correlation.
        FCM = np.corrcoef(run)
        FCMatrix[:,:,i] = FCM

    # Average over the 100 mice. 
    FCAverage = FCMatrix.mean(2)

    # Plot Matrix
    # For Individual slice of FC Matrix: FCMatrix[:,:,0]
    np.fill_diagonal(FCAverage,np.nan)

    # Remove Frontal Pole
    FCAverage = np.delete(FCAverage,obj=7,axis=0)
    FCAverage = np.delete(FCAverage,obj=7,axis=1)

    # Re-arrange the order to a new order, Ben's Functional Grouping
    index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])

    FCAverage = FCAverage[index_fg] [:,index_fg]

    cs=plt.imshow(FCAverage, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none')
    #plt.title('Functional connectivity matrix', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('Correlation', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.xlabel("Region", fontsize=20)
    plt.ylabel("Region", fontsize=20)
    #plt.savefig("do-not-track\\FCM_exp.pdf",bbox_inches='tight')
    plt.show()

    return
