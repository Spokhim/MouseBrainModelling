# The goal of this Python script is to loop through all of the simulation files and obtain the J_e ranges.
# It then  figures out if it goes past the bifurcation point or not. 

# Import Packages
import numpy as np
import pandas as pd
import pylab
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.colors import ListedColormap
from turbo_colormap import *
import inspect
import os
import csv 
import time
import sys
import glob
import pandas as pd
from pprint import pprint
import re
# RegEx module

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *

# Input Simulation Pipeline
from SimulationPipeline import *
from useful_fns import *

#matplotlib.rcParams['figure.figsize'] = (10.0, 5.0)

# Get all the file names. 
TseriesFile = glob.glob(r"do-not-track\\2020_12_26\\LCycle*Tseries*.csv") 

# Import SCM, don't worry about re-ordering. Needed to Calculate J_e

# Empty dict
ParamsDict = { }
# Name of import file/zip - Which contains connectivity data.
ParamsDict["name"] = "MouseCortex"
ParamsDict["REMOVE"] = [7]
ParamsDict["BINARY"] = True

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

Params = []
TheMin = []
TheMax = []
TheMedian = [] 

# Now for Each TseriesFile
for string in TseriesFile:

    # Obtain Parameter Values
    x = re.findall("\[(.*)\].*\[(.*)\]",string)
    #x = re.findall("\[(.*)\].b_e(...)",string)
    Params.append(x[0])

    # Want to get the J_e ranges

    # Import the Tseries from the file. 
    df = np.genfromtxt(TseriesFile,delimiter="\t")

    bold_time = df[0]
    bold_data = df[1:]


    # Calculate J_e
    # External Current Calculator:
    J_e = []
    # j is jth element

    # Make it faster by only checking the first 1/50 as this is 1.2e6 time. 
    for j in np.arange(len(bold_time)/50):       
        t_0 = []
        # Specific column (or time point)
        for i in np.arange(SCM.shape[0]): 
            # Sum over all external currents (May need to do SCM[:,i] instead, but quite sure SCM[i,:] is correct from tseries analysis)  
            t  = sum(bold_data[:,j]*SCM[i,:])*G_value
            # To obtain currents to particular region
            t_0.append(t)
        J_e.append(t_0)

    J_e = np.array(J_e)

    # Get the Min, Max, Median
    TheMin.append(J_e.min())
    TheMax.append(J_e.max())
    TheMedian.append(J_e.median())

df = pd.DataFrame(Params)
df.columns = ['G', 'B_e']
df["TheMin"] = TheMin
df["TheMax"] = TheMax
df["TheMedian"] = TheMedian

# Drop any duplicates if necessary. 
df = df.drop_duplicates(['G','B_e'],keep='first')

# Export the df
df.to_csv('do-not-track\J_e_LCycle.csv',index=True)


#df_pivot_min = df.pivot(index='B_e', columns='G', values='TheMin')
#df_pivot_max = df.pivot(index='B_e', columns='G', values='TheMax')
#df_pivot_median = df.pivot(index='B_e', columns='G', values='TheMedian')