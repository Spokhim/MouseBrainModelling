# This is code for Shuffling the Heterogeneity gradient randomly
# And then simulating.  This is to validate that the heterogenous gradient is meaningful.

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
from useful_fns import *

################################################################################################################################
# These are Default Values for ParamsDict.  

# Empty dict
ParamsDict = { }

# Name of import file/zip - Which contains connectivity data.
ParamsDict["name"] = "MouseCortex"

# Calculate FC-FC for Mouse?
ParamsDict["FCFC"] = True

# Export Simulation Files?
ParamsDict["ExportSim"] = True

# Monitors or not?  (Aka BOLD or not?)
ParamsDict["BOLD"] = False

# Change to Binary Connectome? (If True, will change the connectome into binary)
ParamsDict["BINARY"] = True

# Snip is the the number of elements you wish to snip at the start to remove the initial condition effects.
ParamsDict["Snip"] = 10
# Note, if BOLD = False, Snip gets multiplied by 100, later in the SimulationPipeline code.  Not the actual dictionary element though.  

# Set the Random State/Seed for the Stochastic Integrator:
ParamsDict["RandState"] = 118

# Remove ith indexed region (7 corresponds to Frontal Pole Cerebral Cortex) - Give it a list if removing multiple regions.  Empty list removes nothing. 
ParamsDict["REMOVE"] = [7]

# Set Simulation Length:
ParamsDict["Simul_length"] = 1.2e5

# Set Linear Coupling Constant:
ParamsDict["G"] = np.array([0.47]) 

# Set integrator time step dt.
ParamsDict["dt"] = 0.1

# Set Additive Noise strength
ParamsDict["noise"] = np.array([0.000013])

# Set Wilson Cowan Model Parameters - Hysteresis
ParamsDict["MODEL_c_ee"] = np.array([16.0])
ParamsDict["MODEL_c_ei"] = np.array([12.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([3.0])

# Model is now defined within SimulationPipeline.py
# However if you adjusting parameters other than these Coupling Parameters, then you need to redefine the model in this file per run.

ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([4]),b_i=numpy.array([4]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 

# Params Dict tag (extra note tags for the name - Example to denote what's being changed/looped.)
ParamsDict["tag"] = ""

################################################################################################################################

# i is PBS_ARRAY_INDEX - Allows for creation of multiple jobs 
i = int(sys.argv[1])

# FCFC Shuffle
ParamsDict["ExportSim"] = False 
# Skimp on computation time and power by reducing to 1.2e4.  Should be fine. 
ParamsDict["Simul_length"] = 1.2e4

# First we must shuffle it randomly. 
# Let's Pause by X seconds to make sure the random stuff is working and different across the parallel cpus (based on  sys time after all)
time.sleep(i*0.1)

df = pd.read_csv("CortexDensitiesAlter.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values
E_mean = np.mean(E_pop)
I_mean = np.mean(I_pop)

# E_normalised is (when excluding region 7) -0.28 to 0.54
E_norm = (E_pop-E_mean)/E_mean
# I_normalised is (when excluding region 7) -0.45 to 1.44
I_norm = (I_pop-I_mean)/I_mean

# Shuffled:
perm = numpy.random.permutation(len(E_norm))
E_normalised = E_norm[perm]
I_normalised = I_norm[perm]

# Limit Cycle Params
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])
# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

Best_Score = 0 
Best_G = 0
Best_Sigma = 0  

# Sweep across the range of Sigma values:
for I in np.arange(6):

    ParamsDict["Sigma"] =I*0.2
    sigma = ParamsDict["Sigma"] 

    # Heterogeneous Coupling Constants (array)
    ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma * I_normalised) 
    ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma * I_normalised) 

    
    # Sweep across the range of G values
    for J in np.arange(41):
        ParamsDict["G"] = np.array([J * 0.05])

        ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                            a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([1.5]),b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                            tau_i=numpy.array([65.0])) 
        Score = Simul_Pipeline(ParamsDict=ParamsDict)[2]

        # If the score is the best score, store it. 
        if Score > Best_Score:
            Best_Score = Score
            Best_G = ParamsDict["G"]
            Best_Sigma = ParamsDict["Sigma"] 
            print(Best_Score)  
            print(Best_G)   
            print(Best_Sigma)         

# Now export the information:
time_now = time.strftime("%Y%m%d-%H%M%S")
np.savetxt("do-not-track/" + str(i) + "_" + ParamsDict["name"] + "_Best_" + time_now + "_.csv", [Best_Score,Best_G,Best_Sigma], delimiter="\t")
np.savetxt("do-not-track/" + str(i) + "_" + ParamsDict["name"] + "_EIHet_" + time_now + "_.csv", [E_normalised,I_normalised], delimiter="\t")