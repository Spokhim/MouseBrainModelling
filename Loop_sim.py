# We separated parts of the Original SimulationPipeline.py
# It's cause we don't really want to change SimulationPiepline.py too much.
# In this file, we define the Parameter Dictionary
# As well as set up the looping.
# For now we do the looping thorugh pipeline instead of through the PBS array as I don't get it.

# First Import the packages.
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
from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *

# Input Simulation Pipeline
from SimulationPipeline import *

################################################################################################################################
# These are Default Values for ParamsDict.  

# Empty dict
ParamsDict = { }

# Name of import file/zip - Which contains connectivity data.
ParamsDict["name"] = "MouseCortex"

# Calculate FC-FC for Mouse?
ParamsDict["FCFC"] = True

# Monitors or not?  (Aka BOLD or not?)
ParamsDict["BOLD"] = True

# Change to Binary Connectome? (If True, will change the connectome into binary)
ParamsDict["BINARY"] = True

# Snip is the the number of elements you wish to snip at the start to remove the initial condition effects.
ParamsDict["Snip"] = 10
# Note, if BOLD = False, Snip gets multiplied by 1000, later in the SimulationPipeline code.  Not the actual dictionary element though.  

# Set the Random State/Seed for the Stochastic Integrator:
ParamsDict["RandState"] = 118

# Set Simulation Length:
ParamsDict["Simul_length"] = 1.2e6

# Set Linear Coupling Constant:
ParamsDict["G"] = np.array([0.47]) 

# Set integrator time step dt.
ParamsDict["dt"] = 0.1

# Set Additive Noise strength
ParamsDict["noise"] = np.array([0.000013])

# Set Wilson Cowan Model Parameters
ParamsDict["MODEL_c_ee"] = np.array([12.0])
ParamsDict["MODEL_c_ei"] = np.array([15.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([8.0])

# Model is now defined within SimulationPipeline.py
# However if you adjusting parameters other than these Coupling Parameters, then you need to redefine the model in this file per run.

# Params Dict tag (extra note tags for the name - Example to denote what's being changed/looped.)
ParamsDict["tag"] = ""

################################################################################################################################

# i is PBS_ARRAY_INDEX - Allows for creation of multiple jobs 
i = int(sys.argv[1])

# Try Heterogeneous 
df = pd.read_csv("CortexDensities.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values
E_prop = E_pop / (E_pop + I_pop)

if i == 1: 
    ParamsDict["k_e"] = E_prop
    ParamsDict["tag"] = "k_e"
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),k_e=ParamsDict["k_e"],k_i=1-ParamsDict["k_e"])  
    Simul_Pipeline(ParamsDict=ParamsDict)
    print(ParamsDict["tag"] ,"Completed")

if i == 2:
    ParamsDict["k_e"] = 2*E_prop
    ParamsDict["tag"] = "2k_e"
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),k_e=ParamsDict["k_e"],k_i=1-ParamsDict["k_e"])  
    Simul_Pipeline(ParamsDict=ParamsDict)
    print(ParamsDict["tag"] ,"Completed")

if i == 3:
    E_prop[7] = 1 - E_prop[7]
    ParamsDict["k_e"] = E_prop
    ParamsDict["tag"] = "k_e_alter"
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),k_e=ParamsDict["k_e"],k_i=1-ParamsDict["k_e"])  
    Simul_Pipeline(ParamsDict=ParamsDict)
    print(ParamsDict["tag"] ,"Completed")

if i == 4:
    E_prop[7] = 1 - E_prop[7]
    ParamsDict["k_e"] = 2*E_prop
    ParamsDict["tag"] = "2k_e_alter"
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),k_e=ParamsDict["k_e"],k_i=1-ParamsDict["k_e"])  
    Simul_Pipeline(ParamsDict=ParamsDict)
    print(ParamsDict["tag"] ,"Completed")

if i == 5:
    E_prop[7] = 1 - E_prop[7]
    ParamsDict["k_e"] = 2*E_prop
    ParamsDict["BOLD"] = False
    ParamsDict["tag"] = "2k_e_alter_nonBOLD"
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),k_e=ParamsDict["k_e"],k_i=1-ParamsDict["k_e"])  
    Simul_Pipeline(ParamsDict=ParamsDict)
    print(ParamsDict["tag"] ,"Completed")

print("Happilly Finished All!")