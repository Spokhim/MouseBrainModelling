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
#from useful_fns import *

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

# Set Initial Conditions.  Comment out if you want TVB to make initial conditions for you.  Init_Cons CANT be in the Dict in this case.  
#ParamsDict["Init_Cons"] = initial_conditions=0.25 + numpy.zeros((37*37,2,37,1))
# Must be of form:
# initial_conditions=0.5 + numpy.zeros((con.number_of_regions*con.number_of_regions,2,con.number_of_regions,1))

# Set Wilson Cowan Model Parameters - Limit Cycle Cut
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])

# Model is now defined within SimulationPipeline.py
# However if you adjusting parameters other than these Coupling Parameters, then you need to redefine the model in this file per run.

ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([1.5]),b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 

# Params Dict tag (extra note tags for the name - Example to denote what's being changed/looped.)
ParamsDict["tag"] = ""

################################################################################################################################

# i is PBS_ARRAY_INDEX - Allows for creation of multiple jobs 
I = int(sys.argv[1])

ParamsDict["Simul_length"] = 1.2e5
#ParamsDict["Snip"] = 100 # Gets multiplied by 100 so this is 1 e4 (or 1%)

# Obtain Het Data
df = pd.read_csv("CortexDensitiesAlter.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values
E_mean = np.mean(E_pop)
I_mean = np.mean(I_pop)

# E_normalised is (when excluding region 7) -0.28 to 0.54
E_normalised = (E_pop-E_mean)/E_mean
# I_normalised is (when excluding region 7) -0.45 to 1.44
I_normalised = (I_pop-I_mean)/I_mean


# Gradient - Version 2 
E_prop = E_pop / (E_pop + I_pop)
I_prop = I_pop / (E_pop + I_pop)

Mean_e_prop = np.mean(E_prop)
Mean_i_prop = np.mean(I_prop)

# -0.1 to 0.08 roughly. 
E_prop_norm = (E_prop-Mean_e_prop)/Mean_e_prop
# -0.4 to 0.7
I_prop_norm = (I_prop-Mean_i_prop)/Mean_i_prop

E_normalised = E_prop_norm
I_normalised = I_prop_norm 

# Simulate LCycle Het 1.2 e5 no G sweep and finer sigma. 
# LCycle
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])
ParamsDict["B_e"] = np.array([1.5]) 
ParamsDict["G"] = np.array([0.75]) 
# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

# Simualte V2 het (so using proportions.)

ParamsDict["Sigma_e"] =I*0.1*5
sigma_e = ParamsDict["Sigma_e"] 

ParamsDict["Sigma_i"] = I*0.1
sigma_i = ParamsDict["Sigma_i"] 

        # Heterogeneous Coupling Constants (array)
ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma_i * I_normalised) 
ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma_i * I_normalised) 
        
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycle_HetV2_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) + "_sig_e" + str(sigma_e) + "_sig_i" + str(sigma_i)
Simul_Pipeline(ParamsDict=ParamsDict)

# Simualte V3 het (so using proportions and adjusting a_e, a_i)
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])

ParamsDict["Sigma_e"] =I*0.1*5
sigma_e = ParamsDict["Sigma_e"] 

ParamsDict["Sigma_i"] = I*0.1
sigma_i = ParamsDict["Sigma_i"] 

        # Heterogeneous Coupling Constants (array)
ParamsDict["MODEL_a_e"] = numpy.array([1.0]) * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_a_i"] = numpy.array([1.0])  * (1 + sigma_i * I_normalised) 
        
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=ParamsDict["MODEL_a_e"],a_i=ParamsDict["MODEL_a_i"],b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycle_HetV3_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) + "_sig_e" + str(sigma_e) + "_sig_i" + str(sigma_i)
Simul_Pipeline(ParamsDict=ParamsDict)

# Finally try Het when I remove Nodes 26,19,3 (of the original order)

ParamsDict["REMOVE"] = [26,19,7,3]

df = pd.read_csv("CortexDensitiesAlter.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values

E_pop = np.delete(E_pop,[26,19,3])
I_pop = np.delete(I_pop,[26,19,3])
E_mean = np.mean(E_pop)
I_mean = np.mean(I_pop)

# E_normalised is (when excluding region 7,26,19,3) -0.27 -> 0.41
E_normalised = (E_pop-E_mean)/E_mean
# I_normalised is (when excluding region 7,26,19,3) -0.4 -> 0.63
I_normalised = (I_pop-I_mean)/I_mean

# With Alt LCycle B_e = 1.5
# Simulate LCycle Het 1.2 e5 no G sweep and finer sigma. 
# LCycle
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])
ParamsDict["B_e"] = np.array([1.5]) 
ParamsDict["G"] = np.array([0.75]) 
# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

ParamsDict["Sigma_e"] =I*0.1
sigma_e = ParamsDict["Sigma_e"] 

ParamsDict["Sigma_i"] = I*0.1
sigma_i = ParamsDict["Sigma_i"] 

        # Heterogeneous Coupling Constants (array)
ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma_i * I_normalised) 
ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma_i * I_normalised) 
        
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycleAlt_Het_Removed_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) + "_sig_e" + str(sigma_e) + "_sig_i" + str(sigma_i)
Simul_Pipeline(ParamsDict=ParamsDict)

# With Best LCycle B_e = 0.8
# Simulate LCycle Het 1.2 e5 no G sweep and finer sigma. 
# LCycle
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])
ParamsDict["B_e"] = np.array([0.8]) 
ParamsDict["G"] = np.array([0.5]) 
# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

ParamsDict["Sigma_e"] =I*0.1
sigma_e = ParamsDict["Sigma_e"] 

ParamsDict["Sigma_i"] = I*0.1
sigma_i = ParamsDict["Sigma_i"] 

        # Heterogeneous Coupling Constants (array)
ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma_e * E_normalised) 
ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma_i * I_normalised) 
ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma_i * I_normalised) 
        
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycleB_Het_Removed_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) + "_sig_e" + str(sigma_e) + "_sig_i" + str(sigma_i)
Simul_Pipeline(ParamsDict=ParamsDict)

'''
# B_e: 0 -> 3
for b_e in np.arange(start=0,stop=3,step=0.1): 
    # LCycle
    ParamsDict["MODEL_c_ee"] = np.array([11.0])
    ParamsDict["MODEL_c_ei"] = np.array([10.0])
    ParamsDict["MODEL_c_ie"] = np.array([10.0])
    ParamsDict["MODEL_c_ii"] = np.array([1.0])
    ParamsDict["B_e"] = np.array([b_e]) 
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([65.0])) 
    ParamsDict["tag"] = "LCycle_BOLD_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)

# B_e goes from 2->5
for b_e in np.arange(start=2,stop=5,step=0.1): 

    # Fixed Pt
    ParamsDict["MODEL_c_ee"] = np.array([12.0])
    ParamsDict["MODEL_c_ei"] = np.array([15.0])
    ParamsDict["MODEL_c_ie"] = np.array([10.0])
    ParamsDict["MODEL_c_ii"] = np.array([8.0])
    ParamsDict["B_e"] = np.array([b_e]) 
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([4]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "FixedPt_BOLD_G" + str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)

#B_e : 3- >5
for b_e in np.arange(start=3,stop=5,step=0.1): 

    # Hysteresis
    ParamsDict["MODEL_c_ee"] = np.array([16.0])
    ParamsDict["MODEL_c_ei"] = np.array([12.0])
    ParamsDict["MODEL_c_ie"] = np.array([10.0])
    ParamsDict["MODEL_c_ii"] = np.array([3.0])
    ParamsDict["B_e"] = np.array([b_e]) 
    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                            a_e=numpy.array([1.3]),a_i=numpy.array([2.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([3.7]),tau_e=numpy.array([10.0]),
                                            tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "Hysteresis_BOLD_G"+ str(ParamsDict["G"]) + "_b_e" + str(ParamsDict["B_e"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)
'''

'''
# LCycle-RHopf = LCycleCut
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])
ParamsDict["B_e"] = np.array([0.8])
ParamsDict["G"] = np.array([0.5]) 
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=ParamsDict["B_e"],b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycleCut_i" + str(i) + "Length_e5"
Simul_Pipeline(ParamsDict=ParamsDict)
'''
'''
# Gradient - Version 2 
E_prop = E_pop / (E_pop + I_pop)
I_prop = I_pop / (E_pop + I_pop)

Mean_e_prop = np.mean(E_prop)
Mean_i_prop = np.mean(I_prop)

# -0.1 to 0.08 roughly. 
E_prop_norm = (E_prop-Mean_e_prop)/Mean_e_prop
# -0.4 to 0.7
I_prop_norm = (I_prop-Mean_i_prop)/Mean_i_prop

E_normalised = E_prop_norm
I_normalised = I_prop_norm 
'''

"""
# Obtain Het Data
df = pd.read_csv("CortexDensitiesAlter.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values
E_mean = np.mean(E_pop)
I_mean = np.mean(I_pop)

# E_normalised is (when excluding region 7) -0.28 to 0.54
E_normalised = (E_pop-E_mean)/E_mean
# I_normalised is (when excluding region 7) -0.45 to 1.44
I_normalised = (I_pop-I_mean)/I_mean

# Set Wilson Cowan Model Parameters - LCycleH - LCycle Parameters close to Hysteresis Case. 
ParamsDict["MODEL_c_ee"] = np.array([16.0])
ParamsDict["MODEL_c_ei"] = np.array([12.0])
ParamsDict["MODEL_c_ie"] = np.array([15.0])
ParamsDict["MODEL_c_ii"] = np.array([3.0])

# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

for J in np.arange(6):
    ParamsDict["sigma"] =J*0.2
    sigma = ParamsDict["sigma"] 

    # Heterogeneous Coupling Constants (array)
    ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma * I_normalised) 
    ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma * I_normalised) 

    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.3]),a_i=numpy.array([2.0]),b_e=numpy.array([4]),b_i=numpy.array([3.7]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "LCycleH_G" + str(ParamsDict["G"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)

# Set Wilson Cowan Model Parameters - LCycleHCut - LCycleCut Parameters close to Hysteresis Case, but chopped so hopf bifurcation.
ParamsDict["MODEL_c_ee"] = np.array([16.0])
ParamsDict["MODEL_c_ei"] = np.array([12.0])
ParamsDict["MODEL_c_ie"] = np.array([15.0])
ParamsDict["MODEL_c_ii"] = np.array([3.0])

# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

for J in np.arange(6):
    ParamsDict["sigma"] =J*0.2
    sigma = ParamsDict["sigma"] 

    # Heterogeneous Coupling Constants (array)
    ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma * I_normalised) 
    ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma * I_normalised) 

    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.3]),a_i=numpy.array([2.0]),b_e=numpy.array([3]),b_i=numpy.array([3.7]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "LCycleH_G" + str(ParamsDict["G"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)
"""

"""
# Set Wilson Cowan Model Parameters - Hysteresis
ParamsDict["MODEL_c_ee"] = np.array([16.0])
ParamsDict["MODEL_c_ei"] = np.array([12.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([3.0])

# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

# Hysteresis
for J in np.arange(6):
    ParamsDict["sigma"] =J*0.2
    sigma = ParamsDict["sigma"] 

    # Heterogeneous Coupling Constants (array)
    ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma * I_normalised) 
    ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma * I_normalised) 

    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.3]),a_i=numpy.array([2.0]),b_e=numpy.array([4]),b_i=numpy.array([3.7]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "Hysteresis_G" + str(ParamsDict["G"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)

# FixedPt-eqbm

# Set Wilson Cowan Model Parameters - Fixed pt
ParamsDict["MODEL_c_ee"] = np.array([12.0])
ParamsDict["MODEL_c_ei"] = np.array([15.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([8.0])

# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

for J in np.arange(6):
    ParamsDict["sigma"] =J*0.2
    sigma = ParamsDict["sigma"] 

    # Heterogeneous Coupling Constants (array)
    ParamsDict["MODEL_c_ie"] = h_ie * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ee"] = h_ee  * (1 + sigma * E_normalised) 
    ParamsDict["MODEL_c_ii"] = h_ii  * (1 + sigma * I_normalised) 
    ParamsDict["MODEL_c_ei"] = h_ei  * (1 + sigma * I_normalised) 

    ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                        a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([4]),b_i=numpy.array([4]),tau_e=numpy.array([10.0]),
                                        tau_i=numpy.array([10.0])) 
    ParamsDict["tag"] = "FixedPt_G" + str(ParamsDict["G"]) 
    Simul_Pipeline(ParamsDict=ParamsDict)
"""
"""
# Limit Cycle Params - LCycleCut
ParamsDict["MODEL_c_ee"] = np.array([11.0])
ParamsDict["MODEL_c_ei"] = np.array([10.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([1.0])

# Homogeneous Coupling constants
h_ee = ParamsDict["MODEL_c_ee"] 
h_ei = ParamsDict["MODEL_c_ei"] 
h_ie = ParamsDict["MODEL_c_ie"] 
h_ii = ParamsDict["MODEL_c_ii"] 

ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                    a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([1.5]),b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                    tau_i=numpy.array([65.0])) 
ParamsDict["tag"] = "LCycleCut_G" + str(ParamsDict["G"]) 
Simul_Pipeline(ParamsDict=ParamsDict)
"""