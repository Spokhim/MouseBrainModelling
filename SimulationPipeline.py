# Simulation Pipeline as Python Script
# Intended for use in SuperComputerCluster Automation
# Outputs Time Series Matrix, Parameters, FCM and SCM vs FCM Spearson Correlation 

# First Import the packages.
import numpy as np
import pylab
import matplotlib.pyplot as plt
# matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
from scipy import stats
# matplotlib.style.use('ggplot')
import inspect
import os
import csv 
import time

from tvb.simulator.lab import *
LOG = get_logger('demo')
# from tvb.simulator.plot.tools import *

# ## Main adjustments below, can hopefully get away with only changing cell below!

# Set up input Parameters for entire pipeline as dictionary elements here:

# Empty dict
ParamsDict = { }

# Name of import file/zip - Which contains connectivity data.
ParamsDict["name"] = "MouseCortex"

# Monitors or not?  (Aka BOLD or not?)
ParamsDict["BOLD"] = True

# Snip is the the number of elements you wish to snip at the start to remove the initial condition effects.
ParamsDict["Snip"] = 10
# Note, if BOLD = False, Snip gets multiplied by 1000. 
if ParamsDict["BOLD"] == False:
    ParamsDict["Snip"] = 1000 * ParamsDict["Snip"]

# Set the Random State/Seed for the Stochastic Integrator:
ParamsDict["RandState"] = 118

# Set Simulation Length:
ParamsDict["Simul_length"] = 1.2e6

# Set Linear Coupling Constant:
ParamsDict["G"] = np.array([0.096])

# Set integrator time step dt.
ParamsDict["dt"] = 0.1

# Set Additive Noise strength
ParamsDict["noise"] = np.array([0.000013])

# Set Wilson Cowan Model Parameters
ParamsDict["MODEL_c_ee"] = np.array([12.0])
ParamsDict["MODEL_c_ei"] = np.array([15.0])
ParamsDict["MODEL_c_ie"] = np.array([10.0])
ParamsDict["MODEL_c_ii"] = np.array([8.0])

# Define the model. 
ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"]) 

# Load the connectivity data from a zip file. 
con = connectivity.Connectivity.from_file(os.getcwd() +"/Connectomes/" + ParamsDict["name"] + ".zip")

# Mouse Allen Connectivity (from within TVB)
# con = connectivity.Connectivity.from_file("../mouse/allen_2mm/Connectivity.h5")

# Set the parameter of the resting state simulation
if ParamsDict["BOLD"] == True:
    sim = simulator.Simulator(model=ParamsDict["MODEL"],
                            connectivity=con,
                            coupling=coupling.Linear(a=ParamsDict["G"]),
                            integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                        random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                            monitors=(monitors.Bold(period=2e3),
                                      monitors.TemporalAverage(period=1e3)),
                            simulation_length=ParamsDict["Simul_length"],
                            #initial_conditions=[1.8,1.8,1.8,1.8,1.8]
                            ).configure()
    # Run the resting state simulation
    (bold_time, bold_data), _ = sim.run()

# No Monitors 
else:

  sim = simulator.Simulator(model=ParamsDict["MODEL"],
                          connectivity=con,
                          coupling=coupling.Linear(a=ParamsDict["G"]),
                          integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                        random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                          simulation_length=ParamsDict["Simul_length"]).configure()
  # Run the resting state simulation
  awer = sim.run()
  bold_time = awer[0][0]
  bold_data = awer[0][1]

# Got lazy
Snip = ParamsDict["Snip"]

# Functional Connectivity Matrix. 
# We note that this is a static analysis.  More advanced version would be a DYnamic version with windowing.

# Convert Simulation output into a form usable by Numpy.
TSeriesMatrix = np.empty((bold_data.shape[2], bold_data.shape[0]-Snip))
#Matrix = []

for i in range(bold_data.shape[2]):
    TSeriesMatrix[i] = bold_data[Snip:,0,i].flatten()

# Functional Conenctivity MAtrix = Pearson Correlation.

FCM = np.corrcoef(TSeriesMatrix)

# Comparing SC Vs FC with Spearman Correlation

# Grab Upper triangles
FCM_Upper = FCM[np.triu_indices(FCM.shape[0], k = 1)]
SCM_Upper = con.weights[np.triu_indices(con.weights.shape[0], k = 1)]

# Spearman Correlation
SCorr = stats.spearmanr(a=FCM_Upper,b=SCM_Upper)
print(SCorr)

# Export the simulation

time_now = time.strftime("%Y%m%d-%H%M%S")

# Params Dictionary - Note how we sort the dictionary.
with open("do-not-track/" + ParamsDict["name"] + "_Params_" + time_now + "_.csv", "w") as outfile:
    writer = csv.writer(outfile)
    for key, val in sorted(ParamsDict.items()):
        writer.writerow([key, val])
    
# Create Time Series and save. 
TSeries = np.concatenate((bold_time[Snip:].reshape(1,len(bold_time[Snip:])),TSeriesMatrix))
np.savetxt("do-not-track/" + ParamsDict["name"] + "_Tseries_" + time_now + "_.csv", TSeries, delimiter="\t")
np.savetxt("do-not-track/" + ParamsDict["name"]  + "_FCM_" + time_now + "_.csv", FCM, delimiter = "\t")
np.savetxt("do-not-track/" + ParamsDict["name"]  + "_Scorr_" +  time_now + "_.csv", SCorr, delimiter = "\t")  