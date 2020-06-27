# Simulation Pipeline as Python Script
# Intended for use in SuperComputerCluster Automation
# Outputs Time Series Matrix, Parameters, FCM and SCM vs FCM Spearson Correlation 

# First Import the packages.
import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy import stats
import inspect
import os
import csv 
import time

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *

def Simul_Pipeline(ParamsDict):
    '''
    Simulation Pipeline, Does the simulations and returns Scorr of FC vs SC.  Saves some csvs in do-no-track folder.  Takes in inputs as ParamsDict
    '''

    # Load the connectivity data from a zip file. 
    con = connectivity.Connectivity.from_file(os.getcwd() +"/Connectomes/" + ParamsDict["name"] + ".zip")

    # Mouse Allen Connectivity (from within TVB)
    # con = connectivity.Connectivity.from_file("../mouse/allen_2mm/Connectivity.h5")

    # Change to Connectome to Binary if desired:
    if ParamsDict["BINARY"]==True:
        con.weights = con.weights!=0

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
    # Note, if BOLD = False, Snip gets multiplied by 1000. 
    if ParamsDict["BOLD"] == False:
        Snip = 1000 * Snip

    # Functional Connectivity Matrix. 
    # We note that this is a static analysis.  More advanced version would be a Dynamic version with windowing.

    # Convert Simulation output into a form usable by Numpy.
    TSeriesMatrix = np.empty((bold_data.shape[2], bold_data.shape[0]-Snip))

    for i in range(bold_data.shape[2]):
        TSeriesMatrix[i] = bold_data[Snip:,0,i].flatten()

    # Functional Conenctivity MAtrix = Pearson Correlation.

    FCM = np.corrcoef(TSeriesMatrix)

    # Set diagonals to NaN
    FCM1 = FCM
    np.fill_diagonal(FCM1,np.nan)

    # Comparing SC vs FC with Spearman Corr
    # Check if SCM is symmetric: 
    SCM = con.weights
    Sym_check = numpy.allclose(SCM, SCM.T,equal_nan=True)

    if Sym_check == True:
        #It is a symmetric SCM, so only use upper triangles
        # Grab Upper triangles
        FCM_Upper = FCM[np.triu_indices(FCM.shape[0], k = 1)]
        SCM_Upper = con.weights[np.triu_indices(con.weights.shape[0], k = 1)]

    elif Sym_check == False:
        # If SCM is not symmetric, need to calcualte spearman corr for entire matrix.
        # Set Diagonal to Nans
        np.fill_diagonal(SCM,np.nan)
        # Remove all Nans for SCM and FCM
        SCM_Upper = SCM[~numpy.isnan(SCM)]
        FCM_Upper = FCM1[~numpy.isnan(FCM1)]

    # Spearman Correlation
    SCorr = stats.spearmanr(a=FCM_Upper,b=SCM_Upper)
    print(SCorr)

    # Export the simulation

    time_now = time.strftime("%Y%m%d-%H%M%S")

    # Params Dictionary - Note how we sort the dictionary.
    with open("do-not-track/" + ParamsDict["tag"] + "_" + ParamsDict["name"] + "_Params_" + time_now + "_.csv", "w") as outfile:
        writer = csv.writer(outfile)
        for key, val in sorted(ParamsDict.items()):
            writer.writerow([key, val])
        
    # Create Time Series and save. 
    TSeries = np.concatenate((bold_time[Snip:].reshape(1,len(bold_time[Snip:])),TSeriesMatrix))
    np.savetxt("do-not-track/" + ParamsDict["tag"] + "_" + ParamsDict["name"] + "_Tseries_" + time_now + "_.csv", TSeries, delimiter="\t")
    np.savetxt("do-not-track/" + ParamsDict["tag"] + "_" + ParamsDict["name"]  + "_FCM_" + time_now + "_.csv", FCM, delimiter = "\t")
    np.savetxt("do-not-track/" + ParamsDict["tag"] + "_" + ParamsDict["name"]  + "_Scorr_" +  time_now + "_.csv", SCorr, delimiter = "\t")  

    return

