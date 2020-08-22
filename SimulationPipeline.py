# Simulation Pipeline as Python Script
# Intended for use in SuperComputerCluster Automation
# Outputs Time Series Matrix, Parameters, FCM and SCM vs FCM Spearson Correlation 

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

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *

def Simul_Pipeline(ParamsDict):
    '''
    Simulation Pipeline, Does the simulations and returns Scorr of FC vs SC.  Saves some csvs in do-no-track folder.  Takes in inputs as ParamsDict
    '''

    # Define the model if it is not yet defined.
    if "MODEL" not in ParamsDict:
        ParamsDict["MODEL"] = models.WilsonCowan(c_ee=ParamsDict["MODEL_c_ee"],c_ei=ParamsDict["MODEL_c_ei"],c_ie=ParamsDict["MODEL_c_ie"] ,c_ii=ParamsDict["MODEL_c_ii"],
                                                a_e=numpy.array([1.0]),a_i=numpy.array([1.0]),b_e=numpy.array([2.8]),b_i=numpy.array([2.8]),tau_e=numpy.array([10.0]),
                                                tau_i=numpy.array([65.0]),) 

    # Load the connectivity data from a zip file. 
    con = connectivity.Connectivity.from_file(os.getcwd() +"/Connectomes/" + ParamsDict["name"] + ".zip")

    # Now need to prepare the connectivity data accordingly.  Unfortuantely doesn't load eveyrthing in properly. May need to adjust con.undirected in future.

    # Remove the ith row and column in centres, tract_lengths and weights. i.e. the specified region(s)
    con.centres = np.delete(con.centres,ParamsDict["REMOVE"])
    con.weights = np.delete(con.weights,obj=ParamsDict["REMOVE"],axis=0)
    con.weights = np.delete(con.weights,obj=ParamsDict["REMOVE"],axis=1)
    con.tract_lengths = np.delete(con.tract_lengths,obj=ParamsDict["REMOVE"],axis=0)
    con.tract_lengths = np.delete(con.tract_lengths,obj=ParamsDict["REMOVE"],axis=1)

    # Number of regions
    con.number_of_regions = con.weights.shape[0]

    # Change to Connectome to Binary if desired:
    if ParamsDict["BINARY"]==True:
        con.weights = con.weights!=0

    # Set the parameter of the resting state simulation

    # Bold Simulation, Initial Conditions Not Provided
    if (ParamsDict["BOLD"] == True and ParamsDict["Init_Cons"] == False):
        sim = simulator.Simulator(model=ParamsDict["MODEL"],
                                connectivity=con,
                                coupling=coupling.Linear(a=ParamsDict["G"]),
                                integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                            random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                                monitors=(monitors.Bold(period=1e3),
                                        monitors.TemporalAverage(period=1e3)),
                                simulation_length=ParamsDict["Simul_length"],
                                #initial_conditions=0.5 + numpy.zeros((con.number_of_regions*con.number_of_regions,2,con.number_of_regions,1)),
                                ).configure()
        # Run the resting state simulation
        (bold_time, bold_data), _ = sim.run()

    # Bold Simulation, Initial Conditions Provided
    elif (ParamsDict["BOLD"] == True and ParamsDict["Init_Cons"] != False):
        sim = simulator.Simulator(model=ParamsDict["MODEL"],
                                connectivity=con,
                                coupling=coupling.Linear(a=ParamsDict["G"]),
                                integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                            random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                                monitors=(monitors.Bold(period=1e3),
                                        monitors.TemporalAverage(period=1e3)),
                                simulation_length=ParamsDict["Simul_length"],
                                initial_conditions=ParamsDict["Init_Cons"],
                                ).configure()
        # Run the resting state simulation
        (bold_time, bold_data), _ = sim.run()

    # No Monitors, Initial Conditions Not Provided
    elif (ParamsDict["BOLD"] == False and ParamsDict["Init_Cons"] == False):
        sim = simulator.Simulator(model=ParamsDict["MODEL"],
                                connectivity=con,
                                coupling=coupling.Linear(a=ParamsDict["G"]),
                                integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                                random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                                simulation_length=ParamsDict["Simul_length"],
                                #initial_conditions=0.5 + numpy.zeros((con.number_of_regions*con.number_of_regions,2,con.number_of_regions,1)),
                                ).configure()
        # Run the resting state simulation
        awer = sim.run()
        bold_time = awer[0][0]
        bold_data = awer[0][1]

    # No Monitors, Initial Conditions Provided
    else:
        sim = simulator.Simulator(model=ParamsDict["MODEL"],
                                connectivity=con,
                                coupling=coupling.Linear(a=ParamsDict["G"]),
                                integrator=integrators.EulerStochastic(dt=ParamsDict["dt"],noise=noise.Additive(nsig=ParamsDict["noise"],
                                                random_stream=np.random.RandomState(ParamsDict["RandState"]))),
                                simulation_length=ParamsDict["Simul_length"],
                                initial_conditions=ParamsDict["Init_Cons"],
                                ).configure()
        # Run the resting state simulation
        awer = sim.run()
        bold_time = awer[0][0]
        bold_data = awer[0][1]


    # Got lazy
    Snip = ParamsDict["Snip"]
    # Note, if BOLD = False, Snip gets multiplied by 100. 
    if ParamsDict["BOLD"] == False:
        Snip = 100 * Snip

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

    # Comparing SC vs FC with Spearman Corr  (Known as SCFC)
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
    Scorr = stats.spearmanr(a=FCM_Upper,b=SCM_Upper)
    print(Scorr)

# Calculate FC-FC score for mouse if requested.
    if ParamsDict["FCFC"] == True:
        # FCM_exp
        FCM_exp = np.genfromtxt('FCM_MouseExperimental.csv',delimiter = "\t")
        # Set diagonals to NaN
        np.fill_diagonal(FCM_exp,np.nan)

        # Remove the ith row and column in FCM (i.e. the specified region)
        FCM_exp = np.delete(FCM_exp,obj=ParamsDict["REMOVE"],axis=0)
        FCM_exp = np.delete(FCM_exp,obj=ParamsDict["REMOVE"],axis=1)

        # Comparing FC_experimental Vs FC_Simulation with Spearman Correlation

        FCM_Exp_U = FCM_exp[np.triu_indices(FCM_exp.shape[0], k = 1)]  
        FCM_Upper = FCM[np.triu_indices(FCM.shape[0], k = 1)]

        # FC-FC Spearman Correlation
        FCFC = stats.spearmanr(a=FCM_Exp_U,b=FCM_Upper)
        print(FCFC)
        # Concancatanate to end of Scorr output file. 
        Scorr = Scorr + FCFC

    # Export the simulation if requested

    if ParamsDict["ExportSim"] == True:
        time_now = time.strftime("%Y%m%d-%H%M%S")
        date = time.strftime("%Y_%m_%d/")
        # Create new directory which is the date. 
        os.makedirs("do-not-track/" + date,exist_ok=True)

        # Params Dictionary - Note how we sort the dictionary.
        with open("do-not-track/" + date + ParamsDict["tag"] + "_" + ParamsDict["name"] + "_Params_" + time_now + "_.csv", "w") as outfile:
            writer = csv.writer(outfile)
            for key, val in sorted(ParamsDict.items()):
                writer.writerow([key, val])
            
        # Create Time Series and save. 
        TSeries = np.concatenate((bold_time[Snip:].reshape(1,len(bold_time[Snip:])),TSeriesMatrix))
        np.savetxt("do-not-track/" + date + ParamsDict["tag"] + "_" + ParamsDict["name"] + "_Tseries_" + time_now + "_.csv", TSeries, delimiter="\t")
        np.savetxt("do-not-track/" + date + ParamsDict["tag"] + "_" + ParamsDict["name"]  + "_FCM_" + time_now + "_.csv", FCM, delimiter = "\t")
        np.savetxt("do-not-track/" + date + ParamsDict["tag"] + "_" + ParamsDict["name"]  + "_Scorr_" +  time_now + "_.csv", Scorr, delimiter = "\t")  

    return Scorr

