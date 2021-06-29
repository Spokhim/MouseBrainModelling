# This Python File contains all the functions to generate the appropriate sub-figures in the paper.
# Run Reproducibility.ipynb to see all the figures and how they're generated (before they are edited with inkscape for use in publication)

from tvb.simulator.coupling import Scaling
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
from pprint import pprint
import re

from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *
from matplotlib.colors import ListedColormap
from turbo_colormap import *

from scipy.io import loadmat

import sys
import glob
import pandas as pd
from pprint import pprint
import scipy.cluster.hierarchy as hierarchy

# Input Simulation Pipeline
from SimulationPipeline import *
from useful_fns import *

# Odeint and fsolve
from scipy.integrate import odeint
from scipy.optimize import fsolve

import random
import seaborn as sns
import matplotlib.lines as mlines

index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])

###  Figure 1: 



def Show_SCM():
    # Fig 1. b) - SCM 
    # Originally from Experimental Mouse Data Analysis.ipynb

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
    SCM = 1/SCM
    # Rearrange SCM
    index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])
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
    # Originally from Experimental Mouse Data Analysis.ipynb

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

def Het_Sample_Figures_Fig1(): 

    Regime = "LCycle"
    G_value = 0.7
    sig_e = 0.2
    sig_i = 0.2
    B_e_value = 1.5

    File_start = r"do-not-track\\2021_01_25\\" + Regime +"*_G[[]" + str(G_value) + "[]]*sig_e" + str(sig_e) + "_sig_i" + str(sig_i) # 1st Het Gradient - LCycleCutOld - 1.2e5

    TseriesFile = glob.glob(File_start+"*Tseries*_.csv")[0]
    ScorrFile = glob.glob(File_start + "*SCorr*_.csv")[0]
    FCMFile = glob.glob(File_start + "*FCM*_.csv")[0]
    ParamsFile = glob.glob(File_start + "*Params*_.csv")[0]

    # Read in Params - Unfortuantely doesn't do dtype conversion as well
    with open(ParamsFile, mode='r') as infile:
        reader = csv.reader(infile)
        ParamsDict = {rows[0]:rows[1] for rows in reader}

    ParamsDict["REMOVE"] = eval(ParamsDict["REMOVE"])
    ParamsDict["BINARY"] = eval(ParamsDict["BINARY"])
    G_value = eval(ParamsDict["G"])[0]
    '''
    # Empty dict
    ParamsDict = { }
    ParamsDict["name"] = "MouseCortex"
    ParamsDict["G"] = np.array([G_value]) 
    ParamsDict["REMOVE"] = [7]
    ParamsDict["BINARY"]=True
    '''

    # Read file import data
    #df = pd.read_csv(all_files[11],delimiter="\t",header=None)
    # Genfromtxt gives us a np array. 
    df = np.genfromtxt(TseriesFile,delimiter="\t")

    bold_time = df[0]
    bold_data = df[1:]

    # Re-arrange the order to a new order, Ben's Functional Grouping
    index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])
    bold_data = bold_data[index_fg]

    #end = 120000
    #bold_time = df[0,len(df[0])-120000:]
    #bold_data = df[1:,len(df[0])-120000:]

    # plt.subplots()

    #plt.figure(num=None, figsize=(60, 30), dpi=80, facecolor='w', edgecolor='k')
    # Let's only look at last 1000ms
    for tseries in bold_data:
        #plt.plot(bold_time[10000:11000],tseries[10000:11000])
        plt.plot(bold_time[len(bold_time)-1000:len(bold_time)],tseries[len(bold_time)-1000:len(bold_time)])
        #plt.plot(bold_time[0:1000],tseries[0:1000])

    plt.xlabel('Time (ms)', fontsize=20)
    plt.ylabel('Amplitude (au)', fontsize=20)
    plt.title('Simulated timeseries', fontsize=20)
    #plt.legend(('0','1','2','3','4'))
    plt.legend(range(38))
    plt.show()

    Scorra = np.genfromtxt(ScorrFile)
    print(Scorra)

    # Plot Simulated FCM
    FCM_sim = np.genfromtxt(FCMFile,delimiter="\t")
    # Re-order
    FCM_sim = FCM_sim[index_fg] [:,index_fg]

    # ListedColormap(turbo_colormap_data)
    cs=plt.imshow(FCM_sim, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none')
    #plt.title('Functional connectivity matrix', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('Correlation', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Node",fontsize=20)
    plt.ylabel("Node",fontsize=20)
    plt.show()

    FCM_exp = np.genfromtxt('FCM_MouseExperimental.csv',delimiter = "\t")
    # Set diagonals to NaN
    np.fill_diagonal(FCM_exp,np.nan)

    # Remove the ith row and column in FCM (i.e. the specified region)
    FCM_exp = np.delete(FCM_exp,obj=ParamsDict["REMOVE"],axis=0)
    FCM_exp = np.delete(FCM_exp,obj=ParamsDict["REMOVE"],axis=1)

    #Re_order
    FCM_exp = FCM_exp[index_fg] [:,index_fg]

    FCM_Exp_U = FCM_exp[np.triu_indices(FCM_exp.shape[0], k = 1)]
    FCM_Sim_U = FCM_sim[np.triu_indices(FCM_sim.shape[0], k = 1)]

    # Spearman Correlation
    SCorr = stats.spearmanr(a=FCM_Exp_U,b=FCM_Sim_U)
    #print(SCorr)

    # Scatterplot
    plt.scatter(FCM_Exp_U,FCM_Sim_U)
    plt.xlabel('$FC_{exp}$', fontsize=20)
    plt.ylabel('$FC_{sim}$', fontsize=20)
    #plt.title('$FCM_{sim}$ vs $FCM_{exp}$', fontsize=20)
    plt.annotate("FCFC = " + str(round(SCorr[0],2)), xy=(0.05, 0.95), xycoords='axes fraction')
    plt.show()

    return

def Homog_sweep(Path):
    # Originally from HomogSweep.ipynb

    SCorr_files = glob.glob(Path)

    Params = []
    SCFC = []
    FCFC = []
    for string in SCorr_files:

        # Obtain Parameter Values
        x = re.findall("\[(.*)\].*\[(.*)\]",string)
        #x = re.findall("\[(.*)\].b_e(...)",string)
        Params.append(x[0])

        # FCFC and SCFC
        a = np.genfromtxt(string)
        SCFC.append(a[0])
        FCFC.append(a[2])

    df = pd.DataFrame(Params)
    df.columns = ['G', 'B_e']
    df["SCFC"] = SCFC
    df["FCFC"] = FCFC

    df = df.sort_values(by=["FCFC"]) 
    df_pivot = df.sort_values('FCFC').drop_duplicates(['G','B_e'],keep='last').pivot(index='B_e', columns='G', values='FCFC')
    # Reason why we need to do this is because the clsuter did something bad and simulated a few points multiple times (this is likely from it deciding to stop runs midway and then restarting)

    x = df_pivot.index[::5]
    y = x.astype(np.float)
    X = df_pivot.columns[::10]
    Y = X.astype(np.float)

    #ListedColormap(turbo_colormap_data) ,vmin=0,vmax=0.65
    cs=plt.imshow(df_pivot, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none',origin='lower',vmin=0,vmax=0.6,)
    axcb=plt.colorbar(cs)
    axcb.set_label('FCFC', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.yticks(ticks=np.arange(len(df_pivot.index))[::5],labels=y,fontsize=16)
    plt.xticks(ticks=np.arange(len(df_pivot.columns))[::10],labels=Y,fontsize=16)
    plt.xlabel("$G$", fontsize=20)
    plt.ylabel("$B_e$ (mV)", fontsize=20)

    plt.show()

    return

def Bif_Diagram_Generation():
    # Originally from Code in Progress.ipynb. 

    def sig(v):
    # Sigmoid function
        return 1/(1 + np.exp(-v))

    def func(z,J_e,J_i=0,b_e=5,b_i=3.7,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=10,w_ii=3,):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = -E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e))
        f[1] = -I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i))

        return f


    # Hysteresis 

    J_e = np.arange(0,5,0.05)
    array = []
    array_3 = []

    for i in J_e:
        
        array_2 = []

        for j in np.arange(100):
            # Brute force try random ones:
            sol = fsolve(func,[random.random()/2,random.random()/2],i,maxfev=100000)
            # func(root) should be almost 0.0.
            if all(np.isclose(func(sol,i), [0.0, 0.0])):
                # Get only E value
                array_2.append(sol[0])
            else:
                array_2.append(np.nan)    

        # Then add the median, max, and min to the end of array

        if i < 2.35:
            MAX = np.nan
            MIN = np.nanmin(array_2)
            #MID = np.median(array_2)

        if i > 1.1:
            MAX = np.nanmax(array_2)

        if i > 2:
            MIN= np.nan

        '''
        for k in array_2:
            if (np.isclose(k,MAX)):
                #Do nothing
                dummy = 1
            else:
                MID = k

        '''
        array.append([MAX,MIN])
        array_3.append(array_2)

    array = np.array(array)
    array_3 = np.array(array_3)

    super_threshold_indices = (array_3 < 0.07) 
    array_3[super_threshold_indices] = np.nan
    super_threshold_indices = (array_3 > 0.28)
    array_3[super_threshold_indices] = np.nan

    plt.plot(J_e,array[:,0],'k',)
    plt.plot(J_e,array_3,'k.',markersize=3)
    plt.plot(J_e,array[:,1],'k')
    plt.legend(('Stable','Unstable'),title="Equilibrium",loc="lower right",title_fontsize="x-large",fontsize="x-large")

    plt.vlines(x=0,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(0,0.5,"$B_e$=5",fontsize=14)

    plt.vlines(x=1,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(1,0.5,"$B_e$=4",fontsize=14)

    plt.vlines(x=2,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(2,0.5,"$B_e$=3",fontsize=14)

    #plt.vlines(x=3,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    #plt.text(3,0.5,"$B_e$=0",fontsize=14)

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )

    plt.show()


    # Define the Wilson Cowan Equations

    def sig(v):
        # Sigmoid function
        return 1/(1 + np.exp(-v))

    # Fixed pt - J_i=0,b_e=3,b_i=3.7,a_e=1,a_i=1,w_ee=12,w_ei=15,w_ie=10,w_ii=8,
    # LCycle - J_i=0,b_e=2,b_i=2.8,a_e=1,a_i=1,w_ee=11,w_ei=10,w_ie=10,w_ii=1,
    # LCycleH - J_i=0,b_e=3,b_i=4,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=15,w_ii=3,

    def func(z,t,J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,_tau_i):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = (-E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e)))/tau_e
        f[1] = (-I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i)))/tau_i

        return f

    # Fixed Pt

    The_max = []
    The_min = []

    # Parameters

    # Fixed Pt
    w_ee=12
    w_ei=15
    w_ie=10
    w_ii=8
    b_e=5
    b_i=4
    tau_e=10
    tau_i=10
    a_e=1
    a_i=1
    #J_e=3.5
    J_i=0

    # Integrator Settings
    length = 10000 
    dt = 0.1
    t = np.arange(0, length, dt)
    J = np.arange(0,5,0.05)

    for J_e in J:
        # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
        solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

        # Obtain the max and min of the last 1/10 (tenth) elements
        eqbm_max = max(solut[-int(length/dt/10):,0])
        eqbm_min = min(solut[-int(length/dt/10):,0])

        The_max.append(eqbm_max)
        The_min.append(eqbm_min)

    plt.plot(J,The_max,color='k')
    plt.plot(J,The_min,color='k')

    plt.vlines(x=0,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(0,0.5,"$B_e$=5",fontsize=14)

    plt.vlines(x=1,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(1,0.5,"$B_e$=4",fontsize=14)

    plt.vlines(x=2,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(2,0.5,"$B_e$=3",fontsize=14)

    plt.vlines(x=3,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(3,0.5,"$B_e$=2",fontsize=14)

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )

    plt.show()

    # LCycle 

    The_max = []
    The_min = []

    # LCycle
    w_ee=11
    w_ei=10
    w_ie=10
    w_ii=1
    b_e=3
    b_i=2.8
    tau_e=10
    tau_i=65
    a_e=1
    a_i=1
    #J_e=3.5
    J_i=0

    # Integrator Settings
    length = 10000 
    dt = 0.1
    t = np.arange(0, length, dt)
    J = np.arange(0,5,0.05)

    for J_e in J:
        # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
        solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

        # Obtain the max and min of the last 1/10 (tenth) elements
        eqbm_max = max(solut[-int(length/dt/10):,0])
        eqbm_min = min(solut[-int(length/dt/10):,0])

        The_max.append(eqbm_max)
        The_min.append(eqbm_min)



    plt.plot(J,The_max,color='k')
    plt.plot(J,The_min,color='k')

    plt.vlines(x=0,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(0,0.5,"$B_e$=3",fontsize=14)

    plt.vlines(x=1,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(1,0.5,"$B_e$=2",fontsize=14)

    plt.vlines(x=2,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(2,0.5,"$B_e$=1",fontsize=14)

    plt.vlines(x=3,ymin=0,ymax=0.5,colors='k',linestyles='dashed')
    plt.text(3,0.5,"$B_e$=0",fontsize=14)

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )

    plt.show()

    return

def Single_Run_Plots(File_start,Regime,G_value,B_e_value):

    # Originally from SingleRunAnalysis.ipynb

    TseriesFile = glob.glob(File_start+"*Tseries*_.csv")[0]
    ScorrFile = glob.glob(File_start + "*SCorr*_.csv")[0]
    FCMFile = glob.glob(File_start + "*FCM*_.csv")[0]
    ParamsFile = glob.glob(File_start + "*Params*_.csv")[0]

    # Read in Params - Unfortuantely doesn't do dtype conversion as well
    with open(ParamsFile, mode='r') as infile:
        reader = csv.reader(infile)
        ParamsDict = {rows[0]:rows[1] for rows in reader}

    ParamsDict["REMOVE"] = eval(ParamsDict["REMOVE"])
    ParamsDict["BINARY"] = eval(ParamsDict["BINARY"])

    # Read file import data
    #df = pd.read_csv(all_files[11],delimiter="\t",header=None)
    # Genfromtxt gives us a np array. 
    df = np.genfromtxt(TseriesFile,delimiter="\t")

    bold_time = df[0]
    bold_data = df[1:]

    # Re-arrange the order to a new order, Ben's Functional Grouping
    index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])
    bold_data = bold_data[index_fg]

    # Look at specific region:

    time_region_start = len(bold_time)-1000
    time_region_end = len(bold_time)
    #time_region_start = 45200
    #time_region_end = 46200

    for i in np.arange(6):    
        Region = i + 8
        plt.plot(bold_time[time_region_start:time_region_end],bold_data[Region][time_region_start:time_region_end])

    plt.xlabel('Time (s)', fontsize=20)
    plt.xticks(ticks=np.arange(start=119000,stop=1.202e5,step=200),labels=('0.0','0.2','0.4','0.6','0.8','1.0'))
    plt.ylabel('$E$ (au)', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(np.arange(6)+8,title='Node',title_fontsize=14,fontsize=14,loc='lower right')
    plt.show()

    # Get rough feel for external currents.
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

    # Now-rearrange SCM according to new order
    SCM = SCM[index_fg] [:,index_fg]

    # Takes roughly 1 min.

    # External Current Calculator:
    J_e = []
    #len(bold_time)
    # j is jth element
    for j in np.arange(len(bold_time)):       
        t_0 = []
        # Specific column (or time point)
        for i in np.arange(SCM.shape[0]): 
            # Sum over all external currents (May need to do SCM[:,i] instead, but quite sure SCM[i,:] is correct from tseries analysis)  
            t  = sum(bold_data[:,j]*SCM[i,:])*G_value
            # To obtain currents to particular region
            t_0.append(t)
        J_e.append(t_0)

    J_e = np.array(J_e)
    print("Max",np.max(J_e))
    print("Median",np.median(J_e))
    print("Min",np.min(J_e))

    df_Je = pd.DataFrame(J_e)

    # Define the Wilson Cowan Equations

    def sig(v):
        # Sigmoid function
        return 1/(1 + np.exp(-v))

    # Fixed pt - J_i=0,b_e=3,b_i=3.7,a_e=1,a_i=1,w_ee=12,w_ei=15,w_ie=10,w_ii=8,
    # LCycle - J_i=0,b_e=2,b_i=2.8,a_e=1,a_i=1,w_ee=11,w_ei=10,w_ie=10,w_ii=1,
    # LCycleH - J_i=0,b_e=3,b_i=4,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=15,w_ii=3,

    def func(z,t,J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,_tau_i):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = (-E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e)))/tau_e
        f[1] = (-I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i)))/tau_i

        return f

    The_max = []
    The_min = []

    # Parameters 

    if Regime == "Hysteresis":
        # Hysteresis
        w_ee=16
        w_ei=12
        w_ie=10
        w_ii=3
        b_e=B_e_value
        #b_e = 2.5
        b_i=3.7
        tau_e=10
        tau_i=10
        a_e=1.3
        a_i=2
        #J_e=3.5
        J_i=0


    elif Regime == "FixedPt": 
        # Fixed Pt
        w_ee=12
        w_ei=15
        w_ie=10
        w_ii=8
        b_e=B_e_value
        #b_e = 2.5
        b_i=4
        tau_e=10
        tau_i=10
        a_e=1
        a_i=1
        #J_e=3.5
        J_i=0

    else:    
        # LCycle
        w_ee=11
        w_ei=10
        w_ie=10
        w_ii=1
        b_e=B_e_value
        b_i=2.8
        tau_e=10
        tau_i=65
        a_e=1
        a_i=1
        #J_e=3.5
        J_i=0

    # Integrator Settings
    length = 10000 
    dt = 0.1
    t = np.arange(0, length, dt)
    J = np.arange(0,5,0.05)

    for J_e in J:
        # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
        solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

        # Obtain the max and min of the last 1/10 (tenth) elements
        eqbm_max = max(solut[-int(length/dt/10):,0])
        eqbm_min = min(solut[-int(length/dt/10):,0])

        The_max.append(eqbm_max)
        The_min.append(eqbm_min)

    # Hysteresis needs a second run for the other branch:
    if Regime =="Hysteresis":
        The_max2 = []
        The_min2 = []

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,5,0.05)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0, 1], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max2.append(eqbm_max)
            The_min2.append(eqbm_min)

    plt.plot(J,The_max,color='k')

    if Regime =="Hysteresis":
        plt.plot(J,The_max2,color='k')
        plt.plot(J,The_min2,color='k')

    colour_list = ('pink','pink','pink','pink','pink','pink','pink','pink','lightblue','lightblue','lightblue','lightblue','lightblue','lightblue','gold','gold','gold','gold','gold','gold','plum','plum','plum','plum','darkorange','darkorange','darkorange','darkorange','darkorange','green','green','green','green','green','green','green','green')

    # Get min-max external currents & plot on the diagram
    for i in np.arange(37):
        # Region i
        plt.plot((df_Je.min(axis=0)[i],df_Je.max(axis=0)[i]),np.zeros(2)+0.5/37*i,'-o',color=colour_list[i])
        region_label = str(i) #+ "-" + acronyms[i]
        #plt.text(df_Je.min(axis=0)[i],0.5/37*i,region_label,fontsize=14)
        
    # The two branches of bifurcation
    plt.plot(J,The_max,color='k')
    plt.plot(J,The_min,color='k')

    #region_label = str(12) + "-" + acronyms[12]
    #plt.text(df_Je.min(axis=0)[12],0.5/37*12,region_label,fontsize=14)

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    #plt.title("Bifurcation Diagram - $E$ against $J_e$", fontsize=20)

    Eqbm_Pos= mlines.Line2D([], [], color='black', markersize=15, label='Equilibrium Position')
    Somatomotor = mlines.Line2D([], [], color='pink', marker='o', label='Somatomotor')
    Medial = mlines.Line2D([], [], color='lightblue', marker='o', label='Medial')
    Temporal = mlines.Line2D([], [], color='gold', marker='o', label='Temporal')
    Visual = mlines.Line2D([], [], color='plum', marker='o', label='Visual')
    Anterolateral = mlines.Line2D([], [], color='darkorange', marker='o', label='Anterolateral')
    Prefrontal = mlines.Line2D([], [], color='green', marker='o', label='Prefrontal')
                            
    plt.legend(handles=[Eqbm_Pos,Prefrontal,Anterolateral,Visual,Temporal,Medial,Somatomotor],loc="lower right",fontsize=14)
    #plt.legend(("Equilibrium Position","Somatomotor","Medial","Temporal","Visual","Anterolateral","Prefrontal"), loc="lower right",fontsize=14)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()

    Scorra = np.genfromtxt(ScorrFile)
    print(Scorra)

    # Plot Simulated FCM
    FCM_sim = np.genfromtxt(FCMFile,delimiter="\t")

    # Re-order
    FCM_sim = FCM_sim[index_fg] [:,index_fg]

    # ListedColormap(turbo_colormap_data)
    cs=plt.imshow(FCM_sim, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none')
    #plt.title('Functional connectivity matrix', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('Correlation', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Node",fontsize=20)
    plt.ylabel("Node",fontsize=20)
    plt.show()

    # Unsorted - Square Carpet

    # ListedColormap(turbo_colormap_data) ListedColormap(turbo_colormap_data) [:,::40]
    cs=plt.imshow(bold_data, cmap="gray", aspect='equal',vmin=0,vmax=0.6, interpolation='none',origin='lower')
    #plt.title('Timeseries data - $E$ Unsorted', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('$E$ (au)', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.ylabel('Node', fontsize=20)
    plt.xlabel('Time (s)', fontsize=20)
    plt.xticks(ticks=np.arange(start=0,stop=1.2e5,step=2e4),labels=np.arange(start=0,stop=120,step=20))
    plt.axis("tight")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()

    return

def LCycle_Het_Bif_Diagrams(Regions,Scaling):
# Taken fron Bifurcation_Diagrams_Analysis.ipynb. 


    # Define the Wilson Cowan Equations

    def sig(v):
        # Sigmoid function
        return 1/(1 + np.exp(-v))

    # Fixed pt - J_i=0,b_e=3,b_i=3.7,a_e=1,a_i=1,w_ee=12,w_ei=15,w_ie=10,w_ii=8,
    # LCycle - J_i=0,b_e=2,b_i=2.8,a_e=1,a_i=1,w_ee=11,w_ei=10,w_ie=10,w_ii=1,
    # LCycleH - J_i=0,b_e=3,b_i=4,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=15,w_ii=3,

    def func(z,t,J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,_tau_i):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = (-E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e)))/tau_e
        f[1] = (-I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i)))/tau_i

        return f

    # Set parameters:
    # LCycle
    W_ee=11
    W_ei=10
    W_ie=10
    W_ii=1
    b_e=3
    b_i=2.8
    tau_e=10
    tau_i=65
    a_e=1
    a_i=1
    #J_e=3.5
    J_i=0

    # Homogeneous Case:

    The_max = []
    The_min = []

    # Parameters

    w_ee=W_ee
    w_ei=W_ei
    w_ie=W_ie
    w_ii=W_ii

    # Integrator Settings
    length = 10000 
    dt = 0.1
    t = np.arange(0, length, dt)
    J = np.arange(0,5,0.05)

    for J_e in J:
        # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
        solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

        # Obtain the max and min of the last 1/10 (tenth) elements
        eqbm_max = max(solut[-int(length/dt/10):,0])
        eqbm_min = min(solut[-int(length/dt/10):,0])

        The_max.append(eqbm_max)
        The_min.append(eqbm_min)

    Homog_max = The_max
    Homog_min = The_min

    # Scaling Term:
    R_e = np.linspace(-Scaling,Scaling,num=Regions)
    #R_i = -np.linspace(-0.1,0.1,num=Regions)
    R_i = np.zeros(Regions)

    h_c_ee =  W_ee * (1 + R_e) 
    h_c_ei = W_ei * (1 + R_i) 
    h_c_ie = W_ie * (1 + R_e) 
    h_c_ii = W_ii * (1 + R_i) 

    # Simulate Het BIfurcation Diagrams

    AllRuns_max = []
    AllRuns_min = []

    for i in np.arange(len(h_c_ee)):
        The_max = []
        The_min = []

        # Parameters
        w_ee=h_c_ee[i]
        w_ei=h_c_ei[i]
        w_ie=h_c_ie[i]
        w_ii=h_c_ii[i]

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,5,0.05)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max.append(eqbm_max)
            The_min.append(eqbm_min)

        AllRuns_max.append(The_max)
        AllRuns_min.append(The_min)
    # Plot.

    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    plt.plot(J,Homog_max,color='k')
    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_max[i],color=colors[i],)   

    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_min[i],color=colors[i],)

    plt.plot(J,Homog_min,color='k')  

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    plt.title("$R_e$ only", fontsize=20)
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")
    plt.show()

    # R_e = 0 

    # Scaling Term:
    #R_e = np.linspace(-0.1,0.1,num=Regions)
    R_i = np.linspace(-Scaling,Scaling,num=Regions)
    R_e = np.zeros(Regions)

    h_c_ee =  W_ee * (1 + R_e) 
    h_c_ei = W_ei * (1 + R_i) 
    h_c_ie = W_ie * (1 + R_e) 
    h_c_ii = W_ii * (1 + R_i) 

    # Simulate Het BIfurcation Diagrams

    AllRuns_max = []
    AllRuns_min = []

    for i in np.arange(len(h_c_ee)):
        The_max = []
        The_min = []

        # Parameters
        w_ee=h_c_ee[i]
        w_ei=h_c_ei[i]
        w_ie=h_c_ie[i]
        w_ii=h_c_ii[i]

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,5,0.05)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max.append(eqbm_max)
            The_min.append(eqbm_min)

        AllRuns_max.append(The_max)
        AllRuns_min.append(The_min)

    # Plot.

    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    plt.plot(J,Homog_max,color='k')
    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_max[i],color=colors[i],)   

    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_min[i],color=colors[i],)



    plt.plot(J,Homog_min,color='k')  

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    plt.title("$R_i$ only", fontsize=20)
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")
    plt.show()

    # R_e = R_i

    # Scaling Term:
    R_e = np.linspace(-Scaling,Scaling,num=Regions)
    R_i = np.linspace(-Scaling,Scaling,num=Regions)
    #R_e = np.zeros(Regions)

    h_c_ee =  W_ee * (1 + R_e) 
    h_c_ei = W_ei * (1 + R_i) 
    h_c_ie = W_ie * (1 + R_e) 
    h_c_ii = W_ii * (1 + R_i) 

    # Simulate Het BIfurcation Diagrams

    AllRuns_max = []
    AllRuns_min = []

    for i in np.arange(len(h_c_ee)):
        The_max = []
        The_min = []

        # Parameters
        w_ee=h_c_ee[i]
        w_ei=h_c_ei[i]
        w_ie=h_c_ie[i]
        w_ii=h_c_ii[i]

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,5,0.05)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max.append(eqbm_max)
            The_min.append(eqbm_min)

        AllRuns_max.append(The_max)
        AllRuns_min.append(The_min)

    # Plot.

    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    plt.plot(J,Homog_max,color='k')
    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_max[i],color=colors[i],)   

    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_min[i],color=colors[i],)



    plt.plot(J,Homog_min,color='k')  

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    plt.title("$R_e = R_i$", fontsize=20)
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")
    plt.show()

    # R_e = -R_i

    # Scaling Term:
    R_e = np.linspace(-Scaling,Scaling,num=Regions)
    R_i = -np.linspace(-Scaling,Scaling,num=Regions)
    #R_e = np.zeros(Regions)

    h_c_ee =  W_ee * (1 + R_e) 
    h_c_ei = W_ei * (1 + R_i) 
    h_c_ie = W_ie * (1 + R_e) 
    h_c_ii = W_ii * (1 + R_i) 

    # Simulate Het BIfurcation Diagrams

    AllRuns_max = []
    AllRuns_min = []

    for i in np.arange(len(h_c_ee)):
        The_max = []
        The_min = []

        # Parameters
        w_ee=h_c_ee[i]
        w_ei=h_c_ei[i]
        w_ie=h_c_ie[i]
        w_ii=h_c_ii[i]

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,5,0.05)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max.append(eqbm_max)
            The_min.append(eqbm_min)

        AllRuns_max.append(The_max)
        AllRuns_min.append(The_min)

    # Plot.

    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    plt.plot(J,Homog_max,color='k')
    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_max[i],color=colors[i],)   

    for i in np.arange(len(h_c_ee)):
        plt.plot(J,AllRuns_min[i],color=colors[i],)



    plt.plot(J,Homog_min,color='k')  

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    plt.title("$R_e = -R_i$", fontsize=20)
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")
    plt.show()

    return

def Single_Run_Het_Plots(File_start,Regime,G_value,B_e_value,sig_e,sig_i): 
# Originally from SingleRun_Het_Analysis.ipynb


    TseriesFile = glob.glob(File_start+"*Tseries*_.csv")[0]
    ScorrFile = glob.glob(File_start + "*SCorr*_.csv")[0]
    FCMFile = glob.glob(File_start + "*FCM*_.csv")[0]
    ParamsFile = glob.glob(File_start + "*Params*_.csv")[0]

    # Read in Params - Unfortuantely doesn't do dtype conversion as well
    with open(ParamsFile, mode='r') as infile:
        reader = csv.reader(infile)
        ParamsDict = {rows[0]:rows[1] for rows in reader}

    ParamsDict["REMOVE"] = eval(ParamsDict["REMOVE"])
    ParamsDict["BINARY"] = eval(ParamsDict["BINARY"])
    G_value = eval(ParamsDict["G"])[0]
    '''
    # Empty dict
    ParamsDict = { }
    ParamsDict["name"] = "MouseCortex"
    ParamsDict["G"] = np.array([G_value]) 
    ParamsDict["REMOVE"] = [7]
    ParamsDict["BINARY"]=True
    '''

    # Read file import data
    #df = pd.read_csv(all_files[11],delimiter="\t",header=None)
    # Genfromtxt gives us a np array. 
    df = np.genfromtxt(TseriesFile,delimiter="\t")

    bold_time = df[0]
    bold_data = df[1:]

    # Re-arrange the order to a new order, Ben's Functional Grouping
    index_fg = np.array([13,31,10,8,7,9,11,12,0,15,19,25,26,27,34,33,35,29,20,28,16,14,17,18,21,36,4,6,5,32,1,22,30,24,23,3,2])
    bold_data = bold_data[index_fg]

    # Get rough feel for external currents.
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

    # Now-rearrange SCM according to new order
    '''
    # For some reason, this method of re-arranging is faulty when calcualting J_e's
    new_order = index_fg
    SCM[:] = [SCM[i] for i in new_order]
    for row in SCM:
        row[:] = [row[i] for i in new_order]
        '''
    SCM = SCM[index_fg] [:,index_fg]

    # Takes roughly 1 min.

    # External Current Calculator:
    J_e = []
    #len(bold_time)
    # j is jth element
    for j in np.arange(len(bold_time)):       
        t_0 = []
        # Specific column (or time point)
        for i in np.arange(SCM.shape[0]): 
            # Sum over all external currents (May need to do SCM[:,i] instead)  
            t  = sum(bold_data[:,j]*SCM[i,:])*G_value
            # To obtain currents to particular region
            t_0.append(t)
        J_e.append(t_0)

    J_e = np.array(J_e)
    print("Max",np.max(J_e))
    print("Median",np.median(J_e))
    print("Min",np.min(J_e))

    df_Je = pd.DataFrame(J_e)

    # Define the Wilson Cowan Equations

    def sig(v):
        # Sigmoid function
        return 1/(1 + np.exp(-v))

    # Fixed pt - J_i=0,b_e=3,b_i=3.7,a_e=1,a_i=1,w_ee=12,w_ei=15,w_ie=10,w_ii=8,
    # LCycle - J_i=0,b_e=2,b_i=2.8,a_e=1,a_i=1,w_ee=11,w_ei=10,w_ie=10,w_ii=1,
    # LCycleH - J_i=0,b_e=3,b_i=4,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=15,w_ii=3,

    def func(z,t,J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,_tau_i):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = (-E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e)))/tau_e
        f[1] = (-I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i)))/tau_i

        return f

    Scorra = np.genfromtxt(ScorrFile)
    print(Scorra)

    # Plot Simulated FCM
    FCM_sim = np.genfromtxt(FCMFile,delimiter="\t")
    # Re-order
    FCM_sim = FCM_sim[index_fg] [:,index_fg]

    # ListedColormap(turbo_colormap_data)
    cs=plt.imshow(FCM_sim, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none')
    #plt.title('Functional connectivity matrix', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('Correlation', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Node",fontsize=20)
    plt.ylabel("Node",fontsize=20)
    plt.show()

    # Unsorted
    # ListedColormap(turbo_colormap_data) ListedColormap(turbo_colormap_data) [:,::40]
    cs=plt.imshow(bold_data, cmap="gray", aspect='equal',vmin=0,vmax=0.6, interpolation='none',origin='lower')
    #plt.title('Timeseries data - $E$ Unsorted', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('$E$ (au)', fontsize=20)
    axcb.ax.tick_params(labelsize=16)
    plt.ylabel('Node', fontsize=20)
    plt.xlabel('Time (ms)', fontsize=20)
    plt.axis("tight")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.savefig("do-not-track\\Fig3.png",bbox_inches='tight')
    plt.show()

    # Bifurcation Diagram Generation

    The_max = []
    The_min = []

    # Parameters

    w_ee=11
    w_ei=10
    w_ie=10
    w_ii=1
    b_e=B_e_value
    b_i=2.8
    tau_e=10
    tau_i=65
    a_e=1
    a_i=1
    #J_e=3.5
    J_i=0

    # Integrator Settings
    length = 10000 
    dt = 0.1
    t = np.arange(0, length, dt)
    J = np.arange(0,6.5,0.025)

    for J_e in J:
        # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
        solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

        # Obtain the max and min of the last 1/10 (tenth) elements
        eqbm_max = max(solut[-int(length/dt/10):,0])
        eqbm_min = min(solut[-int(length/dt/10):,0])

        The_max.append(eqbm_max)
        The_min.append(eqbm_min)

    Homog_max = The_max
    Homog_min = The_min

    h_c_ee = np.fromstring("11.14799794 10.39054527 10.99995073 10.7799896  10.59026485 10.65829007 10.59197152 10.9799298  10.7706965  10.44587736 10.52463601 11.37952596  10.5086166  10.98251696 11.3933426  11.13850692 11.23532824 10.96687459 11.37414929 12.18638905 11.49092802 10.86727987 10.57059861 10.76835704 10.73324928 11.08686034 11.32023041 11.1347711  11.50488189 11.14964262 10.79733952 10.45870114 10.37042949 11.75960984 11.44570641 11.83789755 10.658117", dtype=float, sep=" ")

    h_c_ei = np.fromstring("9.60574234  9.96235553 10.72885439 11.38835088  9.52194115 10.29199435 9.38113738  9.12243791  9.16156425  9.2075798   9.2901189   9.61767677   9.0910805   9.34630171  9.99703063  9.84167024  9.44810341  9.13810887   9.97600759 12.87632019  9.86990012  9.59657777  9.88821828 10.95144449  10.66472672 10.63257872 11.99191439 10.20027733 10.75132571 10.18097831  10.61728726  9.40211826 10.04810382  9.65252972  9.49666513  9.77250279  9.28847438", dtype=float, sep=" ")

    h_c_ie = np.fromstring("10.13454359  9.44595024  9.99995521  9.79999055  9.6275135   9.68935461   9.62906502  9.98175436  9.79154227  9.49625214  9.56785091 10.3450236   9.55328782  9.98410633 10.35758419 10.12591538 10.21393476  9.96988599  10.34013572 11.0785355  10.4462982   9.87934533  9.6096351   9.78941549   9.75749934 10.07896395 10.29111856 10.12251918 10.45898354 10.13603875   9.8157632   9.50791012  9.42766317 10.6905544  10.40518765 10.76172505   9.68919727", dtype=float, sep=" ")

    h_c_ii = np.fromstring("0.96057423 0.99623555 1.07288544 1.13883509 0.95219411 1.02919944  0.93811374 0.91224379 0.91615642 0.92075798 0.92901189 0.96176768 0.90910805 0.93463017 0.99970306 0.98416702 0.94481034 0.91381089  0.99760076 1.28763202 0.98699001 0.95965778 0.98882183 1.09514445 1.06647267 1.06325787 1.19919144 1.02002773 1.07513257 1.01809783  1.06172873 0.94021183 1.00481038 0.96525297 0.94966651 0.97725028 0.92884744", dtype=float, sep=" ")

    # Need to re-order this hard-code for coupling parameters 
    # Need to fix df_Je:
    #df_Je = df_Je[index_fg]
    h_c_ee = h_c_ee[index_fg]
    h_c_ei = h_c_ei[index_fg]
    h_c_ie = h_c_ie[index_fg]
    h_c_ii = h_c_ii[index_fg]

    # Note this df_Je fix thing only works for the flattened version for now. 

    AllRuns_max = []
    AllRuns_min = []

    for i in np.arange(len(h_c_ee)):
        The_max = []
        The_min = []

        # Parameters
        w_ee=h_c_ee[i]
        w_ei=h_c_ei[i]
        w_ie=h_c_ie[i]
        w_ii=h_c_ii[i]
        b_e=B_e_value
        b_i=2.8
        tau_e=10
        tau_i=65
        a_e=1
        a_i=1
        #J_e=3.5
        J_i=0

        # Integrator Settings
        length = 10000 
        dt = 0.1
        t = np.arange(0, length, dt)
        J = np.arange(0,6.5,0.025)

        for J_e in J:
            # Solve it! Note that the additional "args" supplied to "odeint" must be in a tuple; "(a,)".
            solut = odeint(func, [0.5, 0.5], t, args=(J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,tau_i) )

            # Obtain the max and min of the last 1/10 (tenth) elements
            eqbm_max = max(solut[-int(length/dt/10):,0])
            eqbm_min = min(solut[-int(length/dt/10):,0])

            The_max.append(eqbm_max)
            The_min.append(eqbm_min)

        AllRuns_max.append(The_max)
        AllRuns_min.append(The_min)

        colour_list = ('pink','pink','pink','pink','pink','pink','pink','pink','lightblue','lightblue','lightblue','lightblue','lightblue','lightblue','gold','gold','gold','gold','gold','gold','plum','plum','plum','plum','darkorange','darkorange','darkorange','darkorange','darkorange','green','green','green','green','green','green','green','green')

    for i in np.arange(len(h_c_ee)):
        plt.plot(J[0:160],AllRuns_max[i][0:160],color=colour_list[i],label='_nolegend_')   

    for i in np.arange(len(h_c_ee)):
        plt.plot(J[0:160],AllRuns_min[i][0:160],color=colour_list[i],label='_nolegend_')

    plt.plot(J[0:160],Homog_max[0:160],color='k')
    #plt.plot(J,AllRuns_max[12],color='g') 
    #plt.plot(J,AllRuns_max[19],color='r')   
    #plt.plot(J,AllRuns_min[36],color='lightcoral',)

    # Get min-max external currents & plot on the diagram
    # Region 12
    #plt.plot((df_Je.min(axis=0)[12],df_Je.max(axis=0)[12]),np.zeros(2)+0.01,'g-o')
    #plt.text(df_Je.min(axis=0)[12],0.02,"12 - SSp-ul",fontsize=14)
    # Region 19 
    #plt.plot((df_Je.min(axis=0)[19],df_Je.max(axis=0)[19]),np.zeros(2)+0.02,'r-o')

    plt.plot(J[0:160],Homog_min[0:160],color='k')  
    #plt.plot(J,AllRuns_min[12],color='g')
    #plt.plot(J,AllRuns_min[19],color='r')

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au)", fontsize=20)
    #plt.title("Bifurcation Diagram - $E$ against $J_e$", fontsize=20)
    #plt.legend(("Homogeneous","Het Node 12","Het Node 19", "Het Other Nodes","Node 12 $J_e$ Range","Node 19 $J_e$ Range"),loc="lower right",fontsize="x-large")
    #plt.legend(("Homogeneous","Heterogeneous Regions"),loc="lower right",fontsize="x-large")
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    #plt.savefig("do-not-track\\Het_Bifurcation_FuncColoured.pdf",bbox_inches='tight')
    plt.show()

    jet= plt.get_cmap('summer')
    w_colors = jet(np.linspace(0,1,37))
    jet= plt.get_cmap('winter')
    c_colors = jet(np.linspace(0,1,37))

    for i in np.arange(len(h_c_ee)):
        # First of all, figure out what the current range is:
        J_max = df_Je.iloc[:,i].max()
        J_min = df_Je.iloc[:,i].min()

        # Next figure out, when it is in a limit cycle and when it is at a fixed pt. 

        # Fixed Pt when True, LCycle when False.:

        Index_FixPt = (np.array(AllRuns_max[i]) - np.array(AllRuns_min[i]) < 1e-10 ) & (J <= J_max) & (J >= J_min)

        Index_LCycle = (np.array(AllRuns_max[i]) - np.array(AllRuns_min[i]) > 1e-10 ) & (J <= J_max) & (J >= J_min)

        # Plot in Red if Fixed pt:
        plt.plot(J[Index_FixPt],i*np.ones(len(J[Index_FixPt])),color='r')

        # Plot in Blue if Lcycle:
        plt.plot(J[Index_LCycle],i*np.ones(len(J[Index_LCycle])),color='b')

    # Homogeneous case flattened:
    Index_H_FixPt = (np.array(Homog_max) - np.array(Homog_min) < 1e-10 )
    Index_H_LCycle = (np.array(Homog_max) - np.array(Homog_min) > 1e-10 )

    # Plot in Red if Fixed pt:
    plt.plot(J[Index_H_FixPt],38+np.zeros(len(J[Index_H_FixPt])), color = 'r', linestyle='dashed')
    # Plot in Blue if Lcycle:
    plt.plot(J[Index_H_LCycle],38+np.zeros(len(J[Index_H_LCycle])), color = 'b', linestyle='dashed')

    # Add background colours
    plt.axhspan(0, 7.5, facecolor='pink', alpha=0.5)
    plt.axhspan(7.5, 13.5, facecolor='lightblue', alpha=0.5)
    plt.axhspan(13.5, 19.5, facecolor='gold', alpha=0.5)
    plt.axhspan(19.5, 23.5, facecolor='plum', alpha=0.5)
    plt.axhspan(23.5, 28.5, facecolor='darkorange', alpha=0.5)
    plt.axhspan(28.5, 36.5, facecolor='green', alpha=0.5)

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("Node", fontsize=20)
    plt.xticks(fontsize=16, )
    plt.yticks(fontsize=16, )
    #plt.title("Regional Attractor Dynamics against $J_e$", fontsize=20)
    plt.legend(("Fixed Point","Limit Cycle"),loc="lower right",fontsize="x-large")
    #plt.savefig("do-not-track\\Het_Bifurcation_Flat.pdf",bbox_inches='tight')
    plt.show()

    return

def Show_Sig_sweep():
    # From Code in Progress.ipynb

    Regime = "LCycle"

    SCorr_files = glob.glob(r"do-not-track\2021_01_25\\" + Regime +"*Scorr*.csv")

    Params = []
    SCFC = []
    FCFC = []
    for string in SCorr_files:

        # Obtain Parameter Values
        #x = re.findall("\[(.*)\]sig_e(...).*sig_i(...)",string)
        #x = re.findall("\[(.*)\]sig_e(.).*sig_i(...)",string)
        x = re.findall("G\[(.*)\]_b_e.*sig_e(...).*sig_i(...)",string)
        Params.append(x[0])

        # FCFC and SCFC
        a = np.genfromtxt(string)
        SCFC.append(a[0])
        FCFC.append(a[2])

    df = pd.DataFrame(Params)
    df.columns = ['G','sig_e', 'sig_i']
    df["SCFC"] = SCFC
    df["FCFC"] = FCFC

    df_pivot = df.sort_values('FCFC').drop_duplicates(['sig_e','sig_i'],keep='last').pivot(index='sig_e', columns='sig_i', values='FCFC')
    df_pivot = df_pivot.sort_values('sig_e',ascending=True)
    # Reason why we need to do this is because the clsuter did something bad and simulated a few points multiple times (this is likely from it deciding to stop runs midway and then restarting)

    array = np.diag(df_pivot)
    plt.plot([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],array)
    plt.xlabel("$\sigma$",fontsize=20)
    plt.ylabel("FCFC",fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    return

def ShuffleValidations():
    # From ShuffleValidaitons.ipynb
    # First Import the runs:
    Best = glob.glob(r"do-not-track\2021_01_06\*Best*.csv")

    Score = []
    G = []
    Sigma = []

    for item in Best:
        a = np.genfromtxt(item)
        Score.append(a[0])
        G.append(a[1])
        Sigma.append(a[2])

    df = pd.DataFrame({"Score":Score,"G":G,"Sigma":Sigma})
    df = df.sort_values(by="Score")

    plt.hist(Score,cumulative=True,density=True,bins=11,)
    plt.hlines(y=0.95,xmin=0.555,xmax=0.667,color='b',linestyles='dashed')
    plt.hlines(y=0.85,xmin=0.555,xmax=0.667,color='r',linestyles='solid')
    #plt.title("CD Scores of Gradient Shuffled Simulations", fontsize=20)
    plt.ylabel("Cumulative Density", fontsize=20)
    plt.xlabel("Max FCFC Score", fontsize=20)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.legend(("95th Percentile","Percentile of Original Run (85th)"),loc=3,prop={'size': 15})
    #plt.grid()

    plt.show()

    return

def BenchmarkingVsLitereature():
# From Code in Progress.ipynb

    N = 5

    SCFC = (0.38, 0.30, 0.28, np.nan, 0.42)
    Homog = (0.46 , np.nan, 0.41, 0.68, 0.57)
    Het = (np.nan, 0.46, 0.56, 0.70, 0.61)
    error = (np.nan,np.nan,np.nan,np.nan,0.04),(np.nan,np.nan,np.nan,np.nan,0.06)
    ind = np.arange(N)    # the x locations for the groups
    width = 0.5       # the width of the bars: can also be len(x) sequence

    p3 = plt.bar(ind, Het, width, yerr=error)
    p2 = plt.bar(ind, Homog, width,)
    p1 = plt.bar(ind, SCFC, width)

    plt.ylabel('FCFC Score', fontsize=20)
    #plt.title('Benchmarking against Literature', fontsize=20)
    plt.xticks(ind, ('Melozzi, 2019', 'Wang, 2019', 'Demirtaş, 2019', 'Deco, 2020', 'Siu, 2020 ',), fontsize=14.5) # rotation=45, ha="right"
    plt.yticks(fontsize=16)
    plt.xlabel("Dynamical Models in Literature", fontsize=20)
    plt.legend((p1[0], p2[0],p3[0]), ('SCFC','Homogeneous','Heterogeneous'),loc=2, prop={'size': 15})
    plt.annotate("p = 0.15",xy=(4.05,0.65), fontsize=10)
    plt.show()

    return

def Hysteresis_Het_Bif_Diagrams(Regions,Scaling):
    # Originally from Bifurcation_Diagrams_Analysis.ipynb

    def sig(v):
    # Sigmoid function
        return 1/(1 + np.exp(-v))

    # Fixed pt - J_i=0,b_e=3,b_i=3.7,a_e=1,a_i=1,w_ee=12,w_ei=15,w_ie=10,w_ii=8,
    # LCycle - J_i=0,b_e=2,b_i=2.8,a_e=1,a_i=1,w_ee=11,w_ei=10,w_ie=10,w_ii=1,
    # LCycleH - J_i=0,b_e=3,b_i=4,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=15,w_ii=3,
    # Jump - J_i=0,b_e=5,b_i=3.7,a_e=1.3,a_i=2,w_ee=11,w_ei=10,w_ie=10,w_ii=1,


    #def func(z,J_e,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,tau_e,_tau_i):
    def func(z,J_e,J_i=0,b_e=5,b_i=3.7,a_e=1.3,a_i=2,w_ee=16,w_ei=12,w_ie=10,w_ii=3,):
        # Initialise the function array.
        f = np.zeros(2)

        E = z[0]
        I = z[1]

        # Wilson Cowan Equations setting derivative = 0 
        f[0] = -E + (1 - E)*sig(a_e*(w_ee*E - w_ei*I - b_e + J_e))
        f[1] = -I + (1 - I)*sig(a_i*(w_ie*E - w_ii*I - b_i + J_i))

        return f

    J_e = np.arange(0,5,0.05)
    array = []
    array_3 = []

    for i in J_e:
        
        array_2 = []

        for j in np.arange(100):
            # Brute force try random ones:
            sol = fsolve(func,[random.random()/2,random.random()/2],i,maxfev=100000,)
            # func(root) should be almost 0.0.
            if all(np.isclose(func(sol,i), [0.0, 0.0])):
                # Get only E value
                array_2.append(sol[0])
            else:
                array_2.append(np.nan)    

        # Then add the median, max, and min to the end of array

        array_3.append(array_2)


    Homog_Hyst = np.array(array_3)

    J_e = np.arange(0,5,0.05)
    array_final = []
    array_3 = []

    #for NUMBER in np.linspace(-0.1,0.1,num=Regions):
    for NUMBER in np.linspace(-Scaling,Scaling,num=Regions):

        # Calculate the parameters
        J_i=0
        b_e=5
        b_i=3.7
        a_e=1.3
        a_i=2
        w_ee=16*(1 + NUMBER)
        w_ei=12*(1 )
        w_ie=10*(1 + NUMBER)
        w_ii=3*(1 )

        array_3 = []

        for i in J_e:        
            
            array_2 = []

            # Get 100 different points
            for j in np.arange(100):
                # Brute force try random ones:
                sol = fsolve(func,[random.random()/2,random.random()/2],args=(i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,),maxfev=100000,)
                # func(root) should be almost 0.0.
                if all(np.isclose(func(sol,i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,), [0.0, 0.0])):
                    # Get only E value
                    array_2.append(sol[0])
                else:
                    array_2.append(np.nan)    
            # Add all 100 points which are attributed to the same J_e
            array_3.append(array_2)

        # Append each het bif onto array_final
        array_final.append(np.array(array_3))


    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    # Plot Homog case:
    plt.plot(J_e,Homog_Hyst,'k.',markersize=3)

    # Plot rest of Bif Diagrams:
    for I in np.arange(Regions):

        piece = array_final[I]
        # Plot 1 bifurcation diagram
        for i in np.arange(100):        
            plt.scatter(J_e[i]*np.ones(100),piece[i],color=colors[I],s=3)

    # Just edit in the legend in inkscape
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au) ", fontsize=20)
    plt.title("$R_e$ only", fontsize=20)
    plt.show()

    J_e = np.arange(0,5,0.05)
    array_final = []
    array_3 = []

    #for NUMBER in np.linspace(-0.1,0.1,num=Regions):
    for NUMBER in np.linspace(-Scaling,Scaling,num=Regions):

        # Calculate the parameters
        J_i=0
        b_e=5
        b_i=3.7
        a_e=1.3
        a_i=2
        w_ee=16*(1 )
        w_ei=12*(1 + NUMBER)
        w_ie=10*(1 )
        w_ii=3*(1 + NUMBER)

        array_3 = []

        for i in J_e:        
            
            array_2 = []

            # Get 100 different points
            for j in np.arange(100):
                # Brute force try random ones:
                sol = fsolve(func,[random.random()/2,random.random()/2],args=(i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,),maxfev=100000,)
                # func(root) should be almost 0.0.
                if all(np.isclose(func(sol,i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,), [0.0, 0.0])):
                    # Get only E value
                    array_2.append(sol[0])
                else:
                    array_2.append(np.nan)    
            # Add all 100 points which are attributed to the same J_e
            array_3.append(array_2)

        # Append each het bif onto array_final
        array_final.append(np.array(array_3))


    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    # Plot Homog case:
    plt.plot(J_e,Homog_Hyst,'k.',markersize=3)

    # Plot rest of Bif Diagrams:
    for I in np.arange(Regions):

        piece = array_final[I]
        # Plot 1 bifurcation diagram
        for i in np.arange(100):        
            plt.scatter(J_e[i]*np.ones(100),piece[i],color=colors[I],s=3)

    # Just edit in the legend in inkscape
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au) ", fontsize=20)
    plt.title("$R_i$ only", fontsize=20)
    plt.show()

    J_e = np.arange(0,5,0.05)
    array_final = []
    array_3 = []

    #for NUMBER in np.linspace(-0.1,0.1,num=Regions):
    for NUMBER in np.linspace(-Scaling,Scaling,num=Regions):

        # Calculate the parameters
        J_i=0
        b_e=5
        b_i=3.7
        a_e=1.3
        a_i=2
        w_ee=16*(1 + NUMBER)
        w_ei=12*(1 + NUMBER)
        w_ie=10*(1 + NUMBER)
        w_ii=3*(1 + NUMBER)

        array_3 = []

        for i in J_e:        
            
            array_2 = []

            # Get 100 different points
            for j in np.arange(100):
                # Brute force try random ones:
                sol = fsolve(func,[random.random()/2,random.random()/2],args=(i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,),maxfev=100000,)
                # func(root) should be almost 0.0.
                if all(np.isclose(func(sol,i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,), [0.0, 0.0])):
                    # Get only E value
                    array_2.append(sol[0])
                else:
                    array_2.append(np.nan)    
            # Add all 100 points which are attributed to the same J_e
            array_3.append(array_2)

        # Append each het bif onto array_final
        array_final.append(np.array(array_3))


    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    # Plot Homog case:
    plt.plot(J_e,Homog_Hyst,'k.',markersize=3)

    # Plot rest of Bif Diagrams:
    for I in np.arange(Regions):

        piece = array_final[I]
        # Plot 1 bifurcation diagram
        for i in np.arange(100):        
            plt.scatter(J_e[i]*np.ones(100),piece[i],color=colors[I],s=3)

    # Just edit in the legend in inkscape
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au) ", fontsize=20)
    plt.title("$R_e = R_i$", fontsize=20)
    plt.show()


    J_e = np.arange(0,5,0.05)
    array_final = []
    array_3 = []

    #for NUMBER in np.linspace(-0.1,0.1,num=Regions):
    for NUMBER in np.linspace(-Scaling,Scaling,num=Regions):

        # Calculate the parameters
        J_i=0
        b_e=5
        b_i=3.7
        a_e=1.3
        a_i=2
        w_ee=16*(1 + NUMBER)
        w_ei=12*(1 - NUMBER)
        w_ie=10*(1 + NUMBER)
        w_ii=3*(1 - NUMBER)

        array_3 = []

        for i in J_e:        
            
            array_2 = []

            # Get 100 different points
            for j in np.arange(100):
                # Brute force try random ones:
                sol = fsolve(func,[random.random()/2,random.random()/2],args=(i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,),maxfev=100000,)
                # func(root) should be almost 0.0.
                if all(np.isclose(func(sol,i,J_i,b_e,b_i,a_e,a_i,w_ee,w_ei,w_ie,w_ii,), [0.0, 0.0])):
                    # Get only E value
                    array_2.append(sol[0])
                else:
                    array_2.append(np.nan)    
            # Add all 100 points which are attributed to the same J_e
            array_3.append(array_2)

        # Append each het bif onto array_final
        array_final.append(np.array(array_3))


    jet= plt.get_cmap('coolwarm')
    colors = jet(np.linspace(0,1,Regions))

    # Plot Homog case:
    plt.plot(J_e,Homog_Hyst,'k.',markersize=3)

    # Plot rest of Bif Diagrams:
    for I in np.arange(Regions):

        piece = array_final[I]
        # Plot 1 bifurcation diagram
        for i in np.arange(100):        
            plt.scatter(J_e[i]*np.ones(100),piece[i],color=colors[I],s=3)

    # Just edit in the legend in inkscape
    #plt.legend(("Homogeneous","--50%","--30%","--10%","+10%","+30%","+50%"),loc="lower right",fontsize="x-large")

    plt.xlabel("$J_e$ (mV)", fontsize=20)
    plt.ylabel("$E$ (au) ", fontsize=20)
    plt.title("$R_e = -R_i$", fontsize=20)
    plt.show()

    return

def J_e_Border_Check():

    # Originally from Code in Progress.ipynb

    # For J_e Bif Checker

    # Import df: 

    df = pd.read_csv('do-not-track/J_e_LCycle_stats.csv')

    df_pivot_min = df.pivot(index='B_e', columns='G', values='TheMin')
    df_pivot_max = df.pivot(index='B_e', columns='G', values='TheMax')
    df_pivot_median = df.pivot(index='B_e', columns='G', values='TheMedian')
    df_pivot_mean = df.pivot(index='B_e', columns='G', values='TheMean')

    df = (df_pivot_max.sub(df_pivot_max.index,axis=0) < -1.45)
    # ( (df_pivot_max.sub(df_pivot_max.index,axis=0) > -1.75)) & 
    df_pivot = df_pivot_max
    x = df_pivot.index[::5]
    y = x.astype(np.float)
    X = df_pivot.columns[::10]
    Y = X.astype(np.float)

    cs=plt.imshow(df, cmap=ListedColormap(turbo_colormap_data), aspect='equal', interpolation='none',origin='lower',)
    #plt.title(Regime + ' Regime - FCFC', fontsize=20)
    axcb=plt.colorbar(cs)
    axcb.set_label('True/False', fontsize=20)
    axcb.ax.tick_params(labelsize=16)

    plt.yticks(ticks=np.arange(len(df_pivot.index))[::5],labels=y,fontsize=16)
    plt.xticks(ticks=np.arange(len(df_pivot.columns))[::10],labels=Y,fontsize=16)
    plt.xlabel("$G$", fontsize=20)
    plt.ylabel("$B_e$ (mV)", fontsize=20)
    plt.show()

    return