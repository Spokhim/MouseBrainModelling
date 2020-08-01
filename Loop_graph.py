# Jupyter notebooks no longer working, so using python scritps instead.

# First Import the packages.
from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *
import numpy as np
import pylab
import matplotlib.pyplot as plt
#%matplotlib inline
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
from matplotlib.colors import ListedColormap
from turbo_colormap import *
from scipy import stats
#matplotlib.style.use('ggplot')
import inspect
import os
import csv 
import time

import glob
import pandas as pd

# Input Simulation Pipeline
from SimulationPipeline import *

# Now Import our data from our data folder:

# Get all csv filenames in a folder
all_files = glob.glob("do-not-track/*.csv")
#all_files

def par_extract(file):
    # Extracts the characters between square brackets.
    s_filter = ""
    x = False
    for i in file:
        if i == "[": x = True
        elif i == "]": x = False
        elif x: s_filter += i

    return s_filter

# Get Scorr csv filenames in a folder
Scorr_files = glob.glob("do-not-track/LCycle*Scorr*.csv")

SCFC = []
FCFC = []
G_value = []

# Loop to populate the empty array with the numbers from the Scorr csv files:

for item in Scorr_files:
    a = np.genfromtxt(item)
    SCFC.append(a[0])
    FCFC.append(a[2])
    G_value.append(par_extract(item))

print(Scorr_files)
print("SCFC =",SCFC)
print(numpy.amax(SCFC))
print(numpy.where(SCFC==np.amax(SCFC)))
print("FCFC =",FCFC)
print(numpy.amax(FCFC))
print(numpy.where(FCFC==np.amax(FCFC)))

# Check File order. 
# print(Scorr_files)

# Sort it:
def sorter(X,Y):
    Z = [x for _,x in sorted(zip(Y,X))]
    return Z

SCFC = sorter(SCFC,G_value)
FCFC = sorter(FCFC,G_value)
G_value = sorted(G_value)

# Graph it: 
plt.plot(G_value,SCFC)
plt.xlabel('G', fontsize=20)
plt.ylabel('SCorr - SCFC', fontsize=20)
plt.title('Scorr vs G (of FC vs SC)', fontsize=20)
plt.grid()
plt.show()

# Graph it: 
plt.plot(G_value,FCFC)
plt.xlabel('G', fontsize=20)
plt.ylabel('SCorr - FCFC', fontsize=20)
plt.title('Scorr vs G (of Sim vs Exp)', fontsize=20)
plt.grid()
plt.show()