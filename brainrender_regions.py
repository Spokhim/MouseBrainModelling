""" 
    This tutorial shows how to create and render a brainrender scene with some brain regions
"""
import brainrender

brainrender.SHADER_STYLE = "cartoon"
from brainrender.scene import Scene

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
import scipy.cluster.hierarchy as hierarchy
#from tvb.simulator.lab import *
#from tvb.simulator.plot.tools import *
# Input Simulation Pipeline
#from SimulationPipeline import *
#from useful_fns import *

from brainrender.colors import *

# Add the whole thalamus in gray (think they meant red)
#scene.add_brain_regions(["TH"], alpha=0.15)
# Alpha is liek transparency. 
# Add VAL nucleus in wireframe style with the allen color
#scene.add_brain_regions(["VAL"], use_original_color=True, wireframe=True)

# Import the E and I densities
df = pd.read_csv(r"C:\Users\Pok Him\Desktop\MouseBrainModelling\CortexDensitiesAlter.csv",delimiter=",")
E_pop = df.excitatory.values
I_pop = df.inhibitory.values
E_mean = np.mean(E_pop)
I_mean = np.mean(I_pop)

# E_normalised is (when excluding region 7) -0.28 to 0.54
E_normalised = (E_pop-E_mean)/E_mean
# I_normalised is (when excluding region 7) -0.45 to 1.44
I_normalised = (I_pop-I_mean)/I_mean

"""
Regime = "LCycle"
G_value = 0.45
B_e_value = 2.8
File_start = r"D:\Simulations\2020_09_23\\"  + Regime + "_G[[]" + str(G_value) + "[]]_b_e[[]" + str(B_e_value) + "[]]"
TseriesFile = glob.glob(File_start+"*Tseries*_.csv")[0]

df = np.genfromtxt(TseriesFile,delimiter="\t")
bold_time = df[0]
bold_data = df[1:]

Region0 = bold_data[0]
print(Region0)
Value = Region0[0]
print(Value)
"""

# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
# Transparency
alpha = 0.15
print(len(E_normalised))

# For just normal colouring of brain regions. 
# Create a scene_e
scene_e = Scene(
    #screenshot_kwargs=dict(folder="Docs/Media/clean_screenshots"),
    title="Brain Regions",
)

for i in np.arange(len(E_normalised)):
    acronym = df.acronym[i]
    scene_e.add_brain_regions([acronym], alpha=0.5,use_original_color=True,wireframe=True)

scene_e.render()

'''
Peace=colorMap(E_normalised, name=ListedColormap(turbo_colormap_data),vmin=min(E_normalised),vmax=max(E_normalised))

# Create a scene_e
scene_e = Scene(
    #screenshot_kwargs=dict(folder="Docs/Media/clean_screenshots"),
    title="Excitatory Neuron Densities",
)

for i in np.arange(len(E_normalised)):
    colour_value = tuple(Peace[i])
    acronym = df.acronym[i]
    scene_e.add_brain_regions([acronym], alpha=alpha,color=colour_value)

scene_e.render()

# Create a scene_i
scene_i = Scene(
    #screenshot_kwargs=dict(folder="Docs/Media/clean_screenshots"),
    title="Inhibitory Neuron Densities",
)

# Transparency
alpha = 0.15
print(len(I_normalised))
Inhib = colorMap(I_normalised, name=ListedColormap(turbo_colormap_data),vmin=min(I_normalised),vmax=max(I_normalised))

for i in np.arange(len(I_normalised)):
    colour_value = tuple(Inhib[i])
    acronym = df.acronym[i]
    scene_i.add_brain_regions([acronym], alpha=alpha,color=colour_value)


scene_i.render()

# Create a scene
scene = Scene(
    #screenshot_kwargs=dict(folder="Docs/Media/clean_screenshots"),
    title="Brain",
)
scene.render()
'''