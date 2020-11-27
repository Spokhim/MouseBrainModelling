""" 
    This tutorial shows how download and rendered afferent mesoscale projection data
    using the AllenBrainAtlas (ABA) and Scene classes
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

# Create a scene
scene = Scene(title="tractography")

#  Get data - using df just for the acronyms
df = pd.read_csv(r"C:\Users\Pok Him\Desktop\MouseBrainModelling\CortexDensitiesAlter.csv",delimiter=",")

name = df.acronym

for i in np.arange(len(name)):
    acronym = name[i]
    
    # Get the center of mass of the region of interest
    p0 = scene.atlas.get_region_CenterOfMass(acronym)

    # Get projections to that point
    tract = scene.atlas.get_projection_tracts_to_target(p0=p0)

    # Add the brain regions and the projections to it
    #scene.add_brain_regions([acronym], alpha=0.4, use_original_color=True)
    scene.add_tractography(tract, color_by="target_region")

scene.render()