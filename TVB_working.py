from tvb.simulator.lab import *
LOG = get_logger('demo')
from tvb.simulator.plot.tools import *
import numpy as np
import pylab
import matplotlib.pyplot as plt

#%matplotlib inline
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

# Load the connectivity data from a zip file. 
con = connectivity.Connectivity.from_file(r'C:\Users\Pok Him\Desktop\paupau.zip')

# Visualize the structural connectivity matrix
plt.subplots()
cs=plt.imshow(np.log10(con.weights), cmap='jet', aspect='equal', interpolation='none')
plt.title('Structural connectivity matrix', fontsize=20)
axcb=plt.colorbar(cs)
axcb.set_label('Log10(weights)', fontsize=20)

plt.show()

# Simulation

# Set the parameter of the resting state simulation
sim = simulator.Simulator(model=models.ReducedWongWang(w=1.0, I_o=0.3),
                        connectivity=con,
                        coupling=coupling.Linear(a=0.096),
                        integrator=integrators.EulerStochastic(dt=0.1, noise=noise.Additive(nsig=0.000013)),
                        monitors=(monitors.Bold(period=2e3),
                                  monitors.TemporalAverage(period=1e3)),
                        simulation_length=1.2e5).configure()
# Run the resting state simulation
(bold_time, bold_data), _ = sim.run()

# BOLD timeseries Plots

# Display the simulated bold timeseries
plt.subplots()
plt.plot(bold_time,bold_data[:,0,:,0])
plt.xlabel('Time (ms)', fontsize=20)
plt.ylabel('Amplitude (au)', fontsize=20)
#plt.xlim(0,2e5)
plt.title('Simulated BOLD timeseries', fontsize=20)

plt.show()