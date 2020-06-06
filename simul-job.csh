# This is the job scripts file to run jobs on the super cluster.

#!/bin/csh
#PBS -N simul-job
#PBS -l select=1:ncpus=12:mem=8GB
# Minimum acceptable walltime: day-hours:mins:secs
#PBS -l walltime=10:00:00
# Send e-mail if job ends or aborts
#PBS -m bea
#PBS -M psiu5120@uni.sydney.edu.au
#PBS -V

# Show the host on which the job ran
hostname
module load Anaconda3-5.1.0

# Load TVB virtual env.
source ~/TVB/bin/activate.csh

# Working directory is where I ran qsub
cd "$PBS_O_WORKDIR"
python Loop_sim.py
exit