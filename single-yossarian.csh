# This is the job scripts file to run a single job on Yossarian.

#!/bin/csh
#PBS -N single-yossarian
#PBS -q yossarian
#PBS -l select=1:ncpus=1:mem=4GB
# Minimum acceptable walltime: day-hours:mins:secs
#PBS -l walltime=48:00:00
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
#python Single_sim.py 
python J_e_Bif_Checker.py
exit