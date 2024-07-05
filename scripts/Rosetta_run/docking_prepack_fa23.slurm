#!/bin/bash

#----------------------------------------------------
# SLURM job script to run RosettaDock MPI application 
# on Comet.
#
# usage: from within the directory where your flag files are,
# submit to cluster using sbatch 2SIC_docking.qsub
#----------------------------------------------------

#SBATCH --job-name=dock_ppk
#SBATCH --time=5:00:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output logs/docking.%j.out       # Name of stdout output file (%j expands to jobId), remember to make outerr directory
#SBATCH --error logs/docking.%j.err       # Name of stdout output file (%j expands to jobId), remember to make outerr directory

ml gcc
# module load git
# module load python/3.6
# module load openmpi/3.1
# module load gcc/5.5.0

# job description
ROSETTABIN=/home/rsamant2/scr16_jgray21/rsamant2/Rosetta/main/source/bin
ROSETTAEXE=docking_prepack_protocol
COMPILER=linuxgccrelease
EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER

# running with a date and time stamp
time $EXE @docking_prepack_fa23

