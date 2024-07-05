#!/bin/bash

#----------------------------------------------------
# SLURM job script to run RosettaDock MPI application 
# on Comet.
#
# usage: from within the directory where your flag files are,
# submit to cluster using sbatch 2SIC_docking.qsub
#----------------------------------------------------

#SBATCH --job-name=frun
#SBATCH --time=20:00:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --error logs/docking.%j.err       # Name of stdout output file (%j expands to jobId), remember to make outerr directory

ml gcc

# job description
ROSETTABIN=/home/aharmal1/scr16_jgray21/aharmal1/Rosetta/main/source/bin
ROSETTAEXE=mp_dock
COMPILER=linuxgccrelease
EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER
echo Starting MPI job running $EXE

#time $ROSETTABIN/docking_prepack_protocol.linuxgccrelease @docking_prepack_fa19

#echo "Done with Prepack"

# running with a date and time stamp
time mpirun -np 24 $EXE @fa19_docking_flags
