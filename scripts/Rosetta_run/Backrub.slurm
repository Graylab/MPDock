#!/bin/bash

#----------------------------------------------------
# SLURM job script to run RosettaDock MPI application 
# on Comet.
#
# usage: from within the directory where your flag files are,
# submit to cluster using sbatch 2SIC_docking.qsub
#----------------------------------------------------

#SBATCH --job-name=br_ligandA
#SBATCH --time=5:00:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --output logs/br_structure.%j.out       # Name of stdout output file (%j expands to jobId), remember to make outerr directory
#SBATCH --error logs/br_structure.%j.err       # Name of stdout output file (%j expands to jobId), remember to make outerr directory

module load compilers/gcc/8.5.0

# job description
ROSETTABIN=/home/r/rituparnasamanta/Software/rosetta_jan23/source/bin
ROSETTAEXE=backrub
COMPILER=linuxgccrelease
EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER

# running with a date and time stamp
time $EXE @backrub_flags

