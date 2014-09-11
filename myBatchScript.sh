#!/bin/bash
#SBATCH -J 100m_minbias            # Job Name
#SBATCH -o 100m_minbias.o%j        # Output and error file name (%j expands to jobID)
#SBATCH -n 4000               # Total number of mpi tasks requested
#SBATCH -p normal           # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:30:00         # Run time (hh:mm:ss)i

ibrun minbiasParallel.o
