#!/bin/bash
#SBATCH -J 10m_histo            # Job Name
#SBATCH -o 10m_histo.o%j        # Output and error file name (%j expands to jobID)
#SBATCH -n 400               # Total number of mpi tasks requested
#SBATCH -p normal           # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:30:00         # Run time (hh:mm:ss)i

ibrun histoParallel.o
