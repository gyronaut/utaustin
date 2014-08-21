#!/bin/bash
#SBATCH -J pTest            # Job Name
#SBATCH -o pTest.o%j        # Output and error file name (%j expands to jobID)
#SBATCH -n 16               # Total number of mpi tasks requested
#SBATCH -p development      # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:45:00         # Run time (hh:mm:ss)i

ibrun parallelTest.o
