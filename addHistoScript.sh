#!/bin/bash
#SBATCH -J 100m_histos            # Job Name
#SBATCH -o 100m_histos.o%j        # Output and error file name (%j expands to jobID)
#SBATCH -n 1                        # Total number of mpi tasks requested
#SBATCH -N 1
#SBATCH -p normal           # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:30:00         # Run time (hh:mm:ss)i

cd /scratch/03094/jtblair/alice/20140912_100m_batchrun_1/histograms 
hadd -k batch4_histos.root histo_*
