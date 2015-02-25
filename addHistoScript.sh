#!/bin/bash
#SBATCH -J 100m_histos            # Job Name
#SBATCH -o 100m_histos.o%j        # Output and error file name (%j expands to jobID)
#SBATCH -n 1                        # Total number of mpi tasks requested
#SBATCH -N 1
#SBATCH -p normal           # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:30:00         # Run time (hh:mm:ss)i

cd /scratch/03094/jtblair/alice/20140911_100m_batchrun_1/invMass_histograms 
hadd -k invbatch1_histos.root invMass_histo_*
