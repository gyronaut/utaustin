TACC Submission/Compilation Instructions:
=========================================

*For the current Event generating code, the following procedure is required to launch the jobs:*

+ source the alice environment script (/work/03093/deepat/alice/alice-env.sh) - alias shortcut= "ali"

+ change the gen_events_TACC.C code (only need to change the output folder)

+ compile the code (mpic++ gen_events_TACC.C -o minbiasParallel.o)

+ change myBatchScript.sh to call correct executable if necessary (also change walltime, number of cores, etc.)

+ launch job - sbatch myBatchScript


Next Steps for analysis:
-------------------------

+ Look at Hadron and Phi-meson phi distribution (make sure it's flat)

+ Try to reconstruct invariant masses from K+ K- pairs (look at mass distrubition for different delta-phi bins)

+ See if the new "nixing daughter particle correlations" actually fixed the delta-phi distribution
