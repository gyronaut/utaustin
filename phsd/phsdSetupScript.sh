#!/bin/bash

for i in {36..40}
do
    mkdir $i
    cd $i
    echo '207,      MASSTA: target mass' >> inputPHSD
    echo '82,       MSTAPR: protons in target' >> inputPHSD
    echo '1,        MASSPR: projectile mass' >> inputPHSD
    echo '1,        MSPRPR: protons in projectile' >> inputPHSD
    echo '13433049.,    ELAB:  (=4060000. Lab energy per nucleon LHC),=21300 RHIC,=13433049 (5 TeV) = 26120000 (7 TeV)' >> inputPHSD
    echo '2.0,      BMIN:   minimal impact parameter in fm ! no effect for p+A' >> inputPHSD
    echo '2.0,       BMAX:   maximal impact parameter in fm ! no effect for p+A' >> inputPHSD
    echo '0.5,      DeltaB: impact parameter step in fm (used only if IBweight_MC=0)' >> inputPHSD
    echo '30,       NUM:    optimized number of parallel ensambles ("events")' >> inputPHSD
    echo '5,        ISUBS:  number of subsequent runs' >> inputPHSD
    echo $((2*i - 1))',    ISEED:  ANY uneven INTEGER number' >> inputPHSD
    echo '1,        IGLUE: =1 with partonic QGP phase (PHSD mode); =0 - HSD mode' >> inputPHSD 
    echo '140.,      FINALT: final time of calculation in fm/c' >> inputPHSD
    echo '10,       ILOW: output level (default=10)' >> inputPHSD
    echo '0,        Idilept: =0 no dileptons; =1 electron pair; =2  muon pair' >> inputPHSD
    echo "0,        ICQ: =0 free rho's, =1 dropping mass, =2 broadening, =3 drop.+broad." >> inputPHSD
    echo '0,        IHARD: =1 with charm and bottom; =0 - without' >> inputPHSD
    echo '0,        IBweight_MC: =0 constant step in B =DBIMP; =1 choose B by Monte-Carlo ISUBS times in [Bmin,Bmax]' >> inputPHSD
    echo '0,        IUSER: =1 for general users : use default /optimized settings; = 0 for PHSD team' >> inputPHSD
    cd ..
done
