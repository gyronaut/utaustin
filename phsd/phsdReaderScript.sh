#!/bin/bash

for i in {1..80}
do
    echo 'Starting next run:'
    cd /Users/jtblair/4Justin_2018/carina_pPboutput2/$i
    eval aliroot -b -q -x \''~/utaustin/phsd/phsdReader.cxx("/Users/jtblair/4Justin_2018/carina_pPboutput2/'$i'/inputPHSD", "/Users/jtblair/4Justin_2018/carina_pPboutput2/'$i'/phsd.dat")'\'
done
