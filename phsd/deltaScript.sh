#!/bin/bash

for i in {1425..1616}
do
    cd /Users/jtblair/4Justin_2018/carina_pPboutput10/$i
    eval aliroot -b -q -x \''~/utaustin/phsd/deltaReader.cxx("/Users/jtblair/4Justin_2018/carina_pPboutput10/'$i'/fort.433")'\'
done
