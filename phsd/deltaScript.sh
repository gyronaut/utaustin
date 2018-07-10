#!/bin/bash

for i in {1..16}
do
    cd /Users/jtblair/4Justin_2018/200testoutput/$i
    eval aliroot -b -q -x \''~/utaustin/phsd/deltaReader.cxx("/Users/jtblair/4Justin_2018/200testoutput/'$i'/fort.433")'\'
done
