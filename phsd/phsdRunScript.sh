#!/bin/bash
phsdpath='/Users/jtblair/4Justin_carina/phsd-2018-06-08-3398578b5b0f'

for i in {21..30}
do
    echo 'Starting next run:'
    date
    cd /Users/jtblair/4Justin_carina/output/$i
    if [ $(($i % 3)) -eq 0 ]
    then
        nohup $phsdpath/phsd.r8
    else
        nohup $phsdpath/phsd.r8 &
    fi
    cd /Users/jtblair/4Justin_carina/output
done
