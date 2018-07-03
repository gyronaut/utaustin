#!/bin/bash
phsdpath='/Users/jtblair/4Justin_carina/phsd-2018-06-08-3398578b5b0f'

for i in {849..1040}
do
    echo 'Starting next run:'
    date
    cd /Users/jtblair/4Justin_carina/output7/$i
    if [ $(($i % 16)) -eq 0 ]
    then
        sleep 120s
        nohup $phsdpath/phsd.r8 &
        wait $!
    else
        nohup $phsdpath/phsd.r8 &
    fi
done
