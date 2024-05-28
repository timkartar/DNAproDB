#!/bin/bash

##command to run ./production.sh 2>&1 | tee -a production_log.txt
while read p;
do
        if grep -q $p don.txt;then
                echo $p already done
        else
                echo "-------------------- RUNNING $p ----------------------"
                python auto_processStructure.py $p.cif 0
                echo $p >> don.txt
        fi
done <ids.txt
