#!/bin/bash

DIR=$1

for i in `ls "$DIR"/*mod`
do 
    j=${i/.mod/}
    j='model'`echo $j | awk -F"model" '{print $2}'`
    ~/Dropbox/Research/SCC/bin/drawpfsa $i 2 $j
    echo $j
done
