#!/bin/bash
#concatenate DI'files for each chromosome to make "whole_DI"
#chromosome X can be called as 23 henceforth
#Yuting Liu   Nov 26, 2016   First realease

dir=$1
cd "$dir/DI"
for i in $(seq 21)
do
  cat "$dir/DI/chr$i.txt" >> $dir/whole_DI.txt
done
