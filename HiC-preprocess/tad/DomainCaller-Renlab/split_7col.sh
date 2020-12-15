#!bin/bash
#Yuting Liu    Nov 26, 2016    First realease
#split whole genome 7col file based on chromosome

dir=$1 #working direcotry
cd $dir
whole_7col="$dir/whole_7col.txt"
mkdir 7col
cd $dir/7col
for i in $(seq 21)
do
  touch "$dir/7col/chr$i.txt"
  fileo="$dir/7col/chr$i.txt"
  grep "chr$i\s" $whole_7col  > $fileo
done
