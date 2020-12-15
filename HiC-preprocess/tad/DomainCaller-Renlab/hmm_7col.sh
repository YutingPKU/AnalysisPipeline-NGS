#!/bin/bash
#convert hmm outfile into 7 col file
#Yuting Liu    Nov 26, 2016   First realease

dir=$1 #working directory
hmmi="$dir/whole_DI.txt" #hmm input file
hmmo="$dir/whole_hmm.txt" #hmm output file
cd $dir
touch "whole_7col.txt"
fileo="$dir/whole_7col.txt"
cd "/lustre/user/liclab/liuyt/TADcalling/script/domaincall_software/perl_scripts"
perl ./file_ends_cleaner.pl $hmmo $hmmi | perl ./converter_7col.pl > $fileo
