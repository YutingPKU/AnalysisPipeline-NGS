#!/bin/bash
cd /lustre/user/liclab/liuyt/monkey-brain/07456AB/basicdata/matrix/07456A-CP/modify_TAD_boundary
indir=/lustre/user/liclab/liuyt/monkey-brain/07456AB/basicdata/matrix/07456A-CP/modify_insulation_mat
files=`find $indir -type f -name "*mat"`
for mat in $files
do
	CMD="perl /lustre/user/liclab/liuyt/software/crane-nature-2015-master/scripts/matrix2insulation.pl -i $mat -is 1000000 -ids 200000 -im mean  -nt 0.25 -v"
	echo $CMD
	eval $CMD
done
