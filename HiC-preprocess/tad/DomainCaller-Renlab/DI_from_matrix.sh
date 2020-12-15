#!/bin/bash
#calculate directory index(DI) of every chromosome annoted contact matrix
#Yuting Liu   Nov 26, 2016  First Realease

dir=$1 #working directory
mkdir $dir/DI
binsize=40000 #change bin size according to your own data
windowsize=2000000 #window size means when calcualte DI from up/down stream of 2M
#genomefile="/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai" #genome size file: you can choose hg18, hg19 and hg38.
genomefile="/lustre/user/liclab/liuyt/monkey-brain/ref_genome/rheMac8/rheMac8.chrom.size"
#genomefile="/lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/chrom_mm10.sizes"
for i in $(seq 21)
do
	filei="/home/liuyt/lustrelyt/monkey-brain/merged-hic/basic-data/matrix/merged-CP/iced/40000/modify_CP_40000_iced_chr${i}_dense.matrix"
	cd "$dir/DI"
	touch "chr${i}.txt"
	fileo="$dir/DI/chr${i}.txt"
	CMD="perl -w "/lustre/user/liclab/liuyt/TADcalling/script/domaincall_software/perl_scripts/DI_from_matrix.pl" $filei $binsize $windowsize $genomefile > $fileo"
	echo $CMD
	eval $CMD
done
#in file of chromosome X the first column is "X"
#we need to replace it with "23" for next process
#filex="$dir/DI/chr23.txt"
#Rscript /lustre/user/liclab/liuyt/TADcalling/script/shellscript/convert_x_23.r $filex
# awk '{$1=23;print ;}' chr23.txt > chr23.tmp
# tr ' ' '\t' < chr23.tmp > chr23.txt
