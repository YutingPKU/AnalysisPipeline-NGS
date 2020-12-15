#!/bin/bash
#Yuting Liu   Nov 26, 2016     First realease
#correct hmm probablity and find domain for each chromosome

dir=$1 #working directory
cd $dir
mkdir domain
min=2 #corrects probability if size of a cluster is <= min
prob=0.99 #checks probability in a cluster
binsize=4000 #change binsize according to your own data
#faifile="/lustre/user/liclab/liuyt/monkey-brain/ref_genome/rheMac8/rheMac8.chrom.size" #genome size file: you can choose hg18, hg19 and hg38.
#faifile="/lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/chrom_mm10.sizes"
faifile="/lustre/user/liclab/liuyt/monkey-brain/ref_genome/rheMac8/rheMac8.chrom.size"

for i in $(seq 21)
do
  colfile="$dir/7col/chr$i.txt"
  cd "$dir/domain"
  touch "chr$i.txt"
  domainfile="$dir/domain/chr$i.txt"
  cd "/lustre/user/liclab/liuyt/TADcalling/script/domaincall_software/perl_scripts"
  perl hmm_probablity_correcter.pl  $colfile $min $prob $binsize | perl hmm-state_caller.pl $faifile chr$i | perl hmm-state_domains.pl > $domainfile
done
