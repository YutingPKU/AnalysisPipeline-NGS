#! /bin/bash
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -o /data/brainSeq/brainSeq2/yingDATA/XSAnno_v2/data/byTranscript/H2P/scripts/liftOver_HP.output

# liftOver parameters
liftOver_minMatch="0.913"; # liftOver min match Sp1 to Sp2
Sp1ToSp2_chain="/home/liuyt/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step1_LiftOver/hg19ToRheMac8.over.chain";
Sp2ToSp1_chain="/home/liuyt/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step1_LiftOver/rheMac8ToHg19.over.chain";

# filter parameters
perl_dir="/home/liuyt/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/XSAnno/bin/perl_lib"
Sp1Anno="/home/liuyt/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step1_LiftOver/Ensembl.GRCh37.exons.new.bed";
output_dir="/home/liuyt/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step1_LiftOver/H2R";
Sp1="hg19";
Sp2="rheMac8";



BlatFilter $Sp1.2bit $Sp2.2bit $Sp1Anno $Sp2Anno $interID $interPL $intraID $interPL $perl_dir $out_dir $Sp1 $Sp2 $Sp1.11occ $Sp2.11occ
