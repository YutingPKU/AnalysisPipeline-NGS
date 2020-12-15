#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Basename;

#my $dir = getcwd;
#mkdir "LoopBatch";
my $fat4 = "gpu_4l";
my $f4 = "g4c";
my $cores = 4;
# my $fat4 = "fat4way";
# my $f4 = "f4w";
open OUT2,">","LoopBatch/LoopBatch_WTGM.sh" or die $!;
print OUT2 "#!/bin/bash\n\n";
open IN,"<","hic.ls" or die $!;
while(<IN>){
	chomp;
	# my @in = split/\,/;
	# my $id = $in[0];
	my $file = $_;
	my $id = basename($file);
	my $out = "$id"."_10kb_loop";
	my $loop = "/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/results/GZ-10k/merged_loops.bedpe";
	# if($fat4 eq "fat4long"){
	# 	$fat4 = "fat4way";
	# 	$f4 = "f4w";
	# }else{
	# 	$fat4 = "fat4long";
	# 	$f4 = "f4l";
	# }
	open OUT,">","LoopBatch/$id.sh" or die "$!";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH -J $id\n";
	print OUT "#SBATCH -o $id.%j.out\n";
	print OUT "#SBATCH -e $id.%j.err\n";
	print OUT "#SBATCH -p $fat4\n";
	print OUT "#SBATCH -N 1\n";
#	print OUT "#SBATCH -c $cores\n";
	print OUT "#SBATCH --gres=gpu:1\n";
	print OUT "#SBATCH --overcommit\n";
	print OUT "#SBATCH --mincpus=7\n";
	print OUT "#SBATCH --no-requeue\n";
	print OUT "#SBATCH -A lch3000_g1\n";
	print OUT "#SBATCH --qos=lch3000$f4\n\n";
#	print OUT "source /appsnew/source/gcc-8.3.0.sh\n\n";
#	print OUT "cd /lustre1/lch3000_pkuhpc/liuyt/ChenQ/Hi-C/results_FL0/$id;\n";
#	print OUT "bash hicpro2juicebox_cw.sh -i $file -g mm10.chrom.sizes -j /home/lch3000_pkuhpc/lustre1/liuyt/software/bin/juicer_tools.1.9.8_jcuda.0.8.jar -t . -o /lustre1/lch3000_pkuhpc/liuyt/ChenQ/Hi-C/hicfiles";
#	print OUT "time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools_1.13.02.jar  hiccups -r 10000 -k KR -f 0.1 -p 2  --ignore_sparsity  $file $out";
	print OUT "java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools_1.13.02.jar  hiccups -r 10000 -k KR -f 0.1 -p 2 -i 5  -d 20000 --ignore_sparsity $file $out $loop --threads 4 ";

	close OUT;

	print OUT2 "pkubatch $id.sh\n";
}
close IN;
close OUT2;

