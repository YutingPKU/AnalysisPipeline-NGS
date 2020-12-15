#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Basename;

my $dir = "/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_190119/Homer"; 
#mkdir "makeTag";
mkdir "findPeaksBatch";
open OUT2,">","findPeaksBatch/findPeaksBatch.sh" or die $!;
print OUT2 "#!/bin/bash\n\n";
open IN,"<","listfiles/tagDir.list" or die $!;
my $fat4 = "cn-short";
my $f4 = "cns";
# my $fat4 = "fat4way";
# my $f4 = "f4w";
my $cores = 20;
while(<IN>){
	chomp;
	# my @in = split/\,/;
	# my $id = $in[0];
	my $path = $_;
	my ($id, $dirs, $suffix) = fileparse($path, qr/\.[^.]*/);
	my $in = $path;
	#my $out = "$dir/find/$id";
	# if($fat4 eq "fat4long"){
	# 	$fat4 = "fat4way";
	# 	$f4 = "f4w";
	# }else{
	# 	$fat4 = "fat4long";
	# 	$f4 = "f4l";
	# }
	open OUT,">","findPeaksBatch/$id.fragL147.fdr0.01.size147.sh" or die "$!";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH -J $id\n";
	print OUT "#SBATCH -o $id.%j.out\n";
	print OUT "#SBATCH -e $id.%j.err\n";
	print OUT "#SBATCH -p $fat4\n";
	print OUT "#SBATCH -N 1\n";
	print OUT "#SBATCH --ntasks-per-node=$cores\n";
	print OUT "#SBATCH --no-requeue\n";
	print OUT "#SBATCH -A lch3000_g1\n";
	print OUT "#SBATCH --qos=lch3000$f4\n\n";
	print OUT "findPeaks  $in   \\\n";
	print OUT "          -o $in/fragL147.fdr0.01.size147.peaks.txt -size 147 -gsize 2.7e09 -norm 1000000 -fragLength 147 -fdr 0.01 \\\n";
	print OUT "          -style factor \\\n";
	close OUT;

	print OUT2 "pkubatch $id.fragL147.fdr0.01.size147.sh\n";
}
close IN;
close OUT2;

