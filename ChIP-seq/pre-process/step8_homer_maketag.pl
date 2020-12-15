#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Basename;

my $dir = "/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_190119/Homer"; 
#mkdir "makeTag";
mkdir "makeTagBatch";
open OUT2,">","makeTagBatch/makeTagBatch.sh" or die $!;
print OUT2 "#!/bin/bash\n\n";
open IN,"<","listfiles/merge.bam.list" or die $!;
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
	my $out = "$dir/makeTag/$id";
	# if($fat4 eq "fat4long"){
	# 	$fat4 = "fat4way";
	# 	$f4 = "f4w";
	# }else{
	# 	$fat4 = "fat4long";
	# 	$f4 = "f4l";
	# }
	open OUT,">","makeTagBatch/$id.sh" or die "$!";
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
	print OUT "makeTagDirectory  $out  \\\n";
	print OUT "          -mapq 30  \\\n";
	print OUT "          $in \\\n";
	close OUT;

	print OUT2 "pkubatch $id.sh\n";
}
close IN;
close OUT2;

