#!/usr/bin/perl

use strict;
use warnings;

my %chr_tad;
open IN,"<","boundary.txt" or die "$!";
while(<IN>){
	# print "$_";
	chomp;
	my ($chr, $start, $end) = split/\t/;
	push @{$chr_tad{$chr}}, $start;
}
close IN;

my %chr_end;
open IN,"<","rheMac8.chrom.sizes.brief" or die "$!";
while(<IN>){
	chomp;
	my @in = split/\t/;
	$chr_end{$in[0]} = $in[1];
}
close IN;

open OUT,">","TAD.txt" or die "$!";
foreach my $chr (sort keys %chr_tad){
	my $index = 0;
	my @starts = @{$chr_tad{$chr}};
	foreach my $i (0 .. $#starts-1) {
		$index++;
	  printf OUT "$chr"."_tad$index\t%.0f\t", $starts[$i];
	  printf OUT "%.0f\n", $starts[$i+1];
	}
}
close OUT;
