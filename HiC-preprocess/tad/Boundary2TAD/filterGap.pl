#!/usr/bin/perl

use strict;
use warnings;

my %chr_gap;
open IN,"<","rheMac8.contig.hgTable.txt" or die "$!";
while(<IN>){
	# print "$_";
	chomp;
	my @in = split/\t/;
	my $chr = $in[1];
	my $start = $in[2];
	my $end = $in[3];
	my $start_end = "$start"."_$end";
	push @{$chr_gap{$chr}}, $start_end;
}
close IN;

open IN,"<","TAD.txt" or die "$!";
open OUT,">","TAD_filtered.txt" or die "$!";
while(<IN>){
	# print "$_";
	chomp;
	my @in = split/\t/;
	my $tad_chr = (split/\_/,$in[0])[0];
	# print "$tad_chr\n";
	my @gap = @{$chr_gap{$tad_chr}};
	my $count = 0;
	my $tad_start = $in[1];
	my $tad_end = $in[2];
	foreach my $i (0 .. $#gap){
		my ($gap_start, $gap_end) = split/\_/,$gap[$i];
		if($tad_start <= $gap_end && $gap_start <= $tad_end){
			$count++;
		}
	}
	if($count > 0){
		print "$tad_chr:$tad_start-$tad_end\n";
	}else{
		print OUT "$_\n";
	}
}
close IN;
close OUT;

