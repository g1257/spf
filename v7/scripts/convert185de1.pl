#!/usr/bin/perl -w
#
use strict;
use warnings;

while(<STDIN>) {
	if (/MonteCarloWindow/ && !/\[/) {
		chomp;
		my @temp=split;
		my $x = $temp[3];
		$_="MonteCarloWindow[SpinPhi]=$x\n";
	}
	print;
}

