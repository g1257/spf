#!/usr/bin/perl -w

use strict;

my ($label)=@ARGV;

my ($sum,$counter)=(0,0);

# Read standard input, while's there's any
while(<STDIN>) {
	if (/^$label(.*$)/) {
		my $value = $1;
		print "$counter $value\n";
		$sum += $value;
		$counter++;
	}
}
$_ = $sum/$counter; # average
print "#$sum $counter $_\n";
	
