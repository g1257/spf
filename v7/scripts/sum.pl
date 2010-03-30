#!/usr/bin/perl -w

use strict;

my ($label)=@ARGV;

while(<STDIN>) {
	last if (/^$label/);
}

my $counter = 0;
my $sum = 0;
$_= <STDIN>;
chomp;
my $total = $_;

while(<STDIN>) {
	chomp;
	last if (/^#/);
	$sum += $_;
	$counter++;
}

print "$sum $counter ".$sum/$counter."\n";
if ($total!=$counter) {
	print STDERR "WARNING: $total different from $counter\n";
}

