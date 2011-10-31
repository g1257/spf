#!/usr/bin/perl -w

use strict;

my @labels = ("#Theta","#Phi","#Phonons");

my $buffer="";
my $lblsCounter=0;

while(<STDIN>) {
	$buffer .= $_;
	my $label="";
	my $found=0;
	while(<STDIN>) {
		$buffer .= $_;
		my @temp = split(/\n/,$buffer);
		$buffer="";
		#print STDERR "$temp[0] $temp[1]\n" if ($#temp>=1);
		foreach (@temp) {
			$label=$_;
			if (foundLabel($label)) {
				$found=1;
				last;
			}
			print "$label\n";
		}
		last if ($found);
	}
	last if (!$found);
	print STDERR "Converting $label...\n";
	$lblsCounter++;
	$label=~s/(^#)//;
	my @x;
	my $counter=0;
	while(<STDIN>) {
		last if (/^#/);
		my @temp=split;
		$x[$counter++]=$temp[1];
	}
	$buffer = $_;

	printVector($label,\@x);
}

print "IsFrozen0\n";
print STDERR "Converted $lblsCounter labels\n";

sub printVector
{
	my ($label,$v)=@_;
	print "$label\n";
	my $n = scalar(@$v);
	print "$n\n";
	foreach (@$v) {
		print "$_\n";
	}
}

sub foundLabel
{
	my ($x)=@_;
	foreach my $l (@labels) {
		return 1 if ($x=~/^\Q$l/);
	}
	return 0;
}

