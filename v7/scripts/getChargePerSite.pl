#!/usr/bin/perl -w

use strict;
use MyVector;

my $file = shift @ARGV;

my $label = "#LocalCharge";
my @v;
MyVector::get(\@v,$file,$label,1);

my $Norbitals = 2; # Number of orbitals
my $total = $#v + 1;
my $N = $total / (2*$Norbitals); # 2 is there for the spin

# indexing is:
# $v[$i + ($orbital + $spin*$Norbitals) * $N]

for (my $i = 0; $i< $N; $i++) {
	my $s = 0;
	for (my $x=0;$x<2*$Norbitals;$x++) {
		$s += $v[$i + $x*$N];
	}
	print "$i $s\n";
}
