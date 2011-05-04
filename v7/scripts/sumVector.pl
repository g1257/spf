#!/usr/bin/perl -w
use strict;
use MyVector;

my $file = shift @ARGV;

my @s;
my $total = 0;
foreach my $label (@ARGV) {
	my @v;
	MyVector::get(\@v,$file,$label);
	MyVector::sum(\@s,\@v);
	my $n1 = @s;
	my $n2 = @v;
	die "$0: $n1 != $n2\n" if ($n1 != $n2);
	$total++;
}
$s[0] /= $total;

my $finalLabel =$ARGV[1];
$finalLabel =~ s/.{1}$//;
print "$finalLabel"."Total\n";
MyVector::print(\@s);



