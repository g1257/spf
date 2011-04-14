#!/usr/bin/perl -w
#
use strict;
use Average;

my ($file) = @ARGV;


my $label = "CombinedMagnetizationSquared=";
my $mag2 = Average::average($label,$file);
print "<M.M>=$mag2\n";

my @m;
for (my $i=0;$i<3;$i++) {
	$label = "CombinedMagnetization$i=";
	$m[$i] = Average::average($label,$file);
}

$_ = scalarProduct(\@m,\@m);

print "<M>.<M>=".$_."\n";
print "<M.M>-<M>.<M> = ".($mag2-$_)."\n";

sub scalarProduct
{
	my ($v1,$v2)=@_;
	my $n = scalar(@$v1);
	my $sum = 0;
	for (my $i=0;$i<$n;$i++) {
		$sum += $v1->[$i] * $v2->[$i];
	}
	return $sum;
}

