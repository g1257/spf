#!/usr/bin/perl -w

package MyVector;
use strict;

sub get
{
	my ($v,$file,$label,$option)=@_;
	open(FILE,$file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		last if (/^$label/);
	}
	my $counter = 0;
	$_ = <FILE> if ($option);
	while(<FILE>) {
		chomp;
		last if (/^#/);
		$v->[$counter++]=$_;
	}
	close(FILE);
}

sub print
{
	my ($v)=@_;
	my $n = @$v;
	for (my $i=0;$i<$n;$i++) {
		print "$v->[$i]\n";
	}
}

sub sum
{
	my ($v1,$v2)=@_;
	my $n = @$v2;
	for (my $i=0;$i<$n;$i++) {
		if (defined($v1->[$i])) {
			$v1->[$i] += $v2->[$i];
		} else {
			$v1->[$i] = $v2->[$i];
		}
	}
}

1;