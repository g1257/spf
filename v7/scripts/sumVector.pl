#!/usr/bin/perl -w
use strict;

my $file = shift @ARGV;

my @s;
my $total = 0;
foreach my $label (@ARGV) {
	my @v;
	getVector(\@v,$file,$label);
	sumVector(\@s,\@v);
	my $n1 = @s;
	my $n2 = @v;
	die "$0: $n1 != $n2\n" if ($n1 != $n2);
	$total++;
}
$s[0] /= $total;

my $finalLabel =$ARGV[1];
$finalLabel =~ s/.{1}$//;
print "$finalLabel"."Total\n";
printVector(\@s);

sub getVector
{
	my ($v,$file,$label)=@_;
	open(FILE,$file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		last if (/^$label/);
	}
	my $counter = 0;
	while(<FILE>) {
		chomp;
		last if (/^#/);
		$v->[$counter++]=$_;
	}
	close(FILE);
}

sub sumVector
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

sub printVector
{
	my ($v)=@_;
	my $n = @$v;
	for (my $i=0;$i<$n;$i++) {
		print "$v->[$i]\n";
	}
}

