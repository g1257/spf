#!/usr/bin/perl
# Inputs the correlation lattice
# and produces the S(q)=\sum e^(iqx) S_w . S_{w+x}
# 
use strict;


my ($filename) = @ARGV;
my $Dim = 2;


my @cc;
loadData(\@cc);
my $L = sqrt($#cc+1);
procData(\@cc,$#cc+1);

sub loadData {
	my ($cc)=@_;
	my $finalLabel = "#ClassicalCorrelations:";
	
	open(FILE,$filename) or die "Cannot open file $filename$!\n";
	while(<FILE>) {
		last if (/^$finalLabel/);
	}
	#$_=<FILE>; # size;
	chomp;
	my $n = 64; #$_;	
	my $counter=0;
	while (<FILE>) {
		chomp;
		$cc->[$counter++]=$_;
		last if ($counter==$n);
	}
	close FILE;
}

sub procData {
	my ($cc,$n)=@_;
	my (@S,@Si);
	for (my $m=0;$m<$n;$m++) {
		$S[$m]=0;
		$Si[$m]=0;
		for (my $x=0;$x<$n;$x++) {
			my @mm;
			calcComponents($m,\@mm);	 
			my @qq;
			calcComponents($x,\@qq);
			my $sp=0;
			for (my $i=0;$i<$Dim;$i++) {
				$sp += $qq[$i]*$mm[$i];
			}
			$sp = $sp * 2 * 3.1415927;
			$sp = $sp / $L;
			$S[$m] += $cc->[$x]*cos($sp);
			$Si[$m] += $cc->[$x]*sin($sp);
		}
		print "$m $S[$m] $Si[$m]\n";
	}
}

sub calcComponents {
	my ($q,$qq)=@_;
				
	my $r=$q;
	for (my $i=$Dim-1;$i>0;$i--) {
		my $nn=$L**$i;
		$qq->[$i]=int($r/$nn);
		$r = $r % $nn;
	}
	$qq->[0]=$r;
}

	
