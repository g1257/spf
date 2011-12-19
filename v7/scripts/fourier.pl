#!/usr/bin/perl
# Inputs the correlation lattice
# and produces the S(q)=\sum e^(iqx) S_w . S_{w+x}
# 
use strict;

my $sqrt2= sqrt(2.0);

my ($filename,$GlobalLatticeName,$label) = @ARGV;
if (!defined($filename) or !defined($GlobalLatticeName)) {
	die "USAGE is: $0 filename latticeName Label\n";
}

print STDERR "#Arguments are $filename $GlobalLatticeName $label\n";
my $Dim = 2;

$label = "#ClassicalCorrelations:" if (!defined($label));

my @cc;
loadData(\@cc,$label);
my $L = $#cc+1;
print STDERR "#Read $L points\n";

$L /= 2 if ($GlobalLatticeName eq "square45degrees");
$L = sqrt($L);
print STDERR "#Computed lattice length is $L\n";

procData(\@cc,$#cc+1);

sub loadData {
	my ($cc,$finalLabel)=@_;
	
	open(FILE,$filename) or die "Cannot open file $filename$!\n";
	while(<FILE>) {
		last if (/^$finalLabel/);
	}
	$_=<FILE>; # size;
	chomp;
	my $n = $_;	
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
			calcComponentsReciprocal($m,\@mm);	 
			my @qq;
			calcComponentsDirect($x,\@qq);
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

sub calcComponentsReciprocal
{
	my ($q,$qq)=@_;
	if ($GlobalLatticeName eq "square") {
		calcComponentsSquare($q,$qq);
		return;
	}
	($GlobalLatticeName eq "square45degrees") or die "$0: Geometry $GlobalLatticeName is not supported yet (sorry)\n";
	calcCompSquare45DegreesReciprocal($q,$qq);
}

sub calcComponentsDirect
{
	my ($q,$qq)=@_;
	if ($GlobalLatticeName eq "square") {
		calcComponentsSquare($q,$qq);
		return;
	}
	($GlobalLatticeName eq "square45degrees") or die "$0: Geometry $GlobalLatticeName is not supported yet (sorry)\n";
	calcCompSquare45DegreesDirect($q,$qq);
}

# only one function needed for both direct and reciprocal
sub calcComponentsSquare
{
	my ($q,$qq)=@_;

	my $r=$q;
	for (my $i=$Dim-1;$i>0;$i--) {
		my $nn=$L**$i;
		$qq->[$i]=int($r/$nn);
		$r = $r % $nn;
	}
	$qq->[0]=$r;
}

# misses an additional sqrt(2), cancels with Reciprocal
sub calcCompSquare45DegreesDirect
{
	my ($q,$qq)=@_;
	
	my $xind =$q % (2*$L);
	my $yind = int($q / (2*$L));
	if ($xind>=$L) {
		$xind -= ($L - 0.5);
		$yind += 0.5;
	}
	$qq->[0] = $xind;
	$qq->[1]= $yind;
}

# has an additional sqrt(2), cancels with Direct
sub calcCompSquare45DegreesReciprocal
{
	my ($q,$qq)=@_;
	my ($Lprime1,$Lprime2)=($L/$sqrt2,$L*$sqrt2/3);
	if ($q<$L) {
		# first sub-lattice:
		$qq->[0] = $q % $Lprime1;
		$qq->[1] = int($q / $Lprime2);
		return;
	}

	# second sub-lattice;
	my $qMl = $q - $L;
	$qq->[0] = $qMl % $Lprime2;
	$qq->[1] = int($qMl / $Lprime1);
}

