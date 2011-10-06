#!/usr/bin/perl -w

use strict;

my ($rootname, $L, $level, $dospins)=@ARGV;
my $filename="$rootname.sav";
my $D=2;
my $N=$L**$D;

my @theta;
my @phi;
my @phonons0;
my @phonons1;

my @GlobalL = ($L,$L);
my @GlobalLc = (0.5*$L,0.5*$L);
my $PI=3.1415927;

my $outputfile;

$outputfile=$rootname.".local.uq" if (!$dospins);
$outputfile=$rootname.".local.sq" if ($dospins);

#Read data
readVector($filename,$level,$N,"#Theta",\@theta) if ($dospins);
readVector($filename,$level,$N,"#Phi",\@phi) if ($dospins);
readVector($filename,$level,$N,"#Phonons0",\@phonons0) if (!$dospins);
readVector($filename,$level,$N,"#Phonons1",\@phonons1) if (!$dospins);

sub readVector
{
	my ($filename,$level,$N,$label,$array)=@_;
	my $i;
	open (FILE,$filename) or die("Cannot open $filename $!\n");
	
	for ($i=0;$i<$level;$i++) {
		while (<FILE>) {
			last if (/^$label/);
		}
	}
	while(<FILE>) {
		chop;
		if (/(^[^ ]+) +(.*$)/) {
			$array->[$1]=$2;
			last if ($1==$N-1);
		}
	}
	close(FILE);
}

#my $x=15;
for (my $x=0;$x<$N;$x++) {
	my @cd; # Cd(x)
	my @d;
	my $nOfDs = calcCdAndD($x,\@cd,\@d);
	#print STDERR  "nOfDs=$nOfDs $N\n";
	my @q;
	my @sq;
	if (!defined($outputfile)) {
		open(FILE,">$outputfile") or die "Cannot open file $outputfile for writing: $!\n";
 	} else {
		open(FILE,">>$outputfile") or die "Cannot open file $outputfile for appending: $!\n";
	}
	print FILE "#Cluster_$x\n";
	close(FILE);
	#print "at site = $x\n";
	my $nOfKs=calcSq(\@sq,\@q,\@cd,\@d,$nOfDs);
}

sub calcCdAndD
{
	my ($x,$cd,$d)=@_;
	my $counter=0;
	for (my $i=0;$i<$N;$i++) {
		for (my $j=0;$j<$N;$j++) {
		        #die "x=$x, i=$i, j=$j\n";
			next if (!isInCluster($x,$i));
			next if (!isInCluster($x,$j));
			
			my $thisD = calcDistance($i,$j);
			if (!isInVector($d,$thisD,$counter)) {
				$d->[$counter]=$thisD;
				$counter++;
			}
			if (!defined($cd->[$thisD])) {
				$cd->[$thisD] = calcCorrelation($i,$j);
			} else {
				$cd->[$thisD] += calcCorrelation($i,$j);
			}
			#print STDERR "\t\t\t $cd->[$thisD]\n";
			#$cd->[$thisD] -= ($phonons0[$j]**2 + $phonons1[$j]**2);
		}
	}	
	return $counter;			
}

sub isInCluster
{
	my ($x,$i)=@_;
	my @rx = index2Site($x);
	my @ri = index2Site($i);
	return isInCluster2(\@rx,\@ri);
}

sub isInCluster2
{
	my ($rx,$ri)=@_;
	return 0 if (!isInCluster3($rx,$ri,0));
	return 0 if (!isInCluster3($rx,$ri,1));
	return 1;
}
	
sub isInCluster3
{
	my ($rx,$ri,$xory)=@_;
	my $tmp = vecDistOneDim($ri,$rx,$xory); # wraps around
	return 0 if ($tmp>=$GlobalLc[$xory]);
	return 1;
}

sub index2Site
{
	my ($ind)=@_;
	my @r;
	$r[0]=$ind % $GlobalL[0];
	$r[1]=int($ind/$GlobalL[1]);
	return @r;
}

sub calcDistance
{
	my ($i,$j)=@_;
	my @ri = index2Site($i);
	my @rj = index2Site($j);
	my @dist;
	#print STDERR "i=$i j=$j distance=";
	# this is the vectorial distance in the x direction
	$dist[0] = vecDistOneDimB(\@ri,\@rj,0);
	# correct for boundary condition;
	#print STDERR " ($dist[0]) ";
	$dist[0] = $GlobalL[0] - $dist[0] if ($dist[0]>=$GlobalLc[0]);
	$dist[0] = -$GlobalL[0] -$dist[0] if ($dist[0]<= -$GlobalLc[0]);
	#print STDERR "$dist[0] ";
	
	# we add this number so that it is non-negative 
	$dist[0] += $GlobalLc[0] - 1;
	
	$dist[1] = vecDistOneDimB(\@ri,\@rj,1);
	#print STDERR " ($dist[1]) ";
	$dist[1] = $GlobalL[1] - $dist[1] if ($dist[1]>=$GlobalLc[1]);
	$dist[1] = -$GlobalL[1] -$dist[1] if ($dist[1]<= -$GlobalLc[1]);
	#print STDERR "$dist[1] ";
	
	# we add this number so that it is non-negative
	$dist[1] += $GlobalLc[1] - 1;
	# The max for $dist[0] is (2*$GlobalLc[0] - 2)
	# and so there are (2*$GlobalLc[0] - 1) of $dist[0]
	my $idx =  $dist[0] + $dist[1]*(2*$GlobalLc[0] - 1);
	#print STDERR " index=$idx\n";
	return $idx;
	
}

sub calcD
{
	my ($indexOfD)=@_;
	my $something = (2*$GlobalLc[0] - 1);
	my @d;
	$d[0] = $indexOfD % $something;
	$d[1] = int($indexOfD / $something);
	$d[0] -= $GlobalLc[0] - 1;
	$d[1] -= $GlobalLc[1] - 1;
	return @d;
}	

sub vecDistOneDim
{
	my ($ri,$rj,$whatDimension)=@_;
	my $tmp = $ri->[$whatDimension]-$rj->[$whatDimension];
	$tmp += $GlobalL[$whatDimension] if ($tmp<0);
	return $tmp;
}

sub vecDistOneDimB
{
	my ($ri,$rj,$whatDimension)=@_;
	my $tmp = $ri->[$whatDimension]-$rj->[$whatDimension];
	#$tmp += $GlobalL[$whatDimension] if ($tmp<0);
	return $tmp;
}

sub isInVector
{
	my ($v,$what,$nn)=@_;
	die "what not defined\n" if (!defined($what));
	for (my $i=0;$i<$nn;$i++) {
		die "$i $nn" if (!defined($v->[$i]));
		return 1 if ($v->[$i]==$what);
	}
	return 0;
}


sub calcCorrelation
{
	my ($i,$j)=@_;	
	return ($phonons0[$i]*$phonons0[$j] + $phonons1[$i]*$phonons1[$j]) - ($phonons0[$i]**2 + $phonons1[$i]**2) if(!$dospins);
	return cos($theta[$i])*cos($theta[$j]) + sin($theta[$i])*sin($theta[$j])*cos($phi[$i]-$phi[$j]);
}

sub calcSq
{
	my ($sq,$q,$cds,$d,$nOfDs)=@_;
	my $meshL = 16;
	#my $meshL = $GlobalLc[0];
	my $nOfKs = $meshL * $meshL; # $GlobalLc[0]*$GlobalLc[1];
	#my $nOfKs = $GlobalLc[0] * $GlobalLc[1];
	my @sqReal;
	my @sqImag;
	my $i=$nOfKs - 1;
	open(FILE,">>$outputfile") or die "Cannot open file $outputfile for appending: $!\n";
	for (my $i=0;$i<$nOfKs;$i++) { # loop over ks
		#print "at $i-th k out of $nOfKs k's\n" if (($i%10000==0));
		$sqReal[$i] = $sqImag[$i] = 0;
		my @tmp=calcKVector($i,$meshL);
		$q->[$i] = @tmp;
		for (my $j=0;$j<$nOfDs;$j++) { # loop over ds
			my @d = calcD($j);
			my $factor = 2 * $PI * scalarProduct(\@tmp,\@d);
			$factor /= $meshL;
			#die "Undef cds for $j\n" if (!defined($cds->[$j]));
			$sqReal[$i] += $cds->[$j]*cos($factor);
			$sqImag[$i] += $cds->[$j]*sin($factor);
			#print "\t $cds->[$j]*cos($factor) + \n";
		}
		#print "\n";
		$sqReal[$i] /= $GlobalLc[0]*$GlobalLc[1];
		$sqImag[$i] /= $GlobalLc[0]*$GlobalLc[1];
		#print STDERR "$i kx=$tmp[0] ky=$tmp[1] $sqReal[$i] $sqImag[$i]\n";
		print FILE "$i $sqReal[$i]\n";
	}
	close(FILE);
	return $nOfKs;
}	

sub calcKVector
{
	my ($kindex,$meshL)=@_;
	my @k;
	$k[0] = $kindex % $meshL;
	$k[1] = int($kindex / $meshL);
	return @k;
}

sub scalarProduct
{
	my ($q,$d)=@_;
	return $q->[0]*$d->[0]+$q->[1]*$d->[1];
	
}
