#!/usr/bin/perl -w
use strict;

my ($filename)=@ARGV;
my @hoppings;
my ($lx,$ly,$carriers);
my $nbands = 2;
my $Pi = 3.1415927;
readInput($filename,\@hoppings,\$lx,\$ly,\$carriers);


$carriers /= 2; # we do only one spin channel here
$carriers = $lx*$ly if ($carriers<=0);
print "$lx $ly $carriers $hoppings[0] $hoppings[15]\n";
my @seigs;
my @eigs;
calcEigenvals(\@eigs);
@seigs = sort {$a <=> $b} @eigs;
my $e = 2*fillEigs(\@seigs,$carriers);
print "Energy e = $e\n";

sub calcEigenvals
{
	my ($eigs)=@_;
	my $counter=0;
	
	for (my $mx=0;$mx<$lx;$mx++) {
		for (my $my=0;$my<$ly;$my++) {
			my @matrix;
			formMatrix(\@matrix,$mx,$my);
			my $delta = sqrt($matrix[0]*$matrix[0]+4*$matrix[1]*$matrix[1]-2*$matrix[0]*$matrix[3]+$matrix[3]*$matrix[3]);
			$eigs->[$counter++] = 0.5*($matrix[0]+$matrix[3]-$delta);
			$eigs->[$counter++] = 0.5*($matrix[0]+$matrix[3]+$delta);
		}
	}
	print STDERR "$counter eigenvalues\n";
	
			
}

sub fillEigs
{
	my ($seigs,$c) = @_;
	my $sum=0;
	for (my $i=0;$i<$c;$i++) {
		print STDERR "$i $seigs->[$i]\n";
		$sum += $seigs->[$i];
	}
	return $sum;
}

sub formMatrix
{
	my ($matrix,$mx,$my)=@_;
	my $kx = 2*$mx*$Pi/$lx;
	my $ky = 2*$my*$Pi/$ly;
	
	for (my $i=0;$i<4;$i++) {
		$matrix->[$i] = 2*($hoppings[$i+0]*cos($kx)+$hoppings[$i+1*$nbands*$nbands]*cos($ky)+
			$hoppings[$i+2*$nbands*$nbands]*cos($kx+$ky) + $hoppings[$i+3*$nbands*$nbands]*cos($kx-$ky));
	}
}

sub readInput
{
	my ($filename,$hoppings,$lx,$ly,$carriers)=@_;
	open(FILE,$filename) or die "Cannot open file $filename: $!\n";
	while(<FILE>) {
		chomp;
		if (/^carriers/i) {
			my @temp=split;
			$$carriers = $temp[1];
			next;
		}
		if (/^latticelengths/i) {
			my @temp=split;
			($$lx,$$ly)=split(/,/,$temp[1]);
			next;
		}
		if (/^hoppings/i) {
			readHoppings($hoppings);
			next;
		}
			
	}
	close(FILE);
}

sub readHoppings
{
	my ($hoppings)=@_;
	$_=<FILE>; # size
	my $tot = $nbands * $nbands * 4; # 4 directions x,y,x+y and x-y
	chomp;
	($_ == $tot) or die "Expecting 16 hoppings, got $_ instead\n"; 
	my $counter=0;
	while (<FILE>) {
		chomp;
		my @temp=split;
		for (my $i=0;$i<=$#temp;$i++) {
			$hoppings->[$counter++]=$temp[$i];
		}
		last if ($counter==$tot);
	}
}


