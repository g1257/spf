# generates the dynvars input file for the MCSTARTTYPE for spfv7
#
# usage: perl mcstart.pl "option" "size" > dynvarsfilename
#        where "option" is a string: either "ferro" or "antiferro"
#              "size" is the number of sites (integer)
#
# In spfv7 input file: MCSTARTTYPE dynvarsfilename
#
#!/usr/bin/perl -w

use strict;

my ($starttype,$size)=@ARGV; # here "size" must be a number whose square root is an integer
my $freezespins=0;
my $freezephonons=0;
my $PI=3.14159265359;
my $D=2;
my $L=sqrt($size);

if ($starttype eq "fm" or $starttype eq "ferro" or $starttype eq "ferromagnetic") {
	generatefileFM($size);
} elsif ($starttype eq "af" or $starttype eq "antiferro" or $starttype eq "antiferromagnetic") {
	generatefileAF($size);
} else {
	die "Cannot generate dynvars input file for unknown starttype: $starttype\n";
}

sub generatefileFM
{
	my ($linsize)=@_;
		
	print "Theta\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		print "0\n";
	}
	print "Phi\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		print "0\n";
	}
	print "IsFrozen $freezespins\n";
	print "Phonon\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		for (my $j=0; $j<3; $j++) {
			print "3\n";
			print "0\n";
			print "0\n";
			print "0\n";
			#print "\n";
		}
	}
	print "IsFrozenPhonon $freezephonons\n";
}

sub generatefileAF
{
	my ($linsize)=@_;
	my @theta;	
		
	print "Theta\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		if(parity($i,$D,$L)==1) {
			$theta[$i]=0;
		} else {
			$theta[$i]=$PI;
		}
		print "$theta[$i]\n";
	}
	print "Phi\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		print "0\n";
	}
	print "IsFrozen $freezespins\n";
	print "Phonon\n";
	print "$linsize\n";
	for (my $i=0; $i<$linsize; $i++) {
		for (my $j=0; $j<3; $j++) {
			print "3\n";
			print "0\n";
			print "0\n";
			print "0\n";
			#print "\n";
		}
	}
	print "IsFrozenPhonon $freezephonons\n";
}

sub parity
{
	my ($i,$d,$size)=@_;
	my ($x,$y,$z)=0;
	
	if ($d==1) {
		$x=$i;
		$y=0;
		$z=0;
	} elsif ($d==2) {
		$x=$i%$size;
		$y=($i-$x)/$size;
		$z=0;
	} elsif ($d==3) {
		$x=$i%$size;
		$z=($i-$x)/$size;
		$y=$z%$size;
		$z=($z-$y)/$size;
	} else {
		die "dim cannot be greater than 3 in parity\n";
	}
	$x += ($y+$z);
	if ($x%2==0) {
		return 1;
	} else {
		return 0;
	}
}
