#!/usr/bin/perl -w
use strict;

my @Betas = (100,80,70,60);
my $NumberOfBetas = $#Betas + 1;

for (my $i=0;$i<$NumberOfBetas;$i++) {
	writeOneInput($i);
}



sub writeOneInput
{
	my ($i)=@_;
	my $beta = $Betas[$i];
	my $previousRun ="none";
	if ($i>0) {
		# Ok this is run at least 1, so ...
		my $j = $i-1; # ... we start from the prev. one
		$previousRun="outputCoolDown$j"
	}
	open(FOUT,">inputCoolDown$i.inp") or die "Cannot open nputCoolDown$i.inp for writing: $!\n";

print FOUT<<EOF;
OPTIONS none
VERSION SPFv7
FILENAME outputCoolDown$i
CARRIERS  0
MU  0
BETA $beta
MCTHERMALIZATIONS 1000
MCEFFECTIVE 1000
MCUNMEASURED 1
MCWINDOW   1 0.5
MCSTARTTYPE $previousRun
MCSTARTLEVEL 0
HISTOGRAMSTEPS 1000
BOUNDARYCONDITIONS periodic
LINSIZE 8
hoppings
16
-0.058 0
0 -0.2196
-0.2196 0
0 -0.058
0.20828 +0.079
0.079 +0.20828
0.20828 -0.079
-0.079 +0.20828
PARAMETERSJ 0.0
POTENTIALV 64  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
JAFNN 0.
JAFNNN 0.1
#EOF
EOF

close(FOUT);
}

