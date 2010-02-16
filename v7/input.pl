#!/usr/bin/perl -w

use strict;

my ($choice)=@ARGV;
my $model="UNDEFINED";

if ($choice==0)  {
	$model = "PhononsTwoOrbitals";
	createInput0();
} else {
	$model = "PnictidesTwoOrbitals";
	createInput1();
}

sub createInput0
{
my $l = 8;
my $n = $l*$l;
my $V = 0.0;
my $potentialV  = "";
for (my $i=0;$i<$n;$i++) {
	$potentialV = $potentialV." $V ";
}
print <<EOF;
OPTIONS none
VERSION SPFv7
FILENAME output0
CARRIERS 48
MU -1
BETA 10
MCTHERMALIZATIONS 0
MCEFFECTIVE 10
MCUNMEASURED 1
MCWINDOW   2 0.5 0.3
MCSTARTTYPE none 
MCSTARTLEVEL 0
HISTOGRAMSTEPS 1000
BOUNDARYCONDITIONS periodic
LINSIZE $l
hoppings
8
1 -0.577350269  -0.577350269 0.33333333 1  0.577350269  0.577350269 0.33333333
POTENTIALV $n $potentialV
JAF 0.0
EJT 3  1.3 1.3 1.3
EDT 3  1 0.5 0.5
#EOF
EOF
}

sub createInput1
{
my $l = 8;
my $n = $l*$l;
my $V = 0.0;
my $potentialV  = "";
for (my $i=0;$i<$n;$i++) {
	$potentialV = $potentialV." $V ";
}
print <<EOF;
OPTIONS none
VERSION SPFv7
FILENAME output
CARRIERS 128
MU 0
BETA 100
MCTHERMALIZATIONS 0
MCEFFECTIVE 2
MCUNMEASURED 1
MCWINDOW   1 0.5
MCSTARTTYPE none 
MCSTARTLEVEL 0
HISTOGRAMSTEPS 1000
BOUNDARYCONDITIONS periodic
LINSIZE $l
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
POTENTIALV $n $potentialV
JAFNN 0.0
JAFNNN 0.0
#EOF
EOF
}

