#!/usr/bin/perl -w

print "************************************\n";
print "* Welcome to input.pl for SPF v6.4\n";
print "************************************\n\n";


print "What type of model would you like to compile? (spell correctly!!)\n";
print "Available: MODEL_KONDO_INF_ONEBAND, MODEL_KONDO_INF_ONEBAND_PHONONS, MODEL_KONDO_INF_TWOBANDS MODEL_KONDO_FINITE ";
print " MODEL_KONDO_DMS_ZINCBLENDE or MODEL_KONDO_DMS_CUBIC or MODEL_KONDO_DMS_THREEBANDS or MODEL_KONDO_BCS\n";
print "Default (press ENTER): MODEL_KONDO_FINITE\n";
$ret=<STDIN>;
$ret = "MODEL_KONDO_FINITE\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$model=$ret;

print "Enter the options you want to use.\n";
print "Available: NOT LISTED HERE YET. SEE THE DOCS.\n";
print "Default (press ENTER): none\n";
$ret=<STDIN>;
$ret = "none\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$options=$ret;


print "Enter the lattice type.\n";
print "Available: 1d square cubic fcc bcc triangular honeycomb prism rectangular \n";
print "Default (press ENTER): square\n";
$ret=<STDIN>;
$ret = "square\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$latticetype=$ret;

if ($latticetype=~/honeycomb/ || $latticetype=~/rectangular/) { 
print "Enter the two lengths for this lattice separated by a comma (don't enter spaces).\n";
print "Default (press ENTER): 8,8\n";
$ret=<STDIN>;
$ret = "8,8\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$latticelengths=$ret;
} elsif ($latticetype=~/prism/) {
	print "Enter the three lengths for this lattice separated by a comma (don't enter spaces).\n";
print "Default (press ENTER): 4,4,4\n";
$ret=<STDIN>;
$ret = "4,4,4\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$latticelengths=$ret;
} else {
print "Enter the length for this lattice\n";
print "Default (press ENTER): 8\n";
$ret=<STDIN>;
$ret = "8\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
if ($latticetype=~/1d/) {
	$latticelengths=$ret;
	$dimension=1;
} elsif ($latticetype=~/square/ || $latticetype=~/triangular/) {
	$latticelengths=$ret.",".$ret;
	$dimension=2;
} elsif ($latticetype=~/cubic/ || $latticetype=~/fcc/ || $latticetype=~/bcc/) {
	$dimension=3;
	$latticelengths="$ret,$ret,$ret";
} else {
	die "input.pl: Invalid lattice type $latticetype\n";
}
} # latticelengths
if ($model=~/MODEL_KONDO_BCS/) {
	$carriers= -1;
	$mu=0;
} else {
	print "Do you want to automatically adjust the chemical potential for this simulation?\n";
	print "Available: Yes or No\n";
	print "Default (press ENTER): Yes\n";
	$ret=<STDIN>;
	$ret = "Yes\n" if ($ret eq "" or $ret eq "\n");
	chop($ret);
	if ($ret=~/yes/i) {
		print "Enter the target value for number of electrons.\n";
		print "Available: Any positive number.\n";
		$ret=<STDIN>;
		$ret = "5\n" if ($ret eq "" or $ret eq "\n");
		chomp($ret);
		$carriers=$ret;
		$ret="initial";
	} else {
		$carriers= -1;
		$ret="";
	}
	print "Enter the $ret value of the chemical potential.\n";
	print "Available: Any.\n";
	$ret=<STDIN>;
	$ret = "0\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$mu=$ret;
}
$beta_label="BETA";
if ($options=~/temperature/) {
	$beta_label="TEMPERATURE";
}
print "Enter the value of the $beta_label.\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "100\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$beta=$ret;

print "Enter the rootname for output files.\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "output\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$rootname=$ret;

print "Enter the number of Monte Carlo thermalization steps.\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "0\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$mctherm=$ret;

print "Enter the number of Monte Carlo measurements steps.\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "1\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$mceffective=$ret;

print "Enter 1 plus the number of steps to be left unmeasured between measurement steps (does not affect thermalization).\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "1\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$mcunmeasured=$ret;

print "Enter Monte Carlo window for the update of classical fields.\n";
print "Suggested (press ENTER): 0.5\n";
print "Available: Any.\n";
$ret=<STDIN>;
$ret = "0.5\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$mcwindow=$ret;

print "Do you want that all Monte Carlo steps be rejected (frozen calculation)?\n";
print "Available: Yes or No\n";
print "Default (press ENTER): No\n";
$ret=<STDIN>;
$ret = "No\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$mcflag=0;
$mcflag=0 if ($ret=~/yes/i);

print "What should be the initial value for the classical fields?\n";
print "Available: Ferromagnetic Random Antiferromagnetic Planar FromFile CE\n";
print "Default (press ENTER): Random\n";
$ret=<STDIN>;
$ret = "Random\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
if ($ret=~/Ferromagnetic/i) {
	$mcstarttype=0;
} elsif ($ret=~/Random/i) {
	$mcstarttype=1;
} elsif ($ret=~/Antiferromagnetic/i) {
	$mcstarttype=2;
} elsif ($ret=~/planar/i) {
	$mcstarttype=3;
} elsif ($ret=~/CE/) {
	$mcstarttype=5;
} elsif ($ret=~/FromFile/) {
	print "Enter a file in sav format from which to read the classical fields.\n";
	print "Available: any\n";
	$ret=<STDIN>;
	chomp($ret);
	$mcstartfile=$ret;
	
	print "If this file has many sets of data, please enter the set number you want to use, otherwise press ENTER.\n";
	$ret=<STDIN>;
	$ret = "0\n" if ($ret eq "" or $ret eq "\n");
	chomp;
	if ($ret==0) {
		$mcstarttype=4;
	} else {
		$mcstartype=6;
		$mcstartlevel=$ret;
	}
}

print "In how many subintervals should the energy interval be divided into to calculate Histograms?\n";
print "Suggested (press ENTER): 1000\n";
$ret=<STDIN>;
$ret = "1000\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$histsteps=$ret;

if ($options=~/histogrambounds/) {
print "Enter the energy bounds for the histograms separated by a comma but do not enter spaces or tabs.\n";
print "Example (press ENTER): -4,5\n";
$ret=<STDIN>;
$ret = "-4,5\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$histbounds=$ret;
}

$n=nofspins();

if ($model=~/MODEL_KONDO_FINITE/) {
print "Enter the number of spins.\n";
print "Default (press ENTER): $n\n";
 $ret=<STDIN>;
$ret = "$n\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$conc=$ret;
} else {
	$conc = $n;
}


#print "Enter the boundary conditions as $dimension words from the set {pbc,abc,open} separated by spaces.\n";
#$n="pbc";
#for ($i=1;$i<$dimension;$i++) { $n="$n pbc"; }
print "Enter the boundary conditions (periodic, open, antiperiodic or ".(2*$D)." numbers.\n";
print "Default (press ENTER): periodic\n";
$ret=<STDIN>;
$ret = "periodic\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
$bc=$ret;


if ($model=~/MODEL_KONDO_FINITE/) {
	print "Enter the value of J or JH.\n";
	print "Default (press ENTER): 8\n";
	$ret = <STDIN>;
	$ret = "8\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$J=$ret;
} elsif ($model=~/MODEL_KONDO_DMS/) {
	print "Enter the two values for J or JH separated by a space.\n";
	print "Default (press ENTER): 8 8\n";
	$ret = <STDIN>;
	$ret = "8 8\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$J=$ret;
}

$n = nofspins();
$Jaf=0;
if ($options=~/jafvector/) {
	print "Enter the file that contains the ".(2*$n)." values for JAFVECTOR.\n";
	print "Default (press ENTER): jafvector.txt\n";
	$ret = <STDIN>;
	$ret = "jafvector.txt\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$Jaf=$ret;
} else {
	if ($model=~/MODEL_KONDO_BCS/) {
		$Jaf=0;
	} elsif ($n==$conc) {
		print "Enter the value for JAF.\n";
		print "Default (press ENTER): 0\n";
		$ret = <STDIN>;
		$ret = "0\n" if ($ret eq "" or $ret eq "\n");
		chomp($ret);
		$Jaf=$ret;
	}
	
}

if ($options=~/havepotential/) {
	print "Enter the file that contains the values of the potential.\n";
	print "Default (press ENTER): potential.txt\n";
	$ret = <STDIN>;
	$ret = "potential.txt\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$potential=$ret;
}

if ($options=~/magneticfield/) {
	print "Enter the value of the magnetic field (Zeeman term).\n";
	print "Default (press ENTER): 0\n";
	$ret = <STDIN>;
	$ret = "0\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	$magneticfield=$ret;
}

print "What kind of diagonalization algorithm do you want to use.\n";
print "Available: DIAG TPEM PEM TPEM_FULLTRACE\n";
print "Default (press ENTER): DIAG\n";
$ret = <STDIN>;
$ret = "DIAG\n" if ($ret eq "" or $ret eq "\n");
chomp($ret);
if ($ret eq "DIAG") {
	$algorithm=0;
} elsif ($ret eq "TPEM") {
	$algorithm=1;
} elsif ($ret eq "PEM") {
	$algorithm=2;
} elsif ($ret eq "TPEM_FULLTRACE") {
	$algorithm=3;
} else {
	die "Algorithm $ret not available.\n";
}

if ($algorithm>0) {
	print "Enter the number of moments and the thresholds for product and trace truncation separated by spaces.\n";
	print "Example (press ENTER): 40 1e-5 1e-6\n";
	 $ret = <STDIN>;
	$ret = "40 1e-5 1e-6\n" if ($ret eq "" or $ret eq "\n");
	$tpemvalues=$ret;
}

if ($model=~/MODEL_KONDO_DMS/) {
	print "ATTENTION: You will need to edit sampleinput.inp by hand and enter the band hoppings!!\n";
	print "Future versions of input.pl will automatize this part.\n";
}

open(FILE,">sampleinput.inp") or die "Cannot open file sampleinput.inp for writing: $!\n";
print FILE<<EOF;
# Created with input.pl for SPF v6.4
# This script is in DEVELOPMENT, CHECK THIS FILE CAREFULLY FOR ERRORS
#
OPTIONS $options
LATTICETYPE $latticetype
LATTICELENGTHS $latticelengths
CARRIERS $carriers
CHEMICALPOTENTIAL $mu
$beta_label $beta
ROOTNAME $rootname
MCTHERMALIZATIONS $mctherm
MCEFFECTIVE $mceffective
MCUNMEASURED $mcunmeasured
MCWINDOW $mcwindow
MCFLAG $mcflag
MCSTARTTYPE $mcstarttype
EOF
if ($mcstarttype==4 || $mcstarttype==6) {
print "MCSTARTFILE $mcstartfile\n";
}
if ($mcstarttype==6) {
print "MCSTARTLEVEL $mcstartlevel\n";
}
print FILE <<EOF;
HISTOGRAMSTEPS $histsteps
EOF
if ($options=~/histogrambounds/) {
print "HISTOGRAMBOUNDS $histbounds\n";
}
print FILE <<EOF; 
HAMILTONIANCONCENTRATION $conc
BOUNDARYCONDITIONS $bc
EOF
if ($model=~/MODEL_KONDO_DMS/ || $model=~/MODEL_KONDO_FINITE/) {
	print FILE "HAMILTONIANJH $J\n";
}
print FILE <<EOF;
JAF $Jaf
EOF

if ($options=~/havepotential/) {
	print FILE "HAMILTONIANPOTENTIAL $potential\n";
}
if ($options=~/magneticfield/) {
	print FILE "HAMILTONIANZEEMAN $magneticfield\n";
}
print FILE "ALGORITHM $algorithm\n";
if ($algorithm>0) {
	print "TPEMVALUES $tpemvalues\n";
}
if ($model=~/MODEL_KONDO_INF_TWOBANDS/) {
	if  ($dimension==2) {
		$bandh="1 -0.577350269  -0.577350269 0.33333333 1  0.577350269  0.577350269 0.33333333 ";
	} elsif ($dimension==3) {
		$bandh=" 1 -0.577350269  -0.577350269 0.33333333 1  0.577350269  0.577350269 0.33333333 0 0 0 1.33333333 ";
	} else {
		$bandh="UNKNOWN";
		print STDERR "Warning: Unknown band hoppings form MODEL_KONDO_INF_TWOBANDS in $dimension dimensions.\n";
	}
	print FILE "BANDHOPPINGS $bandh\n";
	if ($options=~/freezephonons/) {
		$lambda=0;
		$ejt1=$ejt2=10;
	} else {
		print "Enter the value for lambda (electron-phonon coupling).\n";
		print "Default (press ENTER): 0\n";
		$ret = <STDIN>;
		$ret = "0\n" if ($ret eq "" or $ret eq "\n");
		$lambda=$ret;
		chomp($lambda);
		$ejt1 = 1;
		$ejt2 = 0.5;
	}
	print FILE <<EOF;
EJT $lambda $lambda $lambda
EDT $ejt1 $ejt2 $ejt2
EOF
}
if ($model=~/MODEL_KONDO_BCS/) {
	print "Enter the value for Delta (which will remain frozen and homogeneous).\n";
	print "Default (press ENTER): 1\n";
	$ret = <STDIN>;
	$ret = "1\n" if ($ret eq "" or $ret eq "\n");
	chomp($ret);
	print FILE "DELTA $ret 1\n";
}	

print FILE <<EOF;
#EOF
EOF
close(FILE);


print "The file sampleinput.inp has been written.\n";

sub nofspins
{
	$volume = getvolume();
	$ret = $volume;
	if ($latticetype=~/honeycomb/i || $latticetype=~/bcc/) {
		$ret = 2*$volume;
	}
	if ($latticetype=~/fcc/i) {
		$ret = 4*$volume;
	}
	return $ret;
}

sub getvolume
{
	@temp=split(/,/,$latticelengths);
	$volume = $temp[0];
	for ($i=1;$i<=$#temp;$i++) {
		$volume=$volume * $temp[$i];
	}
	return $volume;
} 
