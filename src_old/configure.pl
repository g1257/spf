#!/usr/bin/perl -w

# some definitions
$gsldir="\$\{HOME\}/software/gsl-";


 
# Configure script for spf v6
print "************************************\n";
print "* Welcome to configure.pl for SPF v6.4\n";
print "************************************\n\n";
print "Please enter the name of your platform\n";
print "Available: linux,aix,sgi,osx,crayx1,catamount or tiger\n";
print "Default (press ENTER): linux\n";
$ret=<STDIN>;
$ret = "linux\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$platform=$ret;

print "Will you use the TPEM algorithm as serial or parallel (mpi)?\n";
print "(Note: If you will not use the TPEM then answer serial.)\n";
print "Available: serial or parallel\n";
print "Default (press ENTER): serial\n";
$ret=<STDIN>;
$ret = "serial\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$runtype2=$ret;

$launcher = "serial_launcher";

$runtype="no";

if ($runtype2 eq "serial") {
	print "Do you want to use trivial parallelization (old mpi_wrapper)?\n";
	print "(undocumented).\n";
print "Available: yes or no\n";
print "Default (press ENTER): no\n";
$ret=<STDIN>;
$ret = "no\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$runtype=$ret;
}

if ($runtype=~/yes/i) {
	$runtype = "parallel";
	$launcher = "parallel_launcher";
} else {
	$runtype = "serial";
}

print "What type of model would you like to compile? (spell correctly!!)\n";
print "Available: MODEL_KONDO_INF_ONEBAND, MODEL_KONDO_INF_ONEBAND_PHONONS, MODEL_KONDO_INF_TWOBANDS MODEL_KONDO_FINITE ";
print " MODEL_KONDO_DMS_ZINCBLENDE or MODEL_KONDO_DMS_CUBIC or MODEL_KONDO_DMS_THREEBANDS or ";
print " MODEL_KONDO_BCS or MODEL_KONDO_DMS_FCC \n";
print "Default (press ENTER): MODEL_KONDO_FINITE\n";
$ret=<STDIN>;
$ret = "MODEL_KONDO_FINITE\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$model=$ret;

print "Do you want to use the GSL?\n";
print "Available: yes or no\n";
print "Default (press ENTER): no\n";
$ret=<STDIN>;
$ret = "no\n" if ($ret eq "" or $ret eq "\n");
chop($ret);
$gsl=$ret;
if ($gsl=~/yes/) {
	print "What version of the GSL do you want to use?\n";
	print "Available: any\n";
	print "Default (press ENTER): 1.5\n";
	$ret=<STDIN>;
	$ret="1.5\n"  if ($ret eq "" or $ret eq "\n");
	chop($ret);
	$gsldir="$gsldir$ret";
	print "Where is the GSL installed?\n";
	print "Available: any\n";
	print "Default (press ENTER): $gsldir\n";
	$ret=<STDIN>;
	$ret="$gsldir\n"  if ($ret eq "" or $ret eq "\n");
	chop($ret);	
	$gsldir=$ret;
	print "I will assume your GSL directory is in $gsldir\n";
}

$objectfilelist="filelist.make";



# setup variables
#system("cat header.txt > model.cpp");
#open(MFILE,">model.cpp") or die "Cannot open file model.cpp for writing: $!\n";
#print MFILE<<EOF;
#/* DO NOT EDIT!! CHANGES WILL BE LOST!! Edit kondo_*.cpp instead
#* File created by configure.pl
#* By G.A. */
#EOF


if ($model eq "MODEL_KONDO_INF_TWOBANDS") {
	$exename="spf2b";
	$modelcpp="kondo_inf_twobands";
	
} elsif ($model eq "MODEL_KONDO_INF_ONEBAND") {
	$exename="spf1b";
	$modelcpp="kondo_inf_oneband";
} elsif ($model eq "MODEL_KONDO_FINITE") {
	$exename="spf";
	$modelcpp="kondo_finite";
} elsif ($model eq "MODEL_KONDO_DMS_ZINCBLENDE") { 
	$exename="spf_dms_zb";
	$modelcpp="kondo_dms_zincblende";
} elsif ($model eq "MODEL_KONDO_DMS_CUBIC") {
	$exename="spf_dms_cubic";
	$modelcpp="kondo_dms_cubic";
} elsif ($model eq "MODEL_KONDO_DMS_THREEBANDS") {
	$exename="spf3b";
	$modelcpp="kondo_dms_threebands";
} elsif ($model eq "MODEL_KONDO_INF_ONEBAND_PHONONS") {
	$exename="spf1b_phonons";
	$modelcpp="kondo_inf_oneband_phonons";
} elsif ($model eq "MODEL_KONDO_BCS") {
	$exename = "spfbcs";
	$modelcpp="kondo_bcs";
} elsif ($model eq "MODEL_KONDO_DMS_FCC") {
	$exename = "spfdmsfcc";
	$modelcpp="kondo_dms_fcc";
} else {
	die "Model $model has not been implemented (sorry).\n";
}

$usempi="";
$libs="";
$defines="-I../PartialPsimag ";
$platformid="";
$warning="-Wall";
$gsllib="-lgsl_partial";
$cxx="g++";
$cc="gcc";
$cppoptions="-O2";
$lapack="-llapack -lblas -lm";
$libinc="-L. ";

if ($model eq "MODEL_KONDO_DMS_CUBIC" or $model eq "MODEL_KONDO_DMS_ZINCBLENDE" or $model eq "MODEL_KONDO_DMS_FCC" ) {
	$defines=$defines." -DMODEL_KONDO_DMS_MANYBANDS ";
}

if ($model=~/PHONONS/) {
	$defines="$defines -DMODEL_KONDO_PHONONS ";
}

#if ($objectfilelist eq "mesoscondlist.make") {
#	$exename="conductance";
#}

if ($runtype=~/parallel/i or $runtype2=~/parallel/i) {
	$usempi="-DUSE_MPI" if ($runtype2=~/parallel/i);
	$libs="-lmpi" if ($platform=~/sgi/ || $platform=~/crayx1/i);
	$exename=$exename."_mpi";
}

if ($platform=~/linux/i) {
	$lapack=$lapack." ";
	if ($runtype=~/parallel/i || $runtype2=~/parallel/i) {
		$cxx="mpicxx";
		$cc="mpicc";
	}
} elsif ($platform=~/sgi/i) {
	$lapack="-lscs  -lm";
	$platformid="-D__sgi";
	$cxx=$cc="icpc";
} elsif ($platform=~/aix/i) {
	$cxx="xlC";
	$cc="xlc";
	$lapack=" -lm -lessl";
	$gsllib="-lgsl_partial_aix";
	if ($runtype=~/parallel/i || $runtype2=~/parallel/i) {
		$cxx="mpCC";
		$cc="mpcc";
	}
	$warning="";
} elsif ($platform=~/crayx1/i) {
	$cxx="CC";
	$cc="cc";
	$lapack="-lm -lsci";
	$gsllib="-lgsl_partial_crayx1";
	$warning="";
} elsif ($platform=~/tiger/i) {
	$cxx="/opt/pgi6.0.5/linux86-64/6.0/bin/pgCC";
	$cc="undefined";
	$defines="$defines -I/usr/mpich/mpich-1.2.6/include";
	$libinc="$libinc -L/lib64 -L/opt/pgi6.0/linux86-64/6.0/lib";
	$lapack="$lapack  ";
	$gsllib="-lgsl_partial_tiger";
	$warning="";
	if ($runtype=~/parallel/i || $runtype2=~/parallel/i) {
		$libinc="$libinc -L/usr/mpich/mpich-1.2.6/lib";
		$libs="$libs -lrapl -lpthread -lmpich ";
	}
} elsif ($platform=~/catamount/i) {
	$warning="";
	$gsllib="-lgsl_partial_rizzo";	
	$libs="$libs -lpgftnrtl  ";
	$cxx = "pgCC -fastsse";
	$cc="pgcc -fastsse";
	if ($runtype=~/parallel/i || $runtype2=~/parallel/i) {
		$cxx="mpicxx  -fastsse -Mipa=fast -Minfo=all";
		$cc="mpicc";
		
	}
} else {
	die "Program has not been ported to this platform: $platform\n";
}

# Recap
$libs="$libs $lapack -L../lib ";
$libs="$libs $gsllib " if ($gsl=~/yes/);
$defines="$defines $platformid $usempi -D$model";
$defines="$defines -I$gsldir " if ($gsl=~/yes/);
$defines="$defines -DNO_GSL " if (!($gsl=~/yes/));
if ($platform=~/tiger/i) { # This must be at the end of the libs for tiger
	$libs="$libs -lpgftnrtl ";
}
# Write makefile
if (-e "Makefile") {
	print "Creating back up of file Makefile in Makefile.bak...";
	system("cp Makefile Makefile.bak") if -e "Makefile";
	print "succeded.\n";
}


open(FILE,">Makefile") or die "Cannot open file Makefile for writing: $!\n";
print FILE<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl 
# SPF v6.4 by G.A.
# Platform: $platform
# RunType for TPEM: $runtype2
# RunType for Trivial: $runtype
# Model: $model
# GSL: $gsl
LDFLAGS = $libinc $libs
EXENAME = $exename
CPPFLAGS = $defines $warning $cppoptions
CXX = $cxx 
include $objectfilelist

OBJECTS += $modelcpp.o
OBJECTS += $launcher.o
all: \$(EXENAME)

tpemSample: tpemSample.o tpemplus.o tpemplusSparse.o tpemplusSubspace.o tpemplus.h tpemplusTypes.h
\t\$(CXX) tpemSample.o tpemplus.o tpemplusSparse.o tpemplusSubspace.o -o tpemSample \$(LDFLAGS)


\$(EXENAME): \$(HEADERS) \$(OBJECTS)
\t\$(CXX) -o \$(EXENAME) \$(OBJECTS) \$(LDFLAGS)


clean:
\trm -f core* \$(EXENAME) *.o *.ii *.tt

\$(EXENAME).tar: *.cpp *.h Makefile* *.pl *.inp *.c
\ttar -cvf \$(EXENAME).tar  *.cpp *.h Makefile* *.pl *.inp *.c

######## End of Makefile ########

EOF
close(FILE);

print "configure.pl: File \"Makefile\" has been generated.\n";
$make="make";
$make="gmake" if ($platform=~/aix/i);
print "configure.pl: Cleaning up...\n";
system("$make clean");
#print "configure.pl: Running make...\n";
#system("$make");
print "Now you have to run: \n";
print "$make\n";


