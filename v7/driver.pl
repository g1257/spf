#!/usr/bin/perl -w
use strict;
use File::Temp qw/ tempfile tempdir /;

my $PsimagLiteDir = "../../PsimagLite";
my $model = "PnictidesMultiOrbitals";
print "What model do you want to compile?\n";
print "Available: DmsMultiOrbital  PhononsTwoOrbitals  PnictidesMultiOrbitals\n";
print "Default is: PnictidesMultiOrbitals (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$model = $_ unless ($_ eq "");

my $modelFile = $model;
if ($model eq "PnictidesMultiOrbitals") {
	$modelFile = "PnictidesTwoOrbitals";
	print "How many orbitals?\n";
	print "Available: 2 or 3\n";
	print "Default is: 2 (press ENTER): ";
	$_=<STDIN>;
	s/ //g;
	chomp;
	$modelFile=~s/Two/$_/ unless($_ eq "");
	$modelFile=~s/2/Two/;
	$modelFile=~s/3/Three/;
}

my $geometry = "Square";
print "What geometry do you want to use?\n";
print "Available: Square Cubic Fcc Square45Degrees\n";
print "Default is: Square (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$geometry = $_ unless ($_ eq "");

my $mpi = "n";
print "Do you want to use MPI?\n";
print "Available: y or n\n";
print "Default is: n (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$mpi = $_ unless ($_ eq "");


my $gslDefine = " -DUSE_GSL ";
my $gslLibs = " -lgsl -lgslcblas ";
my $hasGsl = findIfWeHaveGsl($gslDefine,$gslLibs);
if (!$hasGsl) {
	$gslLibs = "  ";
	$gslDefine = "  ";
}

my $compiler = "g++";
$compiler = "mpicxx -DUSE_MPI" if ($mpi eq "y");


createMakefile();
createDriver();

sub createMakefile
{
	my $thisFile = "Makefile";
	system("mv $thisFile $thisFile.bak") if (-r "$thisFile");
	open(FOUT,">$thisFile") or die "Cannot open $thisFile for writing: $!\n";
print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify $0 instead
# This $thisFile was written by $0
# SPF v7 by G.A. and C.S.
# Platform: linux
# Model: $model
# MPI: $mpi
LDFLAGS =  -llapack -lblas -lm $gslLibs
EXENAME = spf
CPPFLAGS =  -Werror -Wall -I${PsimagLiteDir}/src   -IGeometries -IModels/$model -IEngine -IClassicalFields -I../Tpem $gslDefine
CXX = $compiler -DNDEBUG -O3 
#comment out this one for debugging
#CXX = $compiler -g3
#comment out for performance profiling with valgrind
#CXX = $compiler -g3 -O1 -DNDEBUG

all: \$(EXENAME)

\$(EXENAME):  spf.o
	\$(CXX) -o \$(EXENAME) spf.o \$(LDFLAGS)  

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	\$(CXX) \$(CPPFLAGS) -c \$<

spf.cpp: driver.pl
	perl driver.pl

Makefile.dep: spf.cpp  Makefile
	g++ \$(CPPFLAGS) -MM spf.cpp  > Makefile.dep

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt

include Makefile.dep
######## End of $thisFile ########

EOF
	close(FOUT);
	print "$thisFile has been written\n";
}

sub createDriver
{
	my $thisFile = "spf.cpp";
	system("mv $thisFile $thisFile.bak") if (-r "$thisFile");
	open(FOUT,">$thisFile") or die "Cannot open $thisFile for writing: $!\n";
print FOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify $0 instead
 This $thisFile was written by $0
 SPF v7 by G.A. and C.S.
 Platform: linux
 Model: $model
 */

typedef double FieldType;
#include "ParametersEngine.h"
#include "Engine.h"
#include <algorithm>
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<FieldType> MyConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<FieldType> MyConcurrencyType;
#endif
#include "Parameters$modelFile.h"
#include "$model.h"
#include "Geometry$geometry.h"
#include "Random48.h"

typedef PsimagLite::IoSimple::In IoInType;
typedef Spf::ParametersEngine<FieldType,IoInType,MyConcurrencyType> ParametersEngineType;
typedef Spf::Geometry$geometry<FieldType> GeometryType;
typedef Spf::Parameters$modelFile<ParametersEngineType,IoInType> ParametersModelType;
typedef Spf::$model<ParametersEngineType,ParametersModelType,GeometryType,MyConcurrencyType> ModelType;
typedef ModelType::DynVarsType DynVarsType;
typedef PsimagLite::Random48<FieldType> RandomNumberGeneratorType;
typedef Spf::Engine<ParametersEngineType,ModelType,IoInType,RandomNumberGeneratorType> EngineType;

int main(int argc,char *argv[])
{
	MyConcurrencyType concurrency(argc,argv);
	if (argc<2) {
		std::string s = "Usage is: ./" + std::string(argv[0]) +
		" input_file\\n";
		throw std::runtime_error(s.c_str());
	}
	PsimagLite::IoSimple::In io(argv[1]);
	ParametersEngineType engineParams(io,concurrency);
	ParametersModelType mp(io,engineParams);
	
	if (concurrency.root()) std::cout<<EngineType::license();
	GeometryType geometry(engineParams.latticeLength);
	
	ModelType model(engineParams,mp,geometry,concurrency);

	EngineType engine(engineParams,model,io,concurrency);

	engine.main();
}


/* ####### End of $thisFile ######## */

EOF
	close(FOUT);
	print "$thisFile has been written\n";
}

sub findIfWeHaveGsl
{
	my ($gslDefine,$gslLibs)=@_;
	my $dir = tempdir( CLEANUP => 1 );
	my ($fh, $filename) = tempfile( DIR => $dir );

	if (!$fh) {
		print "Do you have the GNU Scientific Library (GSL)?\n";
		print "Available is y or n\n";
		print "Default is n (press ENTER)";
		$_=<STDIN>;
		chomp;
		return 1 if ($_=~/^y/i);
		return 0;
	}

	#$fh or die "Cannot write to temporary filehandle: $!\n";
print $fh <<EOF;
#include "GslWrapper.h"
int main() { return 0;}
EOF
	close($fh);
	my $cppFile = $filename.".cpp";
	system("mv $filename $cppFile");
	unlink("a.out");
# 	print "Temp file is $cppFile\n";
# 	system("cat $cppFile");
	system("g++ -I$PsimagLiteDir/src $gslDefine $cppFile  $gslLibs &> /dev/null"); 
	return 1 if (-x "a.out");
	return 0;
}

