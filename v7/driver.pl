#!/usr/bin/perl -w
use strict;

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
print "Available: Square Cubic Fcc\n";
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
LDFLAGS = -L.   -llapack -lblas -lm -L../lib 
EXENAME = spf
CPPFLAGS = -DNDEBUG -I../../PsimagLite/src   -IGeometries -IModels/$model -IEngine -IClassicalFields 
CXX = $compiler -Werror -Wall -O3 -pg

all: \$(EXENAME)

\$(EXENAME):  main.o
	\$(CXX) -o \$(EXENAME) main.o \$(LDFLAGS)  

# dependencies brought about by Makefile.dep
main.o:
	\$(CXX) \$(CPPFLAGS) -c main.cpp

main.cpp: driver.pl
	perl driver.pl

Makefile.dep: main.cpp
	\$(CXX) \$(CPPFLAGS) -MM main.cpp  > Makefile.dep

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
	my $thisFile = "main.cpp";
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
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<FieldType> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<FieldType> ConcurrencyType;
#endif
#include "Parameters$modelFile.h"
#include "$model.h"
#include "Geometry$geometry.h"
#include "Random48.h"
#include "AlgorithmDiag.h"
#include "GreenFunction.h"

typedef PsimagLite::IoSimple::In IoInType;
typedef Spf::ParametersEngine<FieldType,IoInType> ParametersEngineType;
typedef Spf::Geometry$geometry<FieldType> GeometryType;
typedef Spf::Parameters$modelFile<ParametersEngineType,IoInType> ParametersModelType;
typedef Spf::$model<ParametersEngineType,ParametersModelType,GeometryType,ConcurrencyType> ModelType;
typedef ModelType::DynVarsType DynVarsType;
typedef PsimagLite::Random48<FieldType> RandomNumberGeneratorType;
typedef Spf::AlgorithmDiag<ParametersEngineType,ModelType,RandomNumberGeneratorType> AlgorithmType;
typedef Spf::GreenFunction<ParametersEngineType,AlgorithmType> GreenFunctionType;
typedef Spf::Engine<ParametersEngineType,AlgorithmType,ModelType,ConcurrencyType,RandomNumberGeneratorType,GreenFunctionType> EngineType;

 
int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	if (argc<2) {
		std::string s = "Usage is: ./" + std::string(argv[0]) +
		" input_file\\n";
		throw std::runtime_error(s.c_str());
	}
	PsimagLite::IoSimple::In io(argv[1]);
	ParametersEngineType engineParams(io);
	ParametersModelType mp(io,engineParams);
	// print license
	std::string license = "Copyright (c) 2009-2011, UT-Battelle, LLC\\n"
"All rights reserved\\n"
"\\n"
"[SPF, Version 7.0.0]\\n"
"\\n"
"*********************************************************\\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\\n"
"CONTRIBUTORS \\"AS IS\\" AND ANY EXPRESS OR IMPLIED\\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. \\n"
"\\n"
"Please see full open source license included in file LICENSE.\\n"
"*********************************************************\\n"
"\\n"
"\\n"
"// END LICENSE BLOCK\\n"
;
	if (concurrency.root()) std::cerr<<license;
	GeometryType geometry(engineParams.latticeLength);
	
	ModelType model(engineParams,mp,geometry,concurrency);
	AlgorithmType algorithm(engineParams,model);
	EngineType engine(engineParams,model,algorithm,concurrency);
	
	engine.main();
}


/* ####### End of $thisFile ######## */

EOF
	close(FOUT);
	print "$thisFile has been written\n";
}
