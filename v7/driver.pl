#!/usr/bin/perl -w
use strict;

my $model = "PnictidesTwoOrbitals";
print "What model do you want to compile?\n";
print "Available: DmsMultiOrbital  PhononsTwoOrbitals  PnictidesTwoOrbitals\n";
print "Default is: PnictidesTwoOrbitals (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$model = $_ unless ($_ eq "");

my $geometry = "Square";
print "What geometry do you want to use?\n";
print "Available: Square Cubic Fcc\n";
print "Default is: Square (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$geometry = $_ unless ($_ eq "");

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
LDFLAGS = -L.   -llapack -lblas -lm -L../lib 
EXENAME = spf
CPPFLAGS = -DNDEBUG -I../../PsimagLite/src   -IGeometries -IModels/$model -IEngine -IClassicalFields 
CXX = g++ -Werror -Wall -O3 -pg

all: \$(EXENAME)

\$(EXENAME): clean main.o 
	\$(CXX) -o \$(EXENAME) main.o \$(LDFLAGS) 

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt

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
#include "SimpleReader.h"
#include "ParametersEngine.h"
#include "Engine.h"
#include "ConcurrencySerial.h"
#include "Parameters$model.h"
#include "$model.h"
#include "Geometry$geometry.h"
#include "Random48.h"
#include "AlgorithmDiag.h"
#include "GreenFunction.h"

typedef double FieldType;
typedef Spf::ParametersEngine<FieldType> ParametersEngineType;
typedef PsimagLite::ConcurrencySerial<FieldType> ConcurrencyType;
typedef Spf::Geometry$geometry<FieldType> GeometryType;
typedef Spf::Parameters$model<FieldType> ParametersModelType;
typedef Spf::$model<ParametersEngineType,ParametersModelType,GeometryType> ModelType;
typedef ModelType::DynVarsType DynVarsType;
typedef PsimagLite::Random48<FieldType> RandomNumberGeneratorType;
typedef Spf::AlgorithmDiag<ParametersEngineType,ModelType,RandomNumberGeneratorType> AlgorithmType;
typedef Spf::GreenFunction<ParametersEngineType,AlgorithmType> GreenFunctionType;
typedef Spf::Engine<ParametersEngineType,AlgorithmType,ModelType,ConcurrencyType,RandomNumberGeneratorType,GreenFunctionType> EngineType;

 
int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	ParametersEngineType engineParams;
	ParametersModelType mp;
	Spf::SimpleReader reader(argv[1]);
	reader.load(engineParams);
	reader.load(mp);
	// print license
	std::string license = "Copyright (c) 2009 , UT-Battelle, LLC\\n"
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
	
	ModelType model(engineParams,mp,geometry);
	AlgorithmType algorithm(engineParams,model);
	EngineType engine(engineParams,model,algorithm,concurrency);
	
	engine.main();
}


/* ####### End of $thisFile ######## */

EOF
	close(FOUT);
	print "$thisFile has been written\n";
}
