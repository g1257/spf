#!/usr/bin/perl -w
use strict;
use File::Temp qw/ tempfile tempdir /;

my $PsimagLiteDir = "../../PsimagLite";

my $mpi = "n";
print "Do you want to use MPI?\n";
print "Available: y or n\n";
print "Default is: n (press ENTER): ";
$_=<STDIN>;
s/ //g;
chomp;
$mpi = $_ unless ($_ eq "");

my $blasAndLapack=" -lblas -llapack ";
print "What are the linking flags for BLAS and LAPACK?\n";
print "Available: any\n";
print "Default is: $blasAndLapack (press ENTER): ";
$_=<STDIN>;
s/^ +//g;
s/ +$//g;
chomp;
$blasAndLapack = $_ unless ($_ eq "");

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
# MPI: $mpi
LDFLAGS = $blasAndLapack -lm $gslLibs
EXENAME = spf
CPPFLAGS =  -Werror -Wall -Wstrict-overflow=5 -I${PsimagLiteDir} \\
            -I${PsimagLiteDir}/src   -IGeometries \\
            -IModels/PnictidesMultiOrbitals -IModels/PhononsTwoOrbitals \\
            -IModels/DmsMultiOrbital -IModels/HubbardOneOrbital -IEngine \\
            -IClassicalFields -I../Tpem $gslDefine
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

Makefile.dep: spf.cpp  Makefile
	g++ \$(CPPFLAGS) -MM spf.cpp  > Makefile.dep

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt Makefile.dep

include Makefile.dep
######## End of $thisFile ########

EOF
	close(FOUT);
	print "$thisFile has been written\n";
}

sub findIfWeHaveGsl
{
	my ($gslDefine,$gslLibs)=@_;

	my $slashTmp = "/tmp";
	return askForGslDirectly() unless (-w $slashTmp);

	my $dir = tempdir( CLEANUP => 1 );
	my ($fh, $filename) = tempfile( DIR => $dir );

	if (!$fh) {
		return askForGslDirectly();
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
	system("g++ -I$PsimagLiteDir/src $gslDefine $cppFile  $gslLibs 2>/dev/null"); 
	return 1 if (-x "a.out");
	return 0;
}

sub askForGslDirectly
{
	print "Do you have the GNU Scientific Library (GSL)?\n";
	print "Available is y or n\n";
	print "Default is n (press ENTER)";
	$_=<STDIN>;
	chomp;
	return 1 if ($_=~/^y/i);
	return 0;
}

