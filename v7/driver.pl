#!/usr/bin/perl -w
use strict;
use lib '../../PsimagLite/scripts';
use Make;

my @drivers = ("spf");

my $lapack = Make::findLapack();
my ($gslC,$glsL) = Make::findGsl();
Make::backupMakefile();
writeMakefile();

sub writeMakefile
{
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $libs = "$lapack  $glsL  -lm  -lpthread -lpsimaglite";
	my $cxx = "g++";
	my $cppflags = " -O3 -DNDEBUG  -I../Tpem -IEngine -I../../PsimagLite ";
	$cppflags .= " -I../../PsimagLite/src  $gslC";
	Make::make($fh,\@drivers,"spf","Linux",0,$libs,$cxx,$cppflags,"true"," "," ");

	close($fh);
	print "$0: Done writing Makefile\n";
}

