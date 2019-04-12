#!/usr/bin/perl -w
use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use lib ".";
use PsiTag;

my ($flavor, $generateSources, $su2enabled, $lto) = (NewMake::noFlavor() , 0, 0, 0);
my $usage = "USAGE: $0 [-f flavor] [-s] [-su2] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           's' => \$generateSources,
           'su2' => \$su2enabled,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my @configFiles = ("../TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

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

