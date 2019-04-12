#!/usr/bin/perl -w
use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use lib ".";
use PsiTag;

my ($flavor, $lto) = (NewMake::noFlavor() , 0, 0, 0);
my $usage = "USAGE: $0 [-f flavor] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my @configFiles = ("TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

my @drivers = ("spf");

my %args;
$args{"CPPFLAGS"} = $lto;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "SPFv7";
$args{"configFiles"} = \@configFiles;

createMakefile(\@drivers, \%args);

sub createMakefile
{
        my ($drivers, $args) = @_;
        NewMake::backupMakefile();

        my $fh;
        open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

        NewMake::main($fh, $args, $drivers);

	print "$0: Done writing Makefile\n";
	close($fh);
}


