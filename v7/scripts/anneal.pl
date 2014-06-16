#!/usr/bin/perl

use strict;
use warnings;
my ($templateInput) = @ARGV;
defined($templateInput) or die "USAGE: $0 templateInput\n"; 

my $tmp = readLabel($templateInput,"#Betas");
my @betas = split/,/,$tmp;

for (my $i = 0; $i < scalar(@betas); ++$i) {
	my $input = createInput($i,$betas[$i]);
	runInput($input);
}

sub createInput
{
	my ($i, $beta) = @_;
	my $file = "Input$i.inp";
	open(FOUT,">$file") or die "$0: Cannot write to $file\n";
	
	my $j = $i - 1;
	my $prev = ($i > 0) ? "Output$j" : "none";
	my $output = "Output$i";

	open(FILE,"$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$".$name;
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);
	return $file;
}

sub runInput
{
	my ($input) = @_;
	my $msg = "./spf -f $input";
	print STDERR "$0: About to run $msg\n";
	system("$msg &> /dev/null");
	print STDERR "$0: Finish running $msg\n";
}

	
sub readLabel
{
	my ($file,$label) = @_;
	my $val;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	while(<FILE>) {
		chomp;
		if (/^$label(.*)/) {
			$val = $1;
			last;
		}
	}

	close(FILE);

	defined($val) or die "$0: No $label in $file\n";

	return $val;
}

