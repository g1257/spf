#!/usr/bin/perl -w

package Average;
use strict;

sub average
{
	my ($label,$file) = @_;
	my ($sum,$counter)=(0,0);
	my $fh = *STDIN;
	my $option = 1;
	if (defined($file)) {
		open($fh,$file) or die "Cannot open $file: $!\n";
		$option = 0;
	}
	# Read standard input, while's there's any
	while(<$fh>) {
		if (/$label(.*$)/) {
			my $value = $1;
			$value=~s/,.*$//;
			$value=~s/\(//;
			print "$counter $value\n" if ($option);
			$sum += $value;
			$counter++;
		}
	}
	die "$0: label=$label\n" if ($counter==0);

	$_ = $sum/$counter; # average
	print "#$sum $counter $_\n" if ($option);
	if (defined($file)) {
		close($file);
	}
	return $_;
}
1;

