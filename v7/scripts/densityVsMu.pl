#/usr/bin/perl -w

use strict;

#file that contains the results:
my ($file)=@ARGV;

#location of the average.pl script:
my $average = "scripts/average.pl";

#get electrons first
my @electrons;
my $labelElectrons="\"PnictidesTwoOrbitals: Number_Of_Electrons=\"";
getObservable(\@electrons,$labelElectrons);

#get mus next
my @mus;
my $labelMu="\"Adjustments: mu=\"";
getObservable(\@mus,$labelMu);


#print them side by side
my $total = $#mus+1;
printOneVersusTheOther(\@mus,\@electrons,$total);

sub getObservable
{
	my ($obs,$label)=@_;
	#print "$label\n";
	open(PIPE,"perl $average $label < $file |") or die "$0: Cannot open pipe: $!\n";

	my $counter=0;
	while(<PIPE>) {
		chomp;
		my @temp=split;
		$obs->[$counter++]=$temp[1];
	}
	close(PIPE);
}

sub printOneVersusTheOther
{
	my ($one,$theOther,$total)=@_;
	for (my $i=0;$i<$total;$i++) {
		print "$one->[$i] $theOther->[$i]\n";
	}
}
