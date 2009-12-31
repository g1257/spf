#!/usr/bin/perl

#analyzes and plots the local spin and phononic structure 
#factors in .eps format for the $i-th MC snapshot

use strict;

my ($L, $level, $dospins)=@ARGV;

my $D=2;
my $T;
my @temp;

#arrays for phonons
my @kpipiover2;
my @k0pi;
my @kpi3piover2;

#arrays for spins
my @k00;
my @kpi0;
my @k0piover2;
my @kpipi;
my @k03piover2;
my @ktotal;

my $inputfile;

#files to be used for phonons
my $kpipiover2file;
my $k0pifile;
my $kpi3piover2file;
my $kpipiover2epsfile;
my $k0piepsfile;
my $kpi3piover2epsfile;

#files to be used for spins
my $k00file;
my $kpi0file;
my $k0piover2file;
my $kpipifile;
my $k03piover2file;
my $k00epsfile;
my $kpi0epsfile;
my $k0piover2epsfile;
my $kpipiepsfile;
my $k03piover2epsfile;

my ($i,$j,$k,$l,$m);

#$T=0.035;
for ($T=0.015;$T<0.05;$T+=0.005) {

	system("perl localUk.pl 2dL8n0.75lambda1.2jaf0.033beta${T}_1orb $L $level $dospins");
		
	if (!$dospins) {
		$inputfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.local.uq";
		
		$kpipiover2file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.ukpipiover2.dat";
		$kpipiover2epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.ukpipiover2.eps";
		
		$k0pifile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.uk0pi.dat";
		$k0piepsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.uk0pi.eps";
		
		$kpi3piover2file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.ukpi3piover2.dat";
		$kpi3piover2epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.ukpi3piover2.eps";
		
		$i=$j=$k=0;
		open(FILE,"$inputfile") or die "Cannot open file $inputfile: $!\n";
		while(<FILE>) {
			next if (/^#/);
			@temp = split;
			#for a 16x16 mesh
			$kpipiover2[$i++] = $temp[1] if($temp[0]==72);
			$k0pi[$j++] = $temp[1] if($temp[0]==128);
			$kpi3piover2[$k++] = $temp[1] if ($temp[0]==200);
			#For a 64x64 mesh
			#$kpipiover2[$i++] = $temp[1] if($temp[0]==72);
			#$k0pi[$j++] = $temp[1] if($temp[0]==128);
			#$kpi3piover2[$k++] = $temp[1] if ($temp[0]==200);
		}
		close(FILE);
		
		open(FILE,">$kpipiover2file") or die "Cannot open file $kpipiover2file: $!\n";
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $kpipiover2[$i]\n";
		}
		close(FILE);
		
		open(FILE,">$k0pifile") or die "Cannot open file $k0pifile: $!\n";
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $k0pi[$i]\n";
		}
		close(FILE);
		
		open(FILE,">$kpi3piover2file") or die "Cannot open file $kpi3piover2file: $!\n";
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $kpi3piover2[$i]\n";
		}
		close(FILE);
		
		#now plot
	
		system("perl lmap3d_data.pl $kpipiover2file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $kpipiover2epsfile");
		
		system("perl lmap3d_data.pl $k0pifile $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $k0piepsfile");
		
		system("perl lmap3d_data.pl $kpi3piover2file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $kpi3piover2epsfile");

	} else {
		$inputfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.local.sq";
			
		$k00file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk00.dat";
		$k00epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk00.eps";
		
		$kpi0file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.skpi0.dat";
		$kpi0epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.skpi0.eps";
		
		$k0piover2file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk0piover2.dat";
		$k0piover2epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk0piover2.eps";
		
		$kpipifile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.skpipi.dat";
		$kpipiepsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.skpipi.eps";
		
		$k03piover2file="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk03piover2.dat";
		$k03piover2epsfile="2dL8n0.75lambda1.2jaf0.033beta${T}_1orb.sk03piover2.eps";

        	$i=$j=$k=$l=$m=0;
		open(FILE,"$inputfile") or die "Cannot open file $inputfile: $!\n";
		while(<FILE>) {
			next if (/^#/);
			@temp = split;
			$k00[$i++] = $temp[1] if($temp[0]==0);
			$kpi0[$j++] = $temp[1] if($temp[0]==8);
			$k0piover2[$k++] = $temp[1] if ($temp[0]==64);
			$kpipi[$l++] = $temp[1] if ($temp[0]==136);
			$k03piover2[$m++] = $temp[1] if ($temp[0]==192);
		}
		close(FILE);

		open(FILE,">$k00file") or die "Cannot open file $k00file: $!\n";
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $k00[$i]\n";
		}
		close(FILE);
	
		open(FILE,">$kpi0file") or die "Cannot open file $kpi0file: $!\n"; 
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $kpi0[$i]\n";
		}
		close(FILE);
	
		open(FILE,">$k0piover2file") or die "Cannot open file $k0piover2file: $!\n"; 
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $k0piover2[$i]\n";
		}
		close(FILE);
	
		open(FILE,">$kpipifile") or die "Cannot open file $kpipifile: $!\n"; 
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $kpipi[$i]\n";
		}
		close(FILE);
	
		open(FILE,">$k03piover2file") or die "Cannot open file $k03piover2file: $!\n"; 
		for ($i=0;$i<$L**2;$i++) {
			print FILE "$i $k03piover2[$i]\n";
		}
		close(FILE);
	
		#now plot
	
		system("perl lmap3d_data.pl $k00file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $k00epsfile");
	
		system("perl lmap3d_data.pl $kpi0file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $kpi0epsfile");
	
		system("perl lmap3d_data.pl $k0piover2file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $k0piover2epsfile");

		system("perl lmap3d_data.pl $kpipifile $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $kpipiepsfile");

		system("perl lmap3d_data.pl $k03piover2file $L");
		system("gnuplot cor_map.gpl");
		system("mv temp.eps $k03piover2epsfile");
	}	
}

