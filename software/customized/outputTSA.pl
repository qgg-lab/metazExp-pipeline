#!/usr/bin/perl
use strict;
my $vcfFile = $ARGV[0];
open FF, "<$vcfFile";
while(my $line=<FF>){
#1       207431  rs525766619     AG      A       .       .       dbSNP_150;TSA=deletion;E_Multiple_observations
#1       207453  rs449980494     A       T       .       .       dbSNP_150;TSA=SNV
	chomp($line);
	next if($line=~/##/);
	my @filed = ();
	@filed = split(/\t/, $line);
	my @tt = ();
	@tt = split(/;/, $filed[7]);
	for(my $i=0; $i<=$#tt; $i++){
		if($tt[$i]=~/TSA=(.*)/){
			print $1 . "\n";
			last;
		}
	}
}
close FF;
