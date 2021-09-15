#!/usr/bin/perl
use strict;
if($#ARGV<1){
	print "$0 vcfFile   hisat2SnpFile\n";
	print "This program is used to parse vcf format in ensembl to\n";
	print "generate snp formated file for hisat2-build index\n";
	exit;
}
my ($vcfFile, $hisatSnpFile) = @ARGV;
open FF, "<$vcfFile";
open WW, ">$hisatSnpFile";
while(my $line=<FF>){
	chomp($line);
	next if($line=~/^#/);
	next if($line=~/sequence_alteration/);#discard lines 3946
	next if($line=~/substitution/);#discard lines 88
	next if($line=~/tandem_repeat/);#discard lines 1015
	next if($line=~/indel/);#discard lines 936

	my @field = ();
	@field = split(/\t/, $line);
	if($field[7]=~/=deletion;/){
		#1       341722  rs432843322     TAC     T       .       .       dbSNP_150;TSA=deletion;E_Multiple_observations
		#rs432843322     deletion        1       341721  2	
		#注意:deletion中不存在逗号分割多种类型同时存在的现象
		print WW join("\t", $field[2], "deletion", $field[0], $field[1]-1, length($field[3],) - length($field[4]), "\n");
	}elsif($field[7]=~/=SNV/){
		#1       208003  rs461501953     C       G,T     .       .       dbSNP_150;TSA=SNV
		#rs461501953.0   single  1       208002  G
		#rs461501953.1   single  1       208002  T
		my @tt = ();
		@tt = split(/,/, $field[4]);
		if($#tt>0){
			for(my $i=0; $i<=$#tt; $i++){
				print WW join("\t", $field[2] . "." . $i , "single", $field[0], $field[1]-1, $tt[$i], "\n");	
			}
		}else{
			print WW join("\t", $field[2], "single", $field[0], $field[1]-1, $tt[0], "\n");
		}
	}elsif($field[7]=~/=insertion/){
		#1       510460  rs525579234     C       CA,CTT  .       .       dbSNP_150;TSA=insertion
		#rs525579234.0   insertion       1       510459  A
		#rs525579234.1   insertion       1       510459  TT
		my @tt = ();
		@tt = split(/,/, $field[4]);
		if($#tt>0){
			for(my $i=0; $i<=$#tt; $i++){
				print WW join("\t", $field[2] . "." . $i , "insertion", $field[0], $field[1]-1, substr($tt[$i],1), "\n");
			}
		}else{
			print WW join("\t", $field[2], "insertion", $field[0], $field[1]-1, substr($tt[0],1), "\n");
		}
	}
}
close FF;
close WW;
