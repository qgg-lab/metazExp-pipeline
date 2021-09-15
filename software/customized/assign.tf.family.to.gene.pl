#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--blastRlt \\\n" .
                "--tfFamily \\\n" .
                "--output \n";
	exit;
}

my ($blastRlt, $tfFamily, $output);

GetOptions(
        'blastRlt=s'=>\$blastRlt,
        'tfFamily=s'=>\$tfFamily,
        'output=s'=>\$output,
);

my (%tfInfo, @tt, $seqId, $line);
# 读取转录因子的家族信息
open FF, "<$tfFamily";
# Species                         Protein ID              Family  Genome/EST      Category
# Oryza sativa subsp. japonica    LOC_Os02g27060.1        ARID    Genome          Other
# Oryza sativa subsp. japonica    LOC_Os02g27060.2        ARID    Genome          Other
<FF>;
while($line=<FF>){
	chomp($line);
	@tt = split(/\t/, $line);
	$tfInfo{$tt[1]} = $tt[2];
}
close FF;

#读取blast比对，以转录因子为传递，为每个gene赋予家族信息
# Oryza_sativa_subsp._japonica_LOC_Os01g01290.1   Os01g0102400    100.00  421     0       0       1       421     1       421     0.0       874
# Oryza_sativa_subsp._japonica_LOC_Os01g01312.1   Os01g0102800    100.00  917     0       0       271     1187    1       917     0.0      1913
# Oryza_sativa_subsp._japonica_LOC_Os01g01430.1   Os01g0104200    100.00  339     0       0       2       340     1       339     0.0       711 
my (%geneTfInfo,  $tfId,  $geneId ); 
open FF, "<$blastRlt";
while($line=<FF>){
	@tt = split(/\t/, $line);
	if($tt[0]=~/.*_(LOC_.*)/){
		$tfId = $1;	
		$geneId = $tt[1];
		$geneTfInfo{$geneId} = $tfInfo{$tfId};
	}
}
close FF;


# 输出geneId
my @geneId = keys(%geneTfInfo);
open WW, ">$output";
print WW join("\t", "geneId", "Family") . "\n";
foreach $geneId(@geneId){
	print WW join("\t", $geneId, $geneTfInfo{$geneId}) . "\n";
}
close WW;
