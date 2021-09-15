#!/usr/bin/perl
# plantExpGeneBlastToNcbiProtId ncbiProtId2NcbigeneId ncbigeneId2EntrezId outputPlantExpGene2EntrezId
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--plantExpGeneBlastToNcbiProtId \\\n" .
                "--ncbiProtId2NcbigeneId \\\n" .
                "--ncbigeneId2EntrezId \\\n" .
		"--outputPlantExpGene2EntrezId";
	exit;
}

my ($plantExpGeneBlastToNcbiProtId, $ncbiProtId2NcbigeneId, $ncbigeneId2EntrezId, $outputPlantExpGene2EntrezId);

GetOptions(
        'plantExpGeneBlastToNcbiProtId=s'=>\$plantExpGeneBlastToNcbiProtId,
        'ncbiProtId2NcbigeneId=s'=>\$ncbiProtId2NcbigeneId,
        'ncbigeneId2EntrezId=s'=>\$ncbigeneId2EntrezId,
        'outputPlantExpGene2EntrezId=s'=>\$outputPlantExpGene2EntrezId,
);

my (%protId2NcbiGeneId, %ncbigeneId2EntrezId);
my ($line, @tmp);

# ncbi protId 到 ncbi geneId 的映射
open FF, "<$ncbiProtId2NcbigeneId";
# XP_028961873.1  LOC103428462
# XP_017186679.2  LOC103431315
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	$protId2NcbiGeneId{$tmp[0]} = $tmp[1];
}
close FF;

# ncbi geneId 到 entrezId的映射
open FF, "<$ncbigeneId2EntrezId";
# LOC103446548    103446548
# LOC103406057    103406057
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	$ncbigeneId2EntrezId{$tmp[0]} = $tmp[1];
} 
close FF;

# 读取blast
my (%finalGene2NcbiProtId, $hr, @tmp, $geneId, $protId, $evalue);
$hr = \%finalGene2NcbiProtId;
open FF, "<$plantExpGeneBlastToNcbiProtId";
# MD09G1241100g   XP_028964330.1  100.00  345     0       0       1       345     1       345     0.0       707 
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	($geneId, $protId, $evalue) = ($tmp[0], $tmp[1], $tmp[10]);
	if(not exists($finalGene2NcbiProtId{$geneId})){
		$hr->{$geneId}->{"protId"} = $protId;
		$hr->{$geneId}->{"evalue"} = $evalue;
	}else{
		if($evalue<$hr->{$geneId}->{"evalue"}){
			$hr->{$geneId}->{"protId"} = $protId;
			$hr->{$geneId}->{"evalue"} = $evalue;
		}
	}
}
close FF;

my ($ncbiGeneId, $ncbiProtId);
my @geneId = keys(%finalGene2NcbiProtId);
open WW, ">$outputPlantExpGene2EntrezId";
foreach my $geneId(@geneId){
	$ncbiProtId = $hr->{$geneId}->{"protId"};
#	print $ncbiProtId;
#	<STDIN>;
	$ncbiGeneId = $protId2NcbiGeneId{$ncbiProtId};
#	print $ncbiGeneId;
#	<STDIN>;
	print WW join("\t", $geneId, $ncbigeneId2EntrezId{$ncbiGeneId}) . "\n";
}
close WW;
