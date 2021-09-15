#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--gtfFile \\\n" . 
                "--inputTrsptBlastTrsptRlt \\\n" .
                "--inputTrsptWithGo \\\n" .
                "--outputTrsptGoAnno \\\n" .
		"--outputGeneGoAnno \n";
	exit;
}

my ($gtfFile, $inputTrsptBlastTrsptRlt, $inputTrsptWithGo, $outputTrsptGoAnno, $outputGeneGoAnno);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
        'inputTrsptBlastTrsptRlt=s'=>\$inputTrsptBlastTrsptRlt,
        'inputTrsptWithGo=s'=>\$inputTrsptWithGo,
        'outputTrsptGoAnno=s'=>\$outputTrsptGoAnno,
	'outputGeneGoAnno=s'=>\$outputGeneGoAnno,
);


# set trspt to gene mapping
my (%trsptToGene, %geneGoAnno, $line, $line, @field, @att);
my ($geneId, $trsptId);
open FF, "<$gtfFile";
while($line=<FF>){
	($geneId, $trsptId) = ("", "");
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "exon");
	@att = ();
	@att = split(/;/, $field[8]);
	for(my $i=0; $i<=$#att; $i++){
		if($att[$i]=~/gene_id "(.*)"/){
			$geneId = $1;
		}
		if($att[$i]=~/transcript_id "(.*)"/){
			$trsptId = $1;
		}
	}
	if($trsptId ne "" and $geneId ne ""){
		$trsptToGene{$trsptId} = $geneId;
	}
}
close FF;

my (%targetTrsptGoList);
# read go annot to trspt
open FF, "<$inputTrsptWithGo";
# LOC_Os01g01010.1        GO:0030234      F       enzyme regulator activity       IEA     TAIR:AT3G59570
# LOC_Os01g01010.1        GO:0007165      P       signal transduction     IEA     TAIR:AT3G59570
while($line=<FF>){
	@field = split(/\t/, $line);
	$targetTrsptGoList{$field[0]} .= $field[1] . ",";
}
close FF;


# build queryTrspt to target mapping
my (%queryTrsptToTargetTrspt, $queryTrsptToTargetTrsptHref, $queryTrspt, $targetTrspt, $evalue);
$queryTrsptToTargetTrsptHref = \%queryTrsptToTargetTrspt;
open FF, "<$inputTrsptBlastTrsptRlt";
while($line=<FF>){
	# Os01t0968400-00 LOC_Os01g73730.1        100.000 1181    0       0       1       1181    90      1270    0.0     2182
	# Os03t0231650-00 LOC_Os03g12890.3        96.438  1151    21      12      1       1147    1433    299     0.0     1881
	@field = split(/\t/, $line);
	($queryTrspt, $targetTrspt, $evalue) = ($field[0], $field[1], $field[10]);
	if(not exists($queryTrsptToTargetTrsptHref->{$queryTrspt})){
		$queryTrsptToTargetTrsptHref->{$queryTrspt}->{"targetTrspt"} = $targetTrspt;
		$queryTrsptToTargetTrsptHref->{$queryTrspt}->{"evalue"} = $evalue;
	}elsif(exists($queryTrsptToTargetTrsptHref->{$queryTrspt}) and $queryTrsptToTargetTrsptHref->{$queryTrspt}->{"evalue"} > $evalue){
		$queryTrsptToTargetTrsptHref->{$queryTrspt}->{"targetTrspt"} = $targetTrspt;
		$queryTrsptToTargetTrsptHref->{$queryTrspt}->{"evalue"} = $evalue;
	}
}
close FF;

# ouytput trspt with go anno
my %queryTrsptGoAnno;
open WW, ">$outputTrsptGoAnno";
my @trsptId = keys(%queryTrsptToTargetTrspt);
foreach $queryTrspt(@trsptId){
	$targetTrspt = $queryTrsptToTargetTrsptHref->{$queryTrspt}->{"targetTrspt"};
	print WW join("\t", $queryTrspt, substr($targetTrsptGoList{$targetTrspt}, 0, length($targetTrsptGoList{$targetTrspt}) - 1)) . "\n";
	$queryTrsptGoAnno{$queryTrspt} = substr($targetTrsptGoList{$targetTrspt}, 0, length($targetTrsptGoList{$targetTrspt}) - 1);
}
close FF;

# 获得gene到transcript的列表
@trsptId = ();
@trsptId = keys(%trsptToGene);
my %geneToTrsptList;
foreach $trsptId(@trsptId){
	$geneToTrsptList{$trsptToGene{$trsptId}} .= $trsptId . ",";
}

my @geneId = keys(%geneToTrsptList);
my (%geneGoAnno);
foreach $geneId(@geneId){
	@trsptId = ();
	@trsptId = split(/,/, $geneToTrsptList{$geneId});
	foreach $trsptId(@trsptId){
		if(not exists($queryTrsptGoAnno{$trsptId})){
			next;
		}else{
			if(not exists($geneGoAnno{$geneId})){
				$geneGoAnno{$geneId} = $queryTrsptGoAnno{$trsptId};
			}elsif(length($queryTrsptGoAnno{$trsptId}) > $geneGoAnno{$geneId}){
				$geneGoAnno{$geneId} = $queryTrsptGoAnno{$trsptId};
			}
		}		
	}	
}

open WW, ">$outputGeneGoAnno";
@geneId = ();
@geneId = keys(%geneGoAnno);
foreach $geneId(@geneId){
	next if($geneGoAnno{$geneId} eq "");
	print WW join("\t", $geneId, $geneGoAnno{$geneId}) . "\n";
}
close WW;
