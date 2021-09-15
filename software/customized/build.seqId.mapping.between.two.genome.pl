#!/usr/bin/perl
use strict;
use Getopt::Long; 
if($#ARGV < 0){
	print "\tperl $0 \\\n" .
		"\t\t --firstGenomeFile  ensemblGenome.fa \\\n" . 
		"\t\t --secondGenomeFile refSeqGenome.fa  \\\n" .
		"\t\t --maxDiffBaseNum 100 \\\n" .
		"\t\t --maxDiffBaseRatio 0.0001 \\\n" .
		"\t\t --seqLenChangeRegion 5 \\\n" .
		"\t\t --outputMismatchFile mismatchBase.tsv\\\n" .
		"\t\t --outputSeqIdMappingFile GenomeSeqIdMapping.tsv \\\n\n";

	print "This script is used to find the mapping relationship of seqs in two genome\n";
	print "The mapping sequences in two genomes are not required to perfect match each other.\n\n";

	exit(0);
}
my ($firstGenomeFile, $secondGenomeFile, $outputSeqIdMappingFile, $maxDiffBaseNum, $maxDiffBaseRatio, $seqLenChangeRegion, $outputMismatchFile);
my (%firstGenomeSeq, %secondGenomeSeq);

$maxDiffBaseNum = 100;
$seqLenChangeRegion = 5;

GetOptions(
        'firstGenomeFile=s'=>\$firstGenomeFile,
        'secondGenomeFile=s'=>\$secondGenomeFile,
	'maxDiffBaseNum=s'=>\$maxDiffBaseNum,
	'maxDiffBaseRatio=s'=>\$maxDiffBaseRatio,
	'seqLenChangeRegion=s'=>\$seqLenChangeRegion,
	'outputMismatchFile=s'=>\$outputMismatchFile,
        'outputSeqIdMappingFile=s'=>\$outputSeqIdMappingFile,
);


my ($line, $id, @tt, $seqId, @seqId);
open FF, "<$firstGenomeFile";
#gi|417531923|ref|NC_019467.1|
##>ref|NC_037337.1| Bos taurus ...
##>AmTr_v1.0_scaffold00001 length=15980527
while($line = <FF>){
	chomp($line);
	if($line=~/>/){
		@tt = ();
		@tt = split(/ /, $line);
		$id = substr($tt[0], 1);
		if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
			$id = $1;
		}
	}else{
		${$firstGenomeSeq{$id}}{"seq"} .=uc($line);
	}
}
close FF;
@seqId = ();
@seqId = keys(%firstGenomeSeq);
foreach $seqId(@seqId){
	${$firstGenomeSeq{$seqId}}{"len"} = length(${$firstGenomeSeq{$seqId}}{"seq"});
	${$firstGenomeSeq{$seqId}}{"mapSeqId"} = "-";
}

#print "finsh loading first genome into hash!\n";

open FF, "<$secondGenomeFile";
#gi|417531923|ref|NC_019467.1|
#>ref|NC_037337.1| Bos taurus ...
#>AmTr_v1.0_scaffold00001 length=15980527
while($line = <FF>){
	chomp($line);
	if($line=~/>/){
		@tt = ();
		@tt = split(/ /, $line);
		if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
			$id = $1;
		}
	}else{
		${$secondGenomeSeq{$id}}{"seq"} .=uc($line);
	}
}
close FF;

@seqId = ();
@seqId = keys(%secondGenomeSeq);
foreach $seqId(@seqId){
	${$secondGenomeSeq{$seqId}}{"len"} = length(${$secondGenomeSeq{$seqId}}{"seq"});
	${$secondGenomeSeq{$seqId}}{"mapSeqId"} = "-";
	${$secondGenomeSeq{$seqId}}{"mappingRlt"} = "$seqId\t-\t-\t-\t-";
}



#print "finsh loading second genome into hash!\n";

my ($flag);
my @firstGenomeSeqId = keys %firstGenomeSeq;
my @secondGenomeSeqId = keys %secondGenomeSeq;

# Launch the first round to seek the mapping seqId by perfect match
for(my $i=0; $i<=$#secondGenomeSeqId; $i++){
	$flag = 0;
	for(my $j=0; $j<=$#firstGenomeSeqId; $j++){
		if(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"seq"} eq ${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"seq"} ){

			${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mappingRlt"} = join("\t", $secondGenomeSeqId[$i], $firstGenomeSeqId[$j], "0", "0", "0");

			${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mapSeqId"} = $firstGenomeSeqId[$j];
			${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"mapSeqId"} = $secondGenomeSeqId[$i];

			last;

		}
	}
}

#print "finsh the first round mapping!\n";

open MISMATCH, ">$outputMismatchFile";
my $mismatchTxt = "";
# Launch the second round to seek the mapping seqId 
# in specified different bases and different base ratio
my $diffBaseNum =0;
for(my $i=0; $i<=$#secondGenomeSeqId; $i++){
	for(my $j=0; $j<=$#firstGenomeSeqId; $j++){
		if(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mapSeqId"} eq "-" and ${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"mapSeqId"} eq "-"
			and abs(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"len"} - ${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"len"}) <= $seqLenChangeRegion){
			# calculate diff base num
			$diffBaseNum = 0;
			$mismatchTxt = "";
			for(my $k=0; $k<=${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"len"}-1; $k++){
				if(substr(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"seq"}, $k, 1) ne substr(${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"seq"}, $k, 1)){
					$mismatchTxt .= join("\t", $secondGenomeSeqId[$i], $k+1, substr(${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"seq"}, $k, 1) . " -> " . substr(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"seq"}, $k, 1)) . "\n";
					$diffBaseNum++;
				}
			}
			
			# judge wether build mapping relationship
			my $ratio = $diffBaseNum/${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"len"};
			if($diffBaseNum <= $maxDiffBaseNum and $ratio <= $maxDiffBaseRatio){
				${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mappingRlt"} = join("\t", $secondGenomeSeqId[$i], $firstGenomeSeqId[$j], $diffBaseNum, $ratio, abs(${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"len"} - ${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"len"}));	
				${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mapSeqId"} = $firstGenomeSeqId[$j];
				${$firstGenomeSeq{$firstGenomeSeqId[$j]}}{"mapSeqId"} = $secondGenomeSeqId[$i];
				print MISMATCH $mismatchTxt;
				last;
			}
		}
	}
}

close MISMATCH;

# output mapping relationship
open WW, ">$outputSeqIdMappingFile";
for(my $i=0; $i<=$#secondGenomeSeqId; $i++){
	print WW ${$secondGenomeSeq{$secondGenomeSeqId[$i]}}{"mappingRlt"} . "\n";
}
close WW;
