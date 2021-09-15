#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV <0){
	print "\n\tperl $0 \\\n" . 
	" \t\t--exptInfoList filtered.experiment.list.tsv \\\n" .
	" \t\t--exptReadNumCutoff exptReadNumCutoff.tsv \\\n" .
	" \t\t--jcecList jcec.A3ss,jcec.A5SS,jcec.MXE,jcec.RI,jcec.SE \\\n" . 
	" \t\t--ASCatalogList A3SS.catalog,A5SS.catalog,MXE.catalog,RI.catalog,SE.catalog \\\n" .
	" \t\t--outputInclusionTagFile total.inclusion.with.tissue.expt.tsv \\\n" .
	" \t\t--outputSkipingTagFile total.skiping.with.tissue.expt.tsv \n";
	exit(0);
}

my ($exptInfoList, $exptReadNumCutoff, $jcecList, $ASCatalogList, $outputInclusionTagFile, $outputSkipingTagFile);
my (@jcecFile, @ASCatalogFile, @outputInclusionFile, @outputSkipingFile);
GetOptions(
	'exptInfoList=s'=>\$exptInfoList,,
	'exptReadNumCutoff=s'=>\$exptReadNumCutoff, 
	'jcecList=s'=>\$jcecList,
	'ASCatalogList=s'=>\$ASCatalogList,
	'outputInclusionTagFile=s'=>\$outputInclusionTagFile,
	'outputSkipingTagFile=s'=>\$outputSkipingTagFile,
);

@jcecFile = split(/,/, $jcecList);
@ASCatalogFile = split(/,/, $ASCatalogList);

# 将experiment对应的组织信息读入到hash中
my (%experiment, $exptId);
my ($experimentLine, @fieldName, @field);
open SAMPLE, "<$exptInfoList";
$experimentLine = <SAMPLE>;
chomp($experimentLine);
@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}
while($experimentLine=<SAMPLE>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);
        $exptId = $field[$exptIdPos];
        for(my $i=0; $i<=$#field; $i++){
                ${$experiment{$exptId}}{$fieldName[$i]} = $field[$i];
        }
}
close SAMPLE;

#
# 读取exptReadNumCutoff文件，将对应experiment所对应的cuttoff值读入到hash中
open CUTOFF, "<$exptReadNumCutoff";
# experiment	  Base    inclu   skiping
# SRX1073647      3.48    5.80    5.80
while($experimentLine=<CUTOFF>){
	chomp($experimentLine);
	@field = ();
	@field = split(/\t/, $experimentLine);
	$exptId = $field[0];
	${$experiment{$exptId}}{"inclusionCutoff"} = $field[2]; 
	${$experiment{$exptId}}{"skipingCutoff"} = $field[3];
}
close CUTOFF;

# 将AS的ID分别读进inclusion和skiping两个hash
my (%asInclusion, %asSkiping, $asLine);
foreach my $asCatalogFile(@ASCatalogFile){
	open FF, "<$asCatalogFile";
	<FF>;
	while($asLine=<FF>){
		@field = ();
		@field = split(/\t/, $asLine);
		${$asInclusion{$field[0]}}{"totalExptNum"} = 0;
		${$asSkiping{$field[0]}}{"totalExptNum"} = 0;
	}
	close FF;
}

# 读取jcec文件中AS及其在experiment中的inclusion和skiping对应的read数量
# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
# ZMAYA5SS0000000001      0       27      171     150     SRX2034567
# (1) 如果IJC_SAMPLE ！=0 且 IJC_SAMPLE>=experiment对应的inclusion阈值
#	获得experiment对应的tissue
#	${$asInclusion{ASID}}{"totalExptNum"}++;
#	如果${$asInclusion{ASID}}{tissue}不存在, ${$asInclusion{$field[0]}}{"totalTissueNum"}++
#	${$asInclusion{ASID}}{tissue} .= $experimentId . ","
# (2) 如果SJC_SAMPLE_1 !=0 且 SJC_SAMPLE_1》=experiment对应的skiping阈值
# 	同上。
my ($psiLine, $asId, $inclusionReadNum, $skipingReadNum, $tissue);
foreach my $jcecFile(@jcecFile){
	open FF, "<$jcecFile";
	<FF>;
	while($psiLine=<FF>){
		chomp($psiLine);
		@field = ();
		@field = split(/\t/, $psiLine);
		$asId = $field[0];
		$exptId = $field[$#field];
		$inclusionReadNum = $field[1];
		$skipingReadNum = $field[2];
		$tissue = ${$experiment{$exptId}}{"Tissue"};


		if($inclusionReadNum >= ${$experiment{$exptId}}{"inclusionCutoff"}){
			${$asInclusion{$asId}}{"totalExptNum"}+=1;
			${$asInclusion{$asId}}{$tissue} .= $exptId . ",";
		}

		if($skipingReadNum >= ${$experiment{$exptId}}{"skipingCutoff"}){
			${$asSkiping{$asId}}{"totalExptNum"}+=1;
			${$asSkiping{$asId}}{$tissue} .= $exptId . ",";
		}
	
	}
	close FF;
}


 
# 输出结果文件
open IWW, ">$outputInclusionTagFile";
my @tissue;
my @asId = keys(%asInclusion);
foreach $asId(@asId){

	@tissue = ();
	@tissue = keys($asInclusion{$asId});
	print IWW join("\t", $asId, "total:" . ${$asInclusion{$asId}}{"totalExptNum"});
	foreach $tissue(@tissue){
		next if($tissue eq "totalExptNum");
		${$asInclusion{$asId}}{$tissue} = substr(${$asInclusion{$asId}}{$tissue}, 0, length(${$asInclusion{$asId}}{$tissue})-1);
		print IWW "\t" . $tissue . ":" . ${$asInclusion{$asId}}{$tissue};
	}
	print IWW "\n";

}
close IWW;

open SWW, ">$outputSkipingTagFile";
my @tissue;
my @asId = keys(%asSkiping);
foreach $asId(@asId){
	@tissue = ();
	@tissue = keys($asSkiping{$asId});
	print SWW join("\t", $asId, "total:" . ${$asSkiping{$asId}}{"totalExptNum"});
	foreach $tissue(@tissue){
		next if($tissue eq "totalExptNum");
		${$asSkiping{$asId}}{$tissue} = substr(${$asSkiping{$asId}}{$tissue}, 0, length(${$asSkiping{$asId}}{$tissue})-1);
		print SWW "\t" . $tissue . ":" . ${$asSkiping{$asId}}{$tissue};
	}
	print SWW "\n";

}
close SWW;

