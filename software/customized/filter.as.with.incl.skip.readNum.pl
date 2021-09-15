#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputAsTsv input.as.tsv \\\n" .
                "--jcecTsvFileList jcec.A5SS.tsv,jcec.SE.tsv,jcec.RI.tsv,jcec.MXE.tsv,jcec.A3SS.tsv \\\n" .
                "--experimentTsvFile experiment.tsv \\\n" .
                "--minInclReadNum 3\\\n" .
                "--minSkipReadNum 3 \\\n" .
                "--minExpNum 4 \\\n" .
                "--minPrjNum 3 \\\n" .
		"--ouputAsTsv flitered.as.tsv \n";
	exit;
}

my ($inputAsTsv, $jcecTsvFileList, $experimentTsvFile, $minInclReadNum, $minSkipReadNum, $minExpNum, $minPrjNum, $ouputAsTsv);
my (@jcecTsvFile, $jcecTsvFile);
GetOptions(
        'inputAsTsv=s'=>\$inputAsTsv,
        'jcecTsvFileList=s'=>\$jcecTsvFileList,
        'experimentTsvFile=s'=>\$experimentTsvFile,
        'minInclReadNum=s'=>\$minInclReadNum,
	'minSkipReadNum=s'=>\$minSkipReadNum,
	'minExpNum=s'=>\$minExpNum,
	'minPrjNum=s'=>\$minPrjNum,
	'ouputAsTsv=s'=>\$ouputAsTsv,
);

@jcecTsvFile = split(/,/, $jcecTsvFileList);

my (%as, $asId, @asId);
my (@field, $field, $line);
my ($fieldNameString, $valueString, %expIdToPrj, $expId);

# 读取expId -> prjId的关系
open FF, "<$experimentTsvFile";
while($line = <FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

        @field = split(/, /, $valueString);
        $expId = substr($field[2], 1, length($field[2]) - 2);
        $expIdToPrj{$expId} = $field[11];
}
close FF;

# 扫描jcecTsvFile，构建AsId ->(uniqMinInclExpIdList, uniqMinInclPrjIdList, uniqMinSkipExpIdList, uniqMinSkipPrjIdList,)
foreach $jcecTsvFile(@jcecTsvFile){
	open FF, "<$jcecTsvFile";
	<FF>;
	# ASID    		  IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
	# GGALA5SS0000000001      527     	  31      	  730     	  100     	  SRX1036607
	while($line=<FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);

		if($field[1] >= $minInclReadNum){
			${$as{$field[0]}}{"minInclExpIdList"}.=$field[5] . ",";
		}
		if($field[2] >= $minSkipReadNum){
			${$as{$field[0]}}{"minSkipExpIdList"}.=$field[5] . ",";

		}

	}
	close FF;
}


# 扫描所有AsId，判断其expId数量和prjId数量是否符合要求
open WW, ">$ouputAsTsv";
open FF, "<$inputAsTsv";
# asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach, orthAsId, orthAsTaxon, TissueExpNumList, TissueNum_____"GGALMXE0000034850", "MXE", "Gallus gallus", "20", "-", 719031, 753202, 722489, 721992, 730499, 730443, 753202, 753039, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 719178, 719030, "ENSGALG00000032417", "CHD6", "753040..753202,721993..722489,719031..719178", "ERX1751702.14362.20,SRX893872.18044.18,SRX529337.12111.16", 3, "753040..753202,730444..730499,719031..719178", "", 0, "PF00176,PF00385,PF07533,PF00271", "GO:0005524,GO:0005515,GO:0016817", 1115, "SDI", "", 0, 7, 7, "Novel", "", "", "parietal nerve(1),embryos(1),muscle(2),embryo wing bud(1),gonads(2)", 5
my (@expId, $expId, %prjId, @prjId, $exprEvidence);
while($line=<FF>){
	chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

	# 设置该as的表达证据为1
	$exprEvidence = "Y";

        @field = split(/, /, $valueString);
	$asId = substr($field[0], 1, length($field[0]) - 2);

	# 检测该as的inclusion是否达标
	@expId = ();
	@expId = split(/,/, ${$as{$asId}}{"minInclExpIdList"});
	$exprEvidence = "N" if($#expId +1 < $minExpNum);

	%prjId = ();
	foreach $expId(@expId){
		$prjId{$expIdToPrj{$expId}} = 1;
	}
	@prjId = ();
	@prjId = keys(%prjId);
	$exprEvidence = "N" if($#prjId +1 < $minPrjNum);

	# 检测该as的skiping是否达标
	@expId = ();
        @expId = split(/,/, ${$as{$asId}}{"minSkipExpIdList"});
        $exprEvidence = "N" if($#expId +1 < $minExpNum);

        %prjId = ();
        foreach $expId(@expId){
                $prjId{$expIdToPrj{$expId}} = 1;
        }
        @prjId = ();
        @prjId = keys(%prjId);
        $exprEvidence = "N" if($#prjId +1 < $minPrjNum);

	print WW $fieldNameString . ", " . "splicingExpEVID" . "_____";
	print WW $valueString . ", " . "\"" . $exprEvidence . "\"\n";
}
close FF;
close WW;
