#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--ASCatalogList \\\n" .
                "--ASTypeList \\\n" .
                "--genomeFasta \\\n" .
		"--outputRelatedSeqFile \n";
	exit;
}

my ($ASCatalogList, $ASTypeList, $genomeFasta, $outputRelatedSeqFile);

GetOptions(
        'ASCatalogList=s'=>\$ASCatalogList,
        'ASTypeList=s'=>\$ASTypeList,
        'genomeFasta=s'=>\$genomeFasta,
        'outputRelatedSeqFile=s'=>\$outputRelatedSeqFile,
);

my (%genome, $seq, $id, @tt, $line);
# 将genome序列读入hash
open FF, "<$genomeFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>.*/){
		@tt = split(/ /, $line);
		$id = substr($tt[0],1);
	}else{
		$genome{$id}.=$line;
	}
}
close FF;

# 将所有的AS坐标读入到hash
# 依次将AS的坐标读入到hash中
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
open WW, ">$outputRelatedSeqFile";
my (%as, @asCatalogFile, $asCatalogFile, @asType, $asType, $i, @gg);
my ($asLine, @valueField, @titleField, $fieldId, $asId);
my ($longExon, $shortExon, $flankExon, $riExon, $upstreamExon, $downstreamExon, $skipingExon, $firstExon, $secondExon);
my ($changedSeq, $exchangedSeq1, $exchangedSeq2);

@asCatalogFile = split(/,/, $ASCatalogList);
@asType = split(/,/, $ASTypeList);
for($i=0; $i<=$#asCatalogFile; $i++){
	$asCatalogFile = $asCatalogFile[$i];
	$asType = $asType[$i];

	# 打开一个AS目录文件
	open FF, "<$asCatalogFile";
	# 获得列名
	$asLine = <FF>;
	chomp($asLine);
	@titleField = ();
	@titleField = split(/\t/, $asLine);
	# 逐个读取AS，将其进入hash
	while($asLine=<FF>){
		chomp($asLine);
		@valueField = ();
		@valueField = split(/\t/, $asLine);
		$asId = $valueField[0];
		for($fieldId=0; $fieldId<=$#valueField; $fieldId++){
			${$as{$asId}}{$titleField[$fieldId]} = $valueField[$fieldId];
		}
		
		if($asType eq "A3SS"){
			$longExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"longExonEnd"} - ${$as{$asId}}{"longExonStart_0base"}));
			$shortExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"shortEE"} - ${$as{$asId}}{"shortES"}));
			$flankExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"flankingEE"} - ${$as{$asId}}{"flankingES"}));

			$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"shortES"} -${$as{$asId}}{"longExonStart_0base"} ));

			if(${$as{$asId}}{"strand"} eq "-"){
				$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"longExonEnd"} - ${$as{$asId}}{"shortEE"}));
				$changedSeq=reverse($changedSeq);
				$changedSeq=~tr/ACGT/TGCA/;

				$longExon=reverse($longExon);
				$longExon=~tr/ACGT/TGCA/;
				$shortExon=reverse($shortExon);
				$shortExon=~tr/ACGT/TGCA/;
				$flankExon=reverse($flankExon);
				$flankExon=~tr/ACGT/TGCA/;
			}
			print WW "asId, longExon, shortExon, flankExon, flankExon, changedSeq____" . $asId . ", " . $longExon . ", " . $shortExon . ", " . $flankExon . ", " . $flankExon . ", " . $changedSeq . "\n";
		}

		if($asType eq "A5SS"){
			$longExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"longExonEnd"} - ${$as{$asId}}{"longExonStart_0base"}));
			$shortExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"shortEE"} - ${$as{$asId}}{"shortES"}));
			$flankExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"flankingEE"} - ${$as{$asId}}{"flankingES"}));

			$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"longExonEnd"} - ${$as{$asId}}{"shortEE"} ));

			if(${$as{$asId}}{"strand"} eq "-"){
				
				$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"shortES"} - ${$as{$asId}}{"longExonStart_0base"} ));
				$changedSeq=reverse($changedSeq);
				$changedSeq=~tr/ACGT/TGCA/;

				$longExon=reverse($longExon);
				$longExon=~tr/ACGT/TGCA/;
				$shortExon=reverse($shortExon);
				$shortExon=~tr/ACGT/TGCA/;
				$flankExon=reverse($flankExon);
				$flankExon=~tr/ACGT/TGCA/;
			}
	
			print WW "asId, longExon, shortExon, flankExon, changedSeq____" . $asId . ", " . $longExon . "," . $shortExon . ", " . $flankExon . ", " .  $changedSeq . "\n";		
		}

		if($asType eq "RI"){
			$riExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"riExonEnd"} - ${$as{$asId}}{"riExonStart_0base"}));
			$upstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"upstreamES"}, ${$as{$asId}}{"upstreamEE"} - ${$as{$asId}}{"upstreamES"}));
			$downstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"downstreamEE"} - ${$as{$asId}}{"downstreamES"}));
			$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"downstreamES"} - ${$as{$asId}}{"upstreamEE"} ));
			if(${$as{$asId}}{"strand"} eq "-"){
				$changedSeq=reverse($changedSeq);
				$changedSeq=~tr/ACGT/TGCA/;
				$riExon=reverse($riExon);
				$riExon=~tr/ACGT/TGCA/;
				$upstreamExon=reverse($upstreamExon);
				$upstreamExon=~tr/ACGT/TGCA/;
				$downstreamExon=reverse($downstreamExon);
				$downstreamExon=~tr/ACGT/TGCA/;
			}
			print WW "asId, riExon, upstreamExon, downstreamExon, changedSeq____" . $asId . ", " . $riExon . ", " . $upstreamExon . ", " . $downstreamExon . ", " . $changedSeq . "\n";
		}

		if($asType eq "SE"){
			$skipingExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"exonEnd"} - ${$as{$asId}}{"exonStart_0base"}));
			$upstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"upstreamES"}, ${$as{$asId}}{"upstreamEE"} - ${$as{$asId}}{"upstreamES"}));
			$downstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"downstreamEE"} - ${$as{$asId}}{"downstreamES"}));
			$changedSeq = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"exonEnd"} - ${$as{$asId}}{"exonStart_0base"} ));	
			if(${$as{$asId}}{"strand"} eq "-"){
				$changedSeq=reverse($changedSeq);
				$changedSeq=~tr/ACGT/TGCA/;
				$skipingExon=reverse($skipingExon);
				$skipingExon=~tr/ACGT/TGCA/;
				$upstreamExon=reverse($upstreamExon);
				$upstreamExon=~tr/ACGT/TGCA/;
				$downstreamExon=reverse($downstreamExon);
				$downstreamExon=~tr/ACGT/TGCA/;
			}
			print WW "asId, skipingExon, upstreamExon, downstreamExon, changedSeq____" . $asId . ", " . $skipingExon . ", " . $upstreamExon . ", " . $downstreamExon . ", " . $changedSeq . "\n";
		}

		if($asType eq "MXE"){
			$firstExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"1stExonEnd"} - ${$as{$asId}}{"1stExonStart_0base"}));
			$secondExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"} - ${$as{$asId}}{"2ndExonStart_0base"}));
			$upstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"upstreamES"}, ${$as{$asId}}{"upstreamEE"} - ${$as{$asId}}{"upstreamES"}));
			$downstreamExon = uc(substr($genome{${$as{$asId}}{"chr"}}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"downstreamEE"} - ${$as{$asId}}{"downstreamES"}));
			if(${$as{$asId}}{"strand"} eq "-"){
				$firstExon=reverse($firstExon);
				$firstExon=~tr/ACGT/TGCA/;
				$secondExon=reverse($secondExon);
				$secondExon=~tr/ACGT/TGCA/;
				$upstreamExon=reverse($upstreamExon);
				$upstreamExon=~tr/ACGT/TGCA/;
				$downstreamExon=reverse($downstreamExon);
				$downstreamExon=~tr/ACGT/TGCA/;
			}
			print WW "asId, firstExon, secondExon, upstreamExon, downstreamExon, upstreamExchangedSeq, downstreamExchangedSeq____" . $asId . ", " . $firstExon . ", ". $secondExon . ", " . $upstreamExon . ", " . $downstreamExon . ", " . $changedSeq . ", " . $firstExon . ", ". $secondExon . "\n";
		}
	}
	close FF;
}
close WW;

