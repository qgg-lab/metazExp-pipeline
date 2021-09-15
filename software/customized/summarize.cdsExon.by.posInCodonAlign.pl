#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--exonPosInCodonAlignList \\\n" .
                "--taxonIdList \\\n" .
		"--listOfSummarizedRlt \\\n" .
                "--statisticOfSummarizedRlt \n";
	exit;
}

my ($fileListOfExonPosInCodonAlign, $taxonIdList, $listOfSummarizedRlt, $statisticOfSummarizedRlt);

GetOptions(
        'fileListOfExonPosInCodonAlign=s'=>\$fileListOfExonPosInCodonAlign,
        'taxonIdList=s'=>\$taxonIdList,
        'listOfSummarizedRlt=s'=>\$listOfSummarizedRlt,
	'statisticOfSummarizedRlt=s'=>\$statisticOfSummarizedRlt,
);

# 将exonPosInCodonAlign文件读入，建立$exonPosInCodonAlign --> {$taxonId} --> {$exonId} = $geneId 的hash
my (%orthoCodonAlignPos, $orthoCodonAlignPosHref,  $orthoCodonAlignPosName);
$orthoCodonAlignPosHref = \%orthoCodonAlignPos;
my ($i, $j, @taxonId, $taxonId, @fileOfExonPosInCodonAlign, $fileOfExonPosInCodonAlign);
my ($exonId, @field, $geneId, $posBegin, $posEnd, $ogId, $strand, $chr, $exonBegin, $exonEnd);
@taxonId = split(/,/, $taxonIdList);
@fileOfExonPosInCodonAlign = split(/,/, $fileListOfExonPosInCodonAlign);
for($i=0; $i<=$#taxonId; $i++){
	$taxonId = $taxonId[$i];
	$fileOfExonPosInCodonAlign = $fileOfExonPosInCodonAlign[$i];
	
	open FF, "<$fileOfExonPosInCodonAlign";
	# 
	while(my $line=<FF>){

		($exonId, $chr, $exonBegin, $exonEnd, $strand, $geneId, $ogId, $posBegin, $posEnd) = ("", "", "", "", "", "", "", "", "");
		if($line=~/^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?):(OG\d+):(\d+)\-(\d+)\n/){
			($exonId, $chr, $exonBegin, $exonEnd, $strand, $geneId, $ogId, $posBegin, $posEnd) = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
#			print "$exonId, $chr, $exonBegin, $exonEnd, $strand, $geneId, $ogId, $posBegin, $posEnd\n";
#			<STDIN>;
			$orthoCodonAlignPosName = $ogId . ":" . $posBegin . "-" . $posEnd;
			$orthoCodonAlignPosHref->{$orthoCodonAlignPosName}->{$taxonId}->{$exonId} = $geneId;
		}
	}
	close FF;
}

# 按照codonAlign中的exonPos生成汇总信息
my (@orthoCodonAlignPos, $orthoCodonAlignPos, $listOutputLine, $statOutputLine, %orthoCodonAlignPosToExonNum, %orthoCodonAlignPosToExonIdList);
my (@exonId, $exonId);
open LIST, ">$listOfSummarizedRlt";
open STAT, ">$statisticOfSummarizedRlt";
print LIST join("\t", "orthoCodonAlignPos", @taxonId) . "\n";
print STAT join("\t", "orthoCodonAlignPos", @taxonId) . "\n";
@orthoCodonAlignPos = keys(%orthoCodonAlignPos);
foreach $orthoCodonAlignPos(@orthoCodonAlignPos){

	# 初始化hash
	%orthoCodonAlignPosToExonNum = ();
	%orthoCodonAlignPosToExonIdList = ();

	# 提取每个taxonId下的exon情况
	foreach $taxonId(@taxonId){
#		print $taxonId;
#		<STDIN>;
		@exonId = ();
		if(exists($orthoCodonAlignPosHref->{$orthoCodonAlignPos}->{$taxonId})){
			@exonId = keys(%{$orthoCodonAlignPosHref->{$orthoCodonAlignPos}->{$taxonId}});
		}else{
			@exonId = ();
		}
		$orthoCodonAlignPosToExonNum{$taxonId} = $#exonId + 1;
		if($#exonId >= 0){
			foreach $exonId(@exonId){
				if(not exists($orthoCodonAlignPosToExonIdList{$taxonId})){
					$orthoCodonAlignPosToExonIdList{$taxonId} = $orthoCodonAlignPosHref->{$orthoCodonAlignPos}->{$taxonId}->{$exonId} . ":" . $exonId;
				}else{
					$orthoCodonAlignPosToExonIdList{$taxonId} .= "," . $orthoCodonAlignPosHref->{$orthoCodonAlignPos}->{$taxonId}->{$exonId} . ":" . $exonId;
				}
			}
		}else{
			$orthoCodonAlignPosToExonIdList{$taxonId} = "";
		}
	}

	# 输出结果到list
	print LIST $orthoCodonAlignPos;
	foreach $taxonId(@taxonId){
		print LIST "\t" . $orthoCodonAlignPosToExonIdList{$taxonId};
	}
	print LIST "\n";

	# 输出结果到statistic
	print STAT $orthoCodonAlignPos;
	foreach $taxonId(@taxonId){
		print STAT "\t" . $orthoCodonAlignPosToExonNum{$taxonId};
	}
	print STAT "\n";
}
close LIST;
close STAT;
