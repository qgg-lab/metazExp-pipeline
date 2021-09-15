#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputExonPosInCdsTsv \\\n" .
                "--inputCodonAlignDir \\\n" .
                "--inputOrthoGroupList \\\n" .
                "--species \\\n" .
		"--outputExonPosInCodonAlign \n";
	exit;
}

my ($inputExonPosInCdsTsv, $species, $inputCodonAlignDir, $inputOrthoGroupList, $outputExonPosInCodonAlign);

GetOptions(
        'inputExonPosInCdsTsv=s'=>\$inputExonPosInCdsTsv,
        'inputCodonAlignDir=s'=>\$inputCodonAlignDir,
        'inputOrthoGroupList=s'=>\$inputOrthoGroupList,
	'species=s'=>\$species,
        'outputExonPosInCodonAlign=s'=>\$outputExonPosInCodonAlign,
);

# 将所有orthoGroup的编号读入数组
my (@ogId, $ogIdList, $ogId, $cmd);
$cmd = "grep -v \"Orthogroup\" $inputOrthoGroupList |awk -F \'\\t\' \'{print \$1}\'";
$ogIdList = `$cmd`;
@ogId = split(/\n/, $ogIdList);

# 将codonAlignment读取hash
my (%codonAlignmentSeq, $codonAlignmentSeqHref);
my (@seqId, $seqId, $line);
$codonAlignmentSeqHref=\%codonAlignmentSeq;
foreach $ogId(@ogId){
	open FF, "<$inputCodonAlignDir" . "/" . $ogId . ".codon.align";
	while($line=<FF>){
		chomp($line);
		if($line=~/>(.*)/){
			$seqId = $1;
			$codonAlignmentSeqHref->{$seqId}->{"orthoGroupId"} = $ogId;
		}else{
			$codonAlignmentSeqHref->{$seqId}->{"codonSeq"}.=$line;
		}
	}
	close FF;
}

# 扫描每个mergedInnerCdsExon在对应uniqPep中的位置
my (@field, @exonPosInUniqPepCds, $exonPosInUniqPepCds);
my ($uniqPepId, $exonStartInCds, $exonStopInCds);
my ($exonStartInCodonAlign, $exonStopInCodonAlign, $exonPosInCodonAlignText,  $alignPos, $geneId);
open FF, "<$inputExonPosInCdsTsv";
open WW, ">$outputExonPosInCodonAlign";
# 1       3996    4276    +       Atha_AT1G01010_pep002:155-435,Atha_AT1G01010_pep001:155-435
# 1       4486    4605    +       Atha_AT1G01010_pep002:436-555 
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	@exonPosInUniqPepCds = ();
	@exonPosInUniqPepCds = split(/,/, $field[4]);
	
	$exonPosInCodonAlignText = "";
	foreach $exonPosInUniqPepCds(@exonPosInUniqPepCds){
		if($exonPosInUniqPepCds=~/(.*?):(\d+)\-(\d+)/){
			$uniqPepId = $1;
			$exonStartInCds = $2;
			$exonStopInCds = $3;
		}
		# 检测当前uniqPepId是否存在与codonAlignment中
		if(exists($codonAlignmentSeqHref->{$uniqPepId})){
			if($uniqPepId=~/($species)_(.*?)_(pep\d+)/){
				$geneId = $2;
			}
			$alignPos=&getExonPosInCodonAlign($codonAlignmentSeqHref->{$uniqPepId}->{"orthoGroupId"}, $codonAlignmentSeqHref->{$uniqPepId}->{"codonSeq"}, $geneId, $exonStartInCds, $exonStopInCds);
			if($exonPosInCodonAlignText eq ""){
				$exonPosInCodonAlignText = $alignPos;
			}else{
				$exonPosInCodonAlignText .= "," . $alignPos;
			}
		}
	}
	if($exonPosInCodonAlignText ne ""){
		print WW join("\t", $field[0], $field[1], $field[2], $field[3], $exonPosInCodonAlignText) . "\n";
	}
}
close FF;
close WW;

# 计算exon在codon alignment fasta中的位置
sub getExonPosInCodonAlign{
	my ($ogId, $codonAlignSeq, $geneId, $exonStartInCds, $exonStopInCds) = @_;
	my ($i, $baseNum, $exonStartInCodonAlign, $exonStopInCodonAlign);
	# codonAlignSeq:
	# >Atha_AT3G49601_pep001
	# ATGTACAACGGAATAGGGTTACAGACAGCAAGAGGATCAGGAACTAATGGTTATGTTCAG
	# ACGAATAAGTTTTTTGTGAGA------CCTAGAAATGGTGGTAAGCCTGTTAAAGGTGGG
	# AAAGGTTTTGAG---------------------------GATGATGAAGGTACTGCTGGT
	$baseNum = 0;
	$exonStartInCodonAlign = 0;
	$exonStopInCodonAlign = 0;
	for($i=0; $i<length($codonAlignSeq); $i++){
		if(substr($codonAlignSeq, $i, 1) ne "-"){
			$baseNum++;
			$exonStartInCodonAlign = ($i + 1) if($baseNum == $exonStartInCds);
			$exonStopInCodonAlign = ($i + 1) if($baseNum == $exonStopInCds);
		}
	}
	return $geneId . ":" . $ogId . ":" . $exonStartInCodonAlign . "-" . $exonStopInCodonAlign;
}
