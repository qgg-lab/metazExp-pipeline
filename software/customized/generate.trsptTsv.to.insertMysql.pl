#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf final.anno.gtf \\\n" .
                "--cDNAfasta fina.cDNA.fa \\\n" .
                "--pepFasta fina.pep.fa \\\n" .
                "--orfInCdnaTsv orf.in.cDNA.tsv \\\n" .
                "--goListTsv  go.uniq.term.list.tsv \\\n" .
                "--pfamListTsv pfam.sorted.list.tsv \\\n" .
                "--asToTrsptListTsv as.with.trsptId.list.tsv \\\n" .
		"--trsptTsv trspt.mysql.tsv \n";
	exit;
}

my ($gtf, $cDNAfasta, $pepFasta, $orfInCdnaTsv, $goListTsv, $pfamListTsv, $asToTrsptListTsv, $trsptTsv);

GetOptions(
        'gtf=s'=>\$gtf,
        'cDNAfasta=s'=>\$cDNAfasta,
        'pepFasta=s'=>\$pepFasta,
        'orfInCdnaTsv=s'=>\$orfInCdnaTsv,
	'goListTsv=s'=>\$goListTsv,
	'pfamListTsv=s'=>\$pfamListTsv,
	'asToTrsptListTsv=s'=>\$asToTrsptListTsv,
	'trsptTsv=s'=>\$trsptTsv,
);

my (%trspt, $trsptHref);
$trsptHref = \%trspt;
# 将gtf读入，获得trsptId, geneId, geneSymbol, strand, trsptSymbol, chr, exonSeries和注释来源（通过StringTie）判断
my ($line, @field, $geneId, $trsptId, $source, $geneSymbol, $trsptSymbol, $exonSeries, $chr,  $strand);
open FF, "<$gtf";
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	
	$source = $field[1];
	if($source eq "StringTie"){
		$source = "ImprovedGtf";
	}else{
		$source = "EnsemblGtf";
	}
	$geneId = &getGeneIdInAttrs($field[8]);
	$trsptId = &getTrsptIdInAttrs($field[8]); 
	$geneSymbol = &getGeneSymbolInAttrs($field[8]); 
	$trsptSymbol = &getTrsptSymbolInAttrs($field[8]); 
	$chr = $field[0];
	$strand = $field[6];

	if($field[2] eq "transcript"){
		$trsptHref->{$trsptId}->{"chr"} = $chr;
		$trsptHref->{$trsptId}->{"strand"} = $strand;
		$trsptHref->{$trsptId}->{"geneId"} = $geneId;
		$trsptHref->{$trsptId}->{"geneSymbol"} = $geneSymbol;
		$trsptHref->{$trsptId}->{"trsptSymbol"} = $trsptSymbol;
		$trsptHref->{$trsptId}->{"source"} = $source;
		$trsptHref->{$trsptId}->{"exonSeries"} = "";
	}elsif($field[2] eq "exon"){
		$trsptHref->{$trsptId}->{"exonSeries"} .= $field[3] . ".." . $field[4] . ",";
	}
}
close FF;


# 打开cDNAfasta文件，读取cDNA序列
my (@trsptId);
open FF, "<$cDNAfasta";
while($line = <FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@trsptId = ();
		@trsptId = split(/ /, $trsptId);
		$trsptId = $trsptId[0];
	}else{
		$trsptHref->{$trsptId}->{"cDNAseq"}.=$line;
	}
}
close FF;

# 打开pep文件，读取pep序列文件
my (@trsptId);
open FF, "<$pepFasta";
while($line = <FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@trsptId = ();
		@trsptId = split(/ /, $trsptId);
		$trsptId = $trsptId[0];
	}else{
		$trsptHref->{$trsptId}->{"pepSeq"}.=$line
	}
}
close FF;


# 打开orf文件，读取orf、start, stop 在cDNA上的位置
my ($i, $j, @titleField, $titleField, @valueField, $valueField);
open FF, "<$orfInCdnaTsv";
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonBegin
# Zm00001d013026_T001     391     1038    391     393     1039    1041
$line = <FF>;
chomp($line);
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$trsptId = $valueField[0];
	for($i=0; $i<=$#valueField; $i++){
		$trsptHref->{$trsptId}->{$titleField[$i]} = $valueField[$i];
	}
}
close FF;

# 打开goTermList，读入go
open FF, "<$goListTsv";
# Zm00001d014489_T004     GO:0004553|GO:0005975
# Zm00001d019361_T001     GO:0043531
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptHref->{$field[0]}->{"goList"} = $field[1];
}
close FF;
# 打开pfamList文件，读入pfam
open FF, "<$pfamListTsv";
# Zm00001d014489_T004     PF00232|PF00232
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptHref->{$field[0]}->{"pfamList"} = $field[1];
}
close FF;

# 打开mapping.between.as.and.trspt.tsv，将as按照longAlt和shortAlt分别关联到转录本上
my (%tmpHash, $tmpHashHref);
$tmpHashHref=\%tmpHash;
open FF, "<$asToTrsptListTsv";
# chr     strand  asId    asType  1stExonEnd      1stExonStart_0base      2ndExonEnd      2ndExonStart_0base      downstreamEE    downstreamES
#     exonEnd exonStart_0base flankingEE      flankingES      longExonEnd     longExonStart_0base     riExonEnd       riExonStart_0base       shortEE shortES upstreamEE      upstreamES      residentType    asExonSeries    trsptId trsptExonSeries
# Mt      -       ZMAYA3SS0000023346      A3SS    -1      -1      -1      -1      -1      -1      -1      -1      88286   82759   79334   78911
#    -1      -1      79332   78911   -1      -1      longAlt 82760..88286,78912..79334       SRX2909946.38326.5      82760..88286,78912..79334,76132..76220,
$line = <FF>;
chomp($line);
@titleField = ();
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHashHref->{$titleField[$i]} = $valueField[$i];
	}
	$trsptId = $tmpHashHref->{"trsptId"};
	$trsptHref->{$trsptId}->{$tmpHashHref->{"residentType"} . "AsIdList"} .= $tmpHashHref->{"asId"} . ",";
}
close FF;

# 逐个trspt读取，并且以mysqlTsv输出
# 注意：没有对应的longAltAsIdList和shortAltAsIdList则输出-
open WW, ">$trsptTsv";
my ($fieldString, $valueString);
@trsptId = ();
@trsptId = keys(%trspt);
foreach $trsptId(@trsptId){
	$fieldString = join(", ", "trsptId", "chr", "strand", "trsptSymbol", "geneId", "geneSymbol", "orfStart", "orfStop", "startCodonStart", "startCodonStop", "stopCodonStart", "stopCodonStop", "annoOrigin", "exonSeries", "goTermList", "pfamList", "longAltAsIdList", "shortAltAsIdList", "cDNAseq", "pepSeq");
	
	if(not exists($trsptHref->{$trsptId}->{"goList"})){
		$trsptHref->{$trsptId}->{"goList"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"pfamList"})){
		$trsptHref->{$trsptId}->{"pfamList"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"longAltAsIdList"})){
		$trsptHref->{$trsptId}->{"longAltAsIdList"} = "NA";
	}else{
		$trsptHref->{$trsptId}->{"longAltAsIdList"} = substr($trsptHref->{$trsptId}->{"longAltAsIdList"}, 0, length($trsptHref->{$trsptId}->{"longAltAsIdList"}) -1);
	}
	if(not exists($trsptHref->{$trsptId}->{"shortAltAsIdList"})){
		$trsptHref->{$trsptId}->{"shortAltAsIdList"} = "NA";
	}else{
		$trsptHref->{$trsptId}->{"shortAltAsIdList"} = substr($trsptHref->{$trsptId}->{"shortAltAsIdList"}, 0, length($trsptHref->{$trsptId}->{"shortAltAsIdList"}) -1);
	}

	if(not exists($trsptHref->{$trsptId}->{"pepSeq"})){
		$trsptHref->{$trsptId}->{"pepSeq"} = "-";
	}
	if(substr($trsptHref->{$trsptId}->{"exonSeries"}, length($trsptHref->{$trsptId}->{"exonSeries"})-1, 1) eq ","){
		$trsptHref->{$trsptId}->{"exonSeries"} = substr($trsptHref->{$trsptId}->{"exonSeries"}, 0, length($trsptHref->{$trsptId}->{"exonSeries"}) - 1);
	}
	$valueString = join(", ", 
			$trsptId, $trsptHref->{$trsptId}->{"chr"}, $trsptHref->{$trsptId}->{"strand"},
			$trsptHref->{$trsptId}->{"trsptSymbol"}, $trsptHref->{$trsptId}->{"geneId"}, $trsptHref->{$trsptId}->{"geneSymbol"}, 
			$trsptHref->{$trsptId}->{"orfBegin"}, $trsptHref->{$trsptId}->{"orfEnd"}, 
			$trsptHref->{$trsptId}->{"startCodonBegin"}, $trsptHref->{$trsptId}->{"startCodonEnd"}, 
			$trsptHref->{$trsptId}->{"stopCodonBegin"}, $trsptHref->{$trsptId}->{stopCodonEnd},
			$trsptHref->{$trsptId}->{"source"}, $trsptHref->{$trsptId}->{"exonSeries"},
			$trsptHref->{$trsptId}->{"goList"}, $trsptHref->{$trsptId}->{"pfamList"},
			$trsptHref->{$trsptId}->{"cDNAseq"}, $trsptHref->{$trsptId}->{"pepSeq"},
			);

	print WW $fieldString . "___" . $valueString . "\n";
}
close WW;
sub getTrsptIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_id "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

sub getGeneIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/gene_id "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

sub getGeneSymbolInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/gene_name "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

sub getTrsptSymbolInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_name "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

