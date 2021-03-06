#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf \\\n" .
                "--asProcessTrsptTsv \\\n" .
                "--origCdnaFasta \\\n" .
                "--origPepOfCdna \\\n" .
                "--origOrfInCdna \\\n" .
                "--origPfamInCdna \\\n" .
                "--origGoTermList \\\n" .
		"--trsptMysqlTsv \n";
	exit;
}

my ($gtf, $asProcessTrsptTsv, $origCdnaFasta, $origPepOfCdna, $origOrfInCdna, $origPfamInCdna, $origGoTermList, $trsptMysqlTsv);

GetOptions(
        'gtf=s'=>\$gtf,
        'asProcessTrsptTsv=s'=>\$asProcessTrsptTsv,
        'origCdnaFasta=s'=>\$origCdnaFasta,
        'origPepOfCdna=s'=>\$origPepOfCdna,
	'origOrfInCdna=s'=>\$origOrfInCdna,
	'origPfamInCdna=s'=>\$origPfamInCdna,
	'origGoTermList=s'=>\$origGoTermList,
	'trsptMysqlTsv=s'=>\$trsptMysqlTsv,
);

my (%trspt, $trsptHref);
$trsptHref=\%trspt;

my (@seqId, $seqId, $line);
# 将cDNA序列读入
open FF, "<$origCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@seqId = ();
		@seqId = split(/ /, $1);
		$seqId = $seqId[0];
	}else{
		$trsptHref->{$seqId}->{"cDNAseq"} .= $line;
	}
}
close FF;

# 将pep序列读入
open FF, "<$origPepOfCdna";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@seqId = ();
		@seqId = split(/ /, $1);
		$seqId = $seqId[0];
	}else{
		$trsptHref->{$seqId}->{"pepSeq"} .= $line;
	}
}
close FF;

# 读取gtf中source
my ($cmd, $cmdRltTxt, @line, $source, $trsptId, $geneId, $geneSymbol, $trsptSymbol, $line, @field);
open FF, "<$gtf";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");
	$source = $field[1];
	($geneId, $trsptId, $geneSymbol, $trsptSymbol) = ("", "", "NA", "NA");
	&getGeneIdAndTrspt($field[8], \$geneId, \$trsptId, \$geneSymbol, \$trsptSymbol);
	$trsptHref->{$trsptId}->{"geneId"} = $geneId;
	$trsptHref->{$trsptId}->{"trsptId"} = $trsptId;
	$trsptHref->{$trsptId}->{"geneSymbol"} = $geneSymbol;
	$trsptHref->{$trsptId}->{"trsptSymbol"} = $trsptSymbol;
	if($source eq "StringTie"){
		$trsptHref->{$trsptId}->{"source"} = "improvedAnnoGtf";
	}else{
		$trsptHref->{$trsptId}->{"source"} = "origAnnoGtf";
	}
}

# 将orf信息读入
my ($orfBeg, $orfEnd, $startCodonBeg, $startCodonEnd, $stopCodonBeg, $stopCodonEnd);
open FF, "<$origOrfInCdna";
# trsptId	  orfBeg  orfEnd  start   start   stop    stop
# AT3G53400.1     420     1817    420     422     1818    1820
# AT5G45110.2     515     2062    515     517     2063    2065 
while($line=<FF>){
	chomp($line);
	($trsptId, $orfBeg, $orfEnd, $startCodonBeg, $startCodonEnd, $stopCodonBeg, $stopCodonEnd) = split(/\t/, $line);
	$trsptHref->{$trsptId}->{"orfBeg"} = $orfBeg;
	$trsptHref->{$trsptId}->{"orfEnd"} = $orfEnd;
	$trsptHref->{$trsptId}->{"startCodonBeg"} = $startCodonBeg;
	$trsptHref->{$trsptId}->{"startCodonEnd"} = $startCodonEnd;
	$trsptHref->{$trsptId}->{"stopCodonBeg"} = $stopCodonBeg;
	$trsptHref->{$trsptId}->{"stopCodonEnd"} = $stopCodonEnd;
}
close FF;

my ($pfamPosList);
# 将pfamInCdna读入
open FF, "<$origPfamInCdna";
# AT3G11450.1     PF00226[547,786]|PF00249[2020,2160]
# AT5G58690.3     PF09279[226,450]|PF00388[487,915]|PF00387[1171,1431]|PF00168[1495,1803]
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamPosList) = split(/\t/, $line);
	$trsptHref->{$trsptId}->{"pfamPosList"} = $pfamPosList;
}
close FF;

# 将go注释读入
my ($goList);
open FF, "<$origGoTermList";
# AT3G08840.13    GO:0005524|GO:0008716|GO:0046872
# AT1G49720.2     GO:0003700|GO:0006355 
while($line=<FF>){
	chomp($line);
	($trsptId, $goList) = split(/\t/, $line);
	$trsptHref->{$trsptId}->{"goList"} = $goList;
}
close FF;


# 将as和trspt之间的关系读入
# chr     strand  asId    asType  trsptId residentType    editSite        cutSize insertSize      insertSeq
# 1       -       ATHAA3SS0000002088      A3SS    SRX853408.3.6   longAlt 305     152     0
# 1       -       ATHAA3SS0000005840      A3SS    AT1G01020.3     longAlt 560     139     0
my (@nameField, @valueField, $nameField, $valueField, %tmpHash, $i, $asId);
open FF, "<$asProcessTrsptTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$trsptId = $tmpHash{"trsptId"};
	$asId = $tmpHash{"asId"};

	if($tmpHash{"residentType"} eq "longAlt"){
		if(not exists($trsptHref->{$trsptId}->{"inclusion" . $tmpHash{"asType"}})){
			$trsptHref->{$trsptId}->{"inclusion" . $tmpHash{"asType"}} = $asId ;
		}else{
			$trsptHref->{$trsptId}->{"inclusion" . $tmpHash{"asType"}} .= "," . $asId;
		}
	}else{
		if(not exists($trsptHref->{$trsptId}->{"exclusion" . $tmpHash{"asType"}})){
			$trsptHref->{$trsptId}->{"exclusion" . $tmpHash{"asType"}} = $asId ;
		}else{
			$trsptHref->{$trsptId}->{"exclusion" . $tmpHash{"asType"}} .= "," . $asId;
		}
	}
}
close FF;

 

# 将trspt信息输出到mysqlTSV
my ($nameFieldString, $valueFieldString, @trsptId);
open WW, ">$trsptMysqlTsv";
@trsptId = keys(%trspt);
foreach $trsptId(@trsptId){
	
	next if(not exists($trsptHref->{$trsptId}->{"cDNAseq"}));

######## nameFieldString ##############
	$nameFieldString = join(", ", 
		"trsptId", 
		"trsptSymbol", 
		"geneId", 
		"geneSymbol", 
		"annoSource", 
		"orfStart", 
		"orfStop", 
		"startCodonStart", 
		"startCodonStop", 
		"stopCodonStart", 
		"stopCodonStop", 
		"goTermList", 
		"pfamPosList", 
		"inclusionA3SS", 
		"inclusionA5SS", 
		"inclusionMXE", 
		"inclusionRI", 
		"inclusionSE", 
		"exclusionA3SS", 
		"exclusionA5SS", 
		"exclusionMXE", 
		"exclusionRI", 
		"exclusionSE", 
		"cDNAseq", 
		"pepSeq"
);

	if(not exists($trsptHref->{$trsptId}->{"goList"})){
		$trsptHref->{$trsptId}->{"goList"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"pfamPosList"})){
		$trsptHref->{$trsptId}->{"pfamPosList"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"inclusionA3SS"})){
		$trsptHref->{$trsptId}->{"inclusionA3SS"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"inclusionA5SS"})){
		$trsptHref->{$trsptId}->{"inclusionA5SS"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"inclusionMXE"})){
		$trsptHref->{$trsptId}->{"inclusionMXE"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"inclusionRI"})){
		$trsptHref->{$trsptId}->{"inclusionRI"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"inclusionSE"})){
		$trsptHref->{$trsptId}->{"inclusionSE"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"exclusionA3SS"})){
		$trsptHref->{$trsptId}->{"exclusionA3SS"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"exclusionA5SS"})){
		$trsptHref->{$trsptId}->{"exclusionA5SS"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"exclusionMXE"})){
		$trsptHref->{$trsptId}->{"exclusionMXE"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"exclusionRI"})){
		$trsptHref->{$trsptId}->{"exclusionRI"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"exclusionSE"})){
		$trsptHref->{$trsptId}->{"exclusionSE"} = "NA";
	}

	if(not exists($trsptHref->{$trsptId}->{"pepSeq"})){
		$trsptHref->{$trsptId}->{"pepSeq"} = "-";
	}

	# 如果没有找打orf，那么设置为-1
	if(not exists($trsptHref->{$trsptId}->{"orfBeg"})){
		$trsptHref->{$trsptId}->{"orfBeg"} = -1;
	}
	if(not exists($trsptHref->{$trsptId}->{"orfEnd"})){
		$trsptHref->{$trsptId}->{"orfEnd"} = -1;
	}
	if(not exists($trsptHref->{$trsptId}->{"startCodonBeg"})){
		$trsptHref->{$trsptId}->{"startCodonBeg"} = -1;
	}
	if(not exists($trsptHref->{$trsptId}->{"startCodonEnd"})){
		$trsptHref->{$trsptId}->{"startCodonEnd"} = -1;
	}
	if(not exists($trsptHref->{$trsptId}->{"stopCodonBeg"})){
		$trsptHref->{$trsptId}->{"stopCodonBeg"} = -1;
	}
	if(not exists($trsptHref->{$trsptId}->{"stopCodonEnd"})){
		$trsptHref->{$trsptId}->{"stopCodonEnd"} = -1;
	}

	$trsptHref->{$trsptId}->{"trsptSymbol"}=~tr/ /_/;
	$trsptHref->{$trsptId}->{"geneSymbol"}=~tr/ /_/;

######### valueFieldString ##############################
	$valueFieldString = join(", ", 
		$trsptId, 
		$trsptHref->{$trsptId}->{"trsptSymbol"}, 
		$trsptHref->{$trsptId}->{"geneId"}, 
		$trsptHref->{$trsptId}->{"geneSymbol"}, 
		$trsptHref->{$trsptId}->{"source"}, 
		$trsptHref->{$trsptId}->{"orfBeg"}, 
		$trsptHref->{$trsptId}->{"orfEnd"}, 
		$trsptHref->{$trsptId}->{"startCodonBeg"}, 
		$trsptHref->{$trsptId}->{"startCodonEnd"}, 
		$trsptHref->{$trsptId}->{"stopCodonBeg"}, 
		$trsptHref->{$trsptId}->{"stopCodonEnd"}, 
		$trsptHref->{$trsptId}->{"goList"}, 
		$trsptHref->{$trsptId}->{"pfamPosList"},  
		$trsptHref->{$trsptId}->{"inclusionA3SS"}, 
		$trsptHref->{$trsptId}->{"inclusionA5SS"}, 
		$trsptHref->{$trsptId}->{"inclusionMXE"}, 
		$trsptHref->{$trsptId}->{"inclusionRI"}, 
		$trsptHref->{$trsptId}->{"inclusionSE"}, 
		$trsptHref->{$trsptId}->{"exclusionA3SS"}, 
		$trsptHref->{$trsptId}->{"exclusionA5SS"}, 
		$trsptHref->{$trsptId}->{"exclusionMXE"}, 
		$trsptHref->{$trsptId}->{"exclusionRI"}, 
		$trsptHref->{$trsptId}->{"exclusionSE"}, 
		$trsptHref->{$trsptId}->{"cDNAseq"}, 
		$trsptHref->{$trsptId}->{"pepSeq"}
);

	print WW $nameFieldString . "___" . $valueFieldString . "\n";
}
close WW;

sub getGeneIdAndTrspt{
	my ($attrString, $geneId, $trsptId, $geneSymbol, $trsptSymbol) = @_;
	my (@attr, $attr);

	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}
		if($attr=~/gene_name "(.*)"/){
			$$geneSymbol = $1;
		}
		if($attr=~/transcript_name "(.*)"/){
			$$trsptSymbol = $1;
		}

	
	}
}
