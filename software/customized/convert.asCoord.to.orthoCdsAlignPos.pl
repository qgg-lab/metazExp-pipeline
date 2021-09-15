#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoCdsAlignPosFile \\\n" .
                "--inputAsFileList \\\n" .
                "--AsTypeList \\\n" .
                "--outputAsFileList \n";
	exit;
}

my ($orthoCdsAlignPosFile, $inputAsFileList, $AsTypeList, $outputAsFileList);

GetOptions(
        'orthoCdsAlignPosFile=s'=>\$orthoCdsAlignPosFile,
        'inputAsFileList=s'=>\$inputAsFileList,
        'AsTypeList=s'=>\$AsTypeList,
        'outputAsFileList=s'=>\$outputAsFileList,
);

# 将orthoCdsAlignPos对应的外显子位置读入hash
my (%orthoCdsAlignPos, $orthoCdsAlignPosHref);
$orthoCdsAlignPosHref = \%orthoCdsAlignPos;
my ($line, @titleField, @valueField, $geneId, $exonId, $chr, $start, $stop, $strand, $orthoCdsAlignPos, $conservLevel, $duplication);
# print $orthoCdsAlignPosFile;
# <STDIN>;
open FF, "<$orthoCdsAlignPosFile";
# geneId  	  exonId  	  chr     start   stop    strand  orthoCdsAlignPos        conservLevel    duplication
# AT1G01020       AthaExon000001  1       7157    7232    -       orthoPos108404  	  Capp    	  2
# AT1G01020       AthaExon000002  1       7384    7450    -       orthoPos124639  	  Capp    	  2
<FF>;
while($line=<FF>){
	chomp($line);
	($geneId, $exonId, $chr, $start, $stop, $strand, $orthoCdsAlignPos, $conservLevel, $duplication) = split(/\t/, $line);
	($orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"chr"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"start"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"stop"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"strand"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"orthoCdsAlignPos"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"conservLevel"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"duplication"}) = ($chr, $start, $stop, $strand, $orthoCdsAlignPos, $conservLevel, $duplication);
#	print join("\t", $geneId, $exonId, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"chr"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"start"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"stop"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"strand"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"orthoCdsAlignPos"}) . "\n";
#	<STDIN>;
}
close FF;



# 读取每个AS的坐标，将其转换为orthoCdsAlignPos上
my (@asType, $asType, @asFile, $asFile, @outputAsFile, $outputAsFile, $i, $j);
my (%as, $asHref, %newAs, $newAsHref);
$asHref=\%as;
$newAsHref = \%newAs;
@asType = split(/,/, $AsTypeList);
@asFile = split(/,/, $inputAsFileList);
@outputAsFile = split(/,/, $outputAsFileList);
for($i=0; $i<=$#asType; $i++){
	$asType = $asType[$i];
	$asFile = $asFile[$i];
	$outputAsFile = $outputAsFile[$i];
	open FF, "<$asFile";
	open WW, ">$outputAsFile";
	# A5SS:ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
	# A3SS:ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
	# SE:  ASID GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
	# RI:  ASID GeneID geneSymbol chr strand riExonStart_0base  riExonEnd upstreamES upstreamEE downstreamES downstreamEE
	# MXE: ASID GeneID geneSymbol chr strand 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE 
	$line = <FF>;
	print WW $line;
	chomp($line);
	@titleField = split(/\t/, $line);
	while($line=<FF>){
		# ATHAA5SS0000003780      "AT1G01020"     "ARV1"  1       -       7761    7987    7941    7987    7563    7649
		chomp($line);
		@valueField = split(/\t/, $line);
		# 将当前的as读入到hash as中
		%as = ();
		for($j=0; $j<=$#valueField; $j++){
			if($valueField[$j]=~/"(.*)"/){
				$as{$titleField[$j]} = $1;
			}else{
				$as{$titleField[$j]} = $valueField[$j];
			}
		}
		
		# 清空结果hash
		%newAs = ();
		# 转换结果保存到结果hash newAs中
		&convertAsCoordinates($asHref, $orthoCdsAlignPosHref, $newAsHref);

		# 将newAs中的结果输出
		if($asType eq "A5SS" or $asType eq "A3SS"){
			print WW join("\t", $newAsHref->{"ASID"}, $newAsHref->{"GeneID"}, $newAsHref->{"geneSymbol"}, $newAsHref->{"chr"}, $newAsHref->{"strand"}, $newAsHref->{"longExonStart_0base"}, $newAsHref->{"longExonEnd"}, $newAsHref->{"shortES"}, $newAsHref->{"shortEE"}, $newAsHref->{"flankingES"}, $newAsHref->{"flankingEE"}) . "\n";;
		}elsif($asType eq "SE"){
			print WW join("\t", $newAsHref->{"ASID"}, $newAsHref->{"GeneID"}, $newAsHref->{"geneSymbol"}, $newAsHref->{"chr"}, $newAsHref->{"strand"}, $newAsHref->{"exonStart_0base"}, $newAsHref->{"exonEnd"}, $newAsHref->{"upstreamES"}, $newAsHref->{"upstreamEE"}, $newAsHref->{"downstreamES"}, $newAsHref->{"downstreamES"}) . "\n";
		}elsif($asType eq "RI"){
			print WW join("\t", $newAsHref->{"ASID"}, $newAsHref->{"GeneID"}, $newAsHref->{"geneSymbol"}, $newAsHref->{"chr"}, $newAsHref->{"strand"}, $newAsHref->{"riExonStart_0base"}, $newAsHref->{"riExonEnd"}, $newAsHref->{"upstreamES"}, $newAsHref->{"upstreamEE"}, $newAsHref->{"downstreamES"}, $newAsHref->{"downstreamES"}) . "\n";
		}elsif($asType eq "MXE"){
			print WW join("\t", $newAsHref->{"ASID"}, $newAsHref->{"GeneID"}, $newAsHref->{"geneSymbol"}, $newAsHref->{"chr"}, $newAsHref->{"strand"}, $newAsHref->{"1stExonStart_0base"}, $newAsHref->{"1stExonEnd"}, $newAsHref->{"2ndExonStart_0base"}, $newAsHref->{"2ndExonEnd"}, $newAsHref->{"upstreamES"}, $newAsHref->{"upstreamEE"}, $newAsHref->{"downstreamES"}, $newAsHref->{"downstreamEE"}). "\n";
		}
	}
	close FF;
	close WW;
}

# 转换AS的坐标到orthoCdsAlignPos上
# 1stExonEnd,1stExonStart_0base,2ndExonEnd,2ndExonStart_0base,ASID,chr,downstreamEE,downstreamES,exonEnd,exonStart_0base,flankingEE,flankingES,GeneID,geneSymbol,longExonEnd,longExonStart_0base,riExonEnd,riExonStart_0base,shortEE,shortES,strand,upstreamEE,upstreamES
sub convertAsCoordinates{
	my ($asHref, $orthoCdsAlignPosHref, $newAsHref) = @_;

	# 首先将基本信息填充到newAs中	
	($newAsHref->{"ASID"}, $newAsHref->{"GeneID"}, $newAsHref->{"geneSymbol"}, $newAsHref->{"chr"}, $newAsHref->{"strand"}) = ($asHref->{"ASID"}, $asHref->{"GeneID"}, $asHref->{"geneSymbol"}, $asHref->{"chr"}, $asHref->{"strand"});

	# 将asHref的相关坐标转换到以orthoCdsAlignPos编号为基准的坐标
	$newAsHref->{"1stExonEnd"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"1stExonEnd"}) if(exists($asHref->{"1stExonEnd"}));
	$newAsHref->{"1stExonStart_0base"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"1stExonStart_0base"} + 1) if(exists($asHref->{"1stExonStart_0base"}));

	$newAsHref->{"2ndExonEnd"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"2ndExonEnd"}) if(exists($asHref->{"2ndExonEnd"}));
	$newAsHref->{"2ndExonStart_0base"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"2ndExonStart_0base" + 1}) if(exists($asHref->{"2ndExonStart_0base"}));

	$newAsHref->{"downstreamEE"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"downstreamEE"}) if(exists($asHref->{"downstreamEE"}));
	$newAsHref->{"downstreamES"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"downstreamES"} + 1) if(exists($asHref->{"downstreamES"}));

	$newAsHref->{"exonEnd"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"exonEnd"}) if(exists($asHref->{"exonEnd"}));
	$newAsHref->{"exonStart_0base"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"exonStart_0base"} + 1) if(exists($asHref->{"exonStart_0base"}));

	$newAsHref->{"flankingEE"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"flankingEE"}) if(exists($asHref->{"flankingEE"}));
	$newAsHref->{"flankingES"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"flankingES"} + 1) if(exists($asHref->{"flankingES"}));

	$newAsHref->{"longExonEnd"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"longExonEnd"}) if(exists($asHref->{"longExonEnd"}));
	$newAsHref->{"longExonStart_0base"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"longExonStart_0base"} + 1) if(exists($asHref->{"longExonStart_0base"}));

	$newAsHref->{"riExonEnd"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"riExonEnd"}) if(exists($asHref->{"riExonEnd"}));
	$newAsHref->{"riExonStart_0base"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"riExonStart_0base"} + 1) if(exists($asHref->{"riExonStart_0base"}));

	$newAsHref->{"shortEE"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"shortEE"}) if(exists($asHref->{"shortEE"}));
	$newAsHref->{"shortES"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"shortES"} + 1) if(exists($asHref->{"shortES"}));

	$newAsHref->{"upstreamEE"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"upstreamEE"}) if(exists($asHref->{"upstreamEE"}));
	$newAsHref->{"upstreamES"} = &getCoordinateInOrthoCdsPos($asHref->{"GeneID"}, $orthoCdsAlignPosHref, $asHref->{"upstreamES"} + 1) if(exists($asHref->{"upstreamES"}));
}

# 定位单个坐标在orthoCdsAlignPos中的位置
sub getCoordinateInOrthoCdsPos{
	my ($geneId, $orthoCdsAlignPosHref, $coordinate) = @_;
	my (@exonId, $exonId, $locatedExonId);
	$locatedExonId = "";
	if(not exists($orthoCdsAlignPosHref->{$geneId})){
		return "-";
	}
	# geneId      exonId      chr     start   stop    strand  orthoCdsAlignPos        conservLevel    duplication
	# AT1G01020       AthaExon000001  1       7157    7232    -       orthoPos108404      Capp        2
	@exonId = ();
	@exonId = keys(%{$orthoCdsAlignPosHref->{$geneId}});
#	print "coordinate:$coordinate:\n";
	foreach $exonId(@exonId){
#		print join("\t", $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"start"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"stop"}, $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"orthoCdsAlignPos"}) . "\n";
		if($coordinate >= $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"start"} and $coordinate <= $orthoCdsAlignPosHref->{$geneId}->{$exonId}->{"stop"}){
			$locatedExonId = $exonId;
		}
	}
#	<STDIN>;
	# 未能定位到任何exon上
	if($locatedExonId  eq ""){
		return "-";
	}

	# 定位到某个exon上
	my $posInOrthoCdsAlignPos = 0;
	if($orthoCdsAlignPosHref->{$geneId}->{$locatedExonId}->{"strand"} eq "+"){
		$posInOrthoCdsAlignPos = $coordinate - $orthoCdsAlignPosHref->{$geneId}->{$locatedExonId}->{"start"} + 1;
	}else{
		$posInOrthoCdsAlignPos = $orthoCdsAlignPosHref->{$geneId}->{$locatedExonId}->{"stop"} - $coordinate + 1;
	}

	return $orthoCdsAlignPosHref->{$geneId}->{$locatedExonId}->{"orthoCdsAlignPos"} . ":" . $posInOrthoCdsAlignPos;
}
