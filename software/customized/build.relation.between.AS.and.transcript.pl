#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--ASCatalogList\\\n" .
		"--ASTypeList \\\n" .
		"--gtf " . 
		"--improvementGtf \\\n" .
		"--outputTagFile \n";
        exit;
}

my ($ASCatalogList, $ASTypeList, $gtf, $asMappingTrsptTsv);
GetOptions(
	'ASTypeList=s'=>\$ASTypeList,
	'ASFileList=s'=>\$ASCatalogList,
	'gtf=s'=>\$gtf,
	'asMappingTrsptTsv=s'=>\$asMappingTrsptTsv,
);

my ($featureLine, @cols, $geneId, @isoformId, $isoformId, $isoformSymbol, $residentType);
my (%trspt, $trsptHref);
$trsptHref = \%trspt;

# 将gtf中的转录本exon串读入hash
open FF, "<$gtf";
while($featureLine = <FF>){
        chomp($featureLine);
        next if($featureLine=~/#/);
        @cols = split(/\t/, $featureLine);

        $isoformId = &getIsoformIdInAttrs($cols[8]);
        $geneId = &getGeneIdInAttrs($cols[8]);
        if($cols[2] eq "transcript"){
                $trsptHref->{$geneId}->{$isoformId}->{"chr"} = $cols[0];
                $trsptHref->{$geneId}->{$isoformId}->{"strand"} = $cols[6];
        }elsif($cols[2] eq "exon"){
                $trsptHref->{$geneId}->{$isoformId}->{"exonSeries"} .= $cols[3] . ".." . $cols[4] . ",";
        }
}
close FF;

# 依次将AS的坐标读入到hash中
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
my (%as, @asCatalogFile, $asCatalogFile, @asType, $asType, $i, @gg);
my ($asLine, @valueField, @titleField, $fieldId, $asId);
my ($longAltIsoformIdList, $shortAltIsoformIdList);

open WW, ">$asMappingTrsptTsv";
print WW join("\t", "chr", "strand", "asId", "asType", "1stExonEnd", "1stExonStart_0base", "2ndExonEnd", "2ndExonStart_0base", "downstreamEE", "downstreamES", "exonEnd", "exonStart_0base", "flankingEE", "flankingES", "longExonEnd", "longExonStart_0base", "riExonEnd", "riExonStart_0base", "shortEE", "shortES", "upstreamEE", "upstreamES", "residentType", "asExonSeries", "trsptId", "trsptExonSeries") . "\n";
@asCatalogFile = split(/,/, $ASCatalogList);
@asType = split(/,/, $ASTypeList);
for($i=0; $i<=$#asCatalogFile; $i++){
	$asCatalogFile = $asCatalogFile[$i];
	$asType = $asType[$i];

	# 打开一个AS文件
	open FF, "<$asCatalogFile";

	# 获得字段名
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
		
		# 脱去geneId上的双引号
		@gg = ();
		@gg = split(/\"/, ${$as{$asId}}{"GeneID"});
		${$as{$asId}}{"GeneID"} = $gg[1];

		${$as{$asId}}{"1stExonEnd"} = -1 if(not(exists(${$as{$asId}}{"1stExonEnd"})));
		${$as{$asId}}{"1stExonStart_0base"} = -1 if(not(exists(${$as{$asId}}{"1stExonStart_0base"})));
		${$as{$asId}}{"2ndExonEnd"} = -1 if(not(exists(${$as{$asId}}{"2ndExonEnd"})));
		${$as{$asId}}{"2ndExonStart_0base"} = -1 if(not(exists(${$as{$asId}}{"2ndExonStart_0base"})));
		${$as{$asId}}{"downstreamEE"} = -1 if(not(exists(${$as{$asId}}{"downstreamEE"})));
		${$as{$asId}}{"downstreamES"} = -1 if(not(exists(${$as{$asId}}{"downstreamES"})));
		${$as{$asId}}{"exonEnd"} = -1 if(not(exists(${$as{$asId}}{"exonEnd"})));
		${$as{$asId}}{"exonStart_0base"} = -1 if(not(exists(${$as{$asId}}{"exonStart_0base"})));
		${$as{$asId}}{"flankingEE"} = -1 if(not(exists(${$as{$asId}}{"flankingEE"})));
		${$as{$asId}}{"flankingES"} = -1 if(not(exists(${$as{$asId}}{"flankingES"})));
		${$as{$asId}}{"longExonEnd"} = -1 if(not(exists(${$as{$asId}}{"longExonEnd"})));
		${$as{$asId}}{"longExonStart_0base"} = -1 if(not(exists(${$as{$asId}}{"longExonStart_0base"})));
		${$as{$asId}}{"riExonEnd"} = -1 if(not(exists(${$as{$asId}}{"riExonEnd"})));
		${$as{$asId}}{"riExonStart_0base"} = -1 if(not(exists(${$as{$asId}}{"riExonStart_0base"})));
		${$as{$asId}}{"shortEE"} = -1 if(not(exists(${$as{$asId}}{"shortEE"})));
		${$as{$asId}}{"shortES"} = -1 if(not(exists(${$as{$asId}}{"shortES"})));
		${$as{$asId}}{"upstreamEE"} = -1 if(not(exists(${$as{$asId}}{"upstreamEE"})));
		${$as{$asId}}{"upstreamES"} = -1 if(not(exists(${$as{$asId}}{"upstreamES"})));

########################################################################################
		# AS的longAltExonSeries由 long AltExon 构成
		${$as{$asId}}{"longAltExonSeries"} = &generateLongAltExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});

		# 获得包含AS的long AltExon的isoform列表
		$longAltIsoformIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"longAltExonSeries"}, $trsptHref);
		@isoformId = ();
		@isoformId = split(/,/, $longAltIsoformIdList);

		# 逐个输出和该AS相应的isoform，及其映射关系
		foreach $isoformId(@isoformId){
			$residentType = "longAlt";
			print WW join("\t", ${$as{$asId}}{"chr"}, ${$as{$asId}}{"strand"}, $asId, $asType, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"}, $residentType, ${$as{$asId}}{"longAltExonSeries"}, $isoformId, $trsptHref->{${$as{$asId}}{"GeneID"}}->{$isoformId}->{"exonSeries"}) . "\n";
		}
	
###########################################################################################
		# AS的shortAltExonSeries由 short AltExon 构成
		${$as{$asId}}{"shortAltExonSeries"} = &generateShortAltExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});
	
		# 获得包含AS的short AltExon的isoform列表
		$shortAltIsoformIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"shortAltExonSeries"}, $trsptHref);
		@isoformId = ();
		@isoformId = split(/,/, $shortAltIsoformIdList);
		foreach $isoformId(@isoformId){
			$residentType ="shortAlt";;
			print WW join("\t", ${$as{$asId}}{"chr"}, ${$as{$asId}}{"strand"}, $asId, $asType, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"}, $residentType, ${$as{$asId}}{"shortAltExonSeries"}, $isoformId, $trsptHref->{${$as{$asId}}{"GeneID"}}->{$isoformId}->{"exonSeries"}) . "\n";
		}
	}
	close FF;
}
close WW;



# 在gtf中获得isoform的编号
sub getIsoformIdInAttrs{
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

# 在gtf中获得gene的编号
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


sub getIsoformIdList{
        my ($geneId, $asExonSeries, $geneTrsptExonSeriesHash) = @_;
        my (@isoformId, $isoformId, $isoformIdList);

        return "" if(not exists($geneTrsptExonSeriesHash->{$geneId}));

        @isoformId = keys(%{$geneTrsptExonSeriesHash->{$geneId}});
        foreach $isoformId(@isoformId){
                if(index($geneTrsptExonSeriesHash->{$geneId}->{$isoformId}->{"exonSeries"}, $asExonSeries)>=0){
                        $isoformIdList.=$isoformId . ",";
                }
        }
        return substr($isoformIdList, 0, length($isoformIdList) - 1);
}



# 默认includion(A5SS,A3SS,RI,SE）为第1链，默认MXE为包含小坐标的备选外显子
sub generateLongAltExonSeries{
	my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		# 默认包含长外显子为第1链
		$longExonStart_0base = $longExonStart_0base + 1;
		$flankingES = $flankingES + 1;		
		return $longExonStart_0base . ".." . $longExonEnd . "," . $flankingES . ".." . $flankingEE;
	}elsif($asType eq "A3SS"){
		# 默认包含长外显子为第1链
		$longExonStart_0base = $longExonStart_0base + 1;
		$flankingES = $flankingES + 1;
		return $flankingES . ".." . $flankingEE . "," . $longExonStart_0base . ".." . $longExonEnd;	
	}elsif($asType eq "RI"){
		# 默认内含子保留为第1链
		$riExonStart_0base = $riExonStart_0base + 1;
		return $riExonStart_0base . ".." . $riExonEnd;	
	}elsif($asType eq "SE"){
		# 默认包含中间外显子为第1链
		$exonStart_0base = $exonStart_0base + 1;
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $exonStart_0base . ".." . $exonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $exonStart_0base . ".." . $exonEnd . "," . $upstreamES . ".." .  $upstreamEE;
		}
	}elsif($asType eq "MXE"){
		# 默认包含小坐标备用外显子为第1链
		$firstExonStart_0base = $firstExonStart_0base + 1;
		$upstreamES = $upstreamES + 1 ;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}
	
}


# 默认包含skiping(A5SS,A3SS,RI,SE)为第2链，MXE则默认包含大坐标备选外显子为第2链
sub generateShortAltExonSeries{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		# 默认包含小外显子为第2链
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $shortES . ".." . $shortEE . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		# 默认包含小外显子为第2链
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;	
		return $flankingES . ".." . $flankingEE . "," . $shortES . ".." . $shortEE;	
	}elsif($asType eq "RI"){
		# 默认内含子不保留为第2链
		$upstreamES = $upstreamES + 1;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." . $upstreamEE;
		}
	}elsif($asType eq "SE"){
		# 默认不包含中间外显子为第2链
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE;
		}
	}elsif($asType eq "MXE"){
		# 默认包含大坐标备用外显子为第2链
		$secondExonStart_0base = $secondExonStart_0base + 1;
		$upstreamES = $upstreamES + 1 ;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $secondExonStart_0base . ".." . $secondExonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $secondExonStart_0base . ".." . $secondExonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}		
}



