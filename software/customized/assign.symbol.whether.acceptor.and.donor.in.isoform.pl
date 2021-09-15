#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--ASCatalogList\\\n" .
		"--ASTypeList \\\n" .
		"--improvementGtf \\\n" .
		"--outputTagFile \n";
        exit;
}

my ($ASCatalogList, $ASTypeList, $ensemblGtf, $improvementGtf, $outputTagFile);
GetOptions(
	'ASTypeList=s'=>\$ASTypeList,
	'ASCatalogList=s'=>\$ASCatalogList,
	'improvementGtf=s'=>\$improvementGtf,
	'outputTagFile=s'=>\$outputTagFile,
);

my ($featureLine, @cols, $geneId, @isoformId, $isoformId, $isoformSymbol);
my (%ensembl, $ensemblHref);
$ensemblHref = \%ensembl;

open WW, ">$outputTagFile";
print WW join("\t", "ASID", "coordinateInGtf") . "\n";
# 将improvement对应cDNA坐标读入到hash中
my (%improvement, $chr, $strand, $start, $stop);
open FF, "<$improvementGtf";
while($featureLine = <FF>){
        chomp($featureLine);
        next if($featureLine=~/#/);
        @cols = split(/\t/, $featureLine);
	$chr = $cols[0];
	$start = $cols[3];
	$stop = $cols[4];
	$strand = $cols[6];
	$improvement{$chr . "_" . $strand . "_" . $start} = 1;
	$improvement{$chr . "_" . $strand . "_" . $stop} = 1;
#	print $chr . "_" . $strand . "_" . $start . "\n";
#	print $chr . "_" . $strand . "_" . $stop . "\n";
#	<STDIN>;
}
close FF;

#open ASW, ">as.coordinate.txt";
#print "Finish load improvement\n";
# 依次将AS的坐标读入到hash中
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
my (%as, @asCatalogFile, $asCatalogFile, @asType, $asType, $i, @gg);
my ($asLine, @valueField, @titleField, $fieldId, $asId);
my ($firstEnsembIsoformIdList, $firstImprovementIsoformIdList, $secondEnsembIsoformIdList, $secondImprovementIsoformIdList);
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
		
		# 脱去geneId上的双引号
		@gg = ();
		@gg = split(/\"/, ${$as{$asId}}{"GeneID"});
		${$as{$asId}}{"GeneID"} = $gg[1];
	
		my $flag = "Y";
		my $coordinate = 0;
		my $position = "";
		if(exists(${$as{$asId}}{"1stExonEnd"})){
			$coordinate = ${$as{$asId}}{"1stExonEnd"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"1stExonStart_0base"})){
			$coordinate = ${$as{$asId}}{"1stExonStart_0base"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"2ndExonEnd"})){
			$coordinate = ${$as{$asId}}{"2ndExonEnd"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"2ndExonStart_0base"})){
			$coordinate = ${$as{$asId}}{"2ndExonStart_0base"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"downstreamEE"})){
			$coordinate = ${$as{$asId}}{"downstreamEE"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"downstreamES"})){
			$coordinate = ${$as{$asId}}{"downstreamES"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"exonEnd"})){
			$coordinate = ${$as{$asId}}{"exonEnd"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"exonStart_0base"})){
			$coordinate = ${$as{$asId}}{"exonStart_0base"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"flankingEE"})){
			$coordinate = ${$as{$asId}}{"flankingEE"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"flankingES"})){
			$coordinate = ${$as{$asId}}{"flankingES"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"longExonEnd"})){
			$coordinate = ${$as{$asId}}{"longExonEnd"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"longExonStart_0base"})){
			$coordinate = ${$as{$asId}}{"longExonStart_0base"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"riExonEnd"})){
			$coordinate = ${$as{$asId}}{"riExonEnd"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"riExonStart_0base"})){
			$coordinate = ${$as{$asId}}{"riExonStart_0base"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"shortEE"})){
			$coordinate = ${$as{$asId}}{"shortEE"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"shortES"})){
			$coordinate = ${$as{$asId}}{"shortES"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"upstreamEE"})){
			$coordinate = ${$as{$asId}}{"upstreamEE"} + 0;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		if(exists(${$as{$asId}}{"upstreamES"})){
			$coordinate = ${$as{$asId}}{"upstreamES"} + 1;
			$position = ${$as{$asId}}{"chr"} . "_" . ${$as{$asId}}{"strand"} . "_" . $coordinate;
			#print ASW $position . "\n";
			if(not exists($improvement{$position})){
				$flag = "N";
				goto OUTPUT;
			}
		}

		
OUTPUT:
		print WW join("\t", $asId, $flag) . "\n";

	}
	close FF;
}
close WW;
#close ASW;

# 从improvementIsoformIdList中减去ensemblIsoformIList
sub removeEnsemblFromImprovement{
	my ($ensemblIdList, $improvementIdList) =@_;
	my ($ensemblId, @ensemblId, $improvementId, @improvementId, %ensemblId, $returnId);

	$returnId = "";
	@ensemblId = split(/,/, $ensemblIdList);
	@improvementId = split(/,/, $improvementIdList);
	foreach $ensemblId(@ensemblId){
		$ensemblId{$ensemblId}=1;
	}
	foreach $improvementId(@improvementId){
		$returnId.= $improvementId . "," if(not exists($ensemblId{$improvementId}));
	}
	return $returnId;
}


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


# 检查AS的两条链是否都出现在指定的gtf中
sub getIsoformIdList{
	my ($geneId, $asExonSeries, $gtf) = @_;
	my (@isoformId, $isoformId, $isoformIdList);

	return "" if(not exists($gtf->{$geneId}));

	@isoformId = keys(%{$gtf->{$geneId}});
	foreach $isoformId(@isoformId){
		if(index($gtf->{$geneId}->{$isoformId}->{"exonSeries"}, $asExonSeries)>=0){
			$isoformIdList.=$isoformId . ",";
		}
	}
	return substr($isoformIdList, 0, length($isoformIdList) - 1);
}

# 根据AS的坐标信息生成AS的外显子坐标串
# 正链时从小到大，负链时从大到小
sub generateFirstAltAsExonSeries{
	my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $longExonStart_0base . ".." . $longExonEnd . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;
		return $flankingES . ".." . $flankingEE . "," . $longExonStart_0base . ".." . $longExonEnd;	
	}elsif($asType eq "RI"){
		$riExonStart_0base = $riExonStart_0base + 1;
		$upstreamES = $upstreamES + 1;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE;
		}
	}elsif($asType eq "SE"){
		$exonStart_0base = $exonStart_0base + 1;
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $exonStart_0base . ".." . $exonEnd;
		}
	}elsif($asType eq "MXE"){
		$firstExonStart_0base = $firstExonStart_0base + 1;
		$secondExonStart_0base = $secondExonStart_0base + 1;
		$upstreamES = $upstreamES + 1 ;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}
	
}

sub generateSecondAltAsExonSeries{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $shortES . ".." . $shortEE . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;	
		return $flankingES . ".." . $flankingEE . "," . $shortES . ".." . $shortEE ;	
	}elsif($asType eq "RI"){
		$riExonStart_0base = $riExonStart_0base + 1;
		return $riExonStart_0base . ".." . $riExonEnd;
	}elsif($asType eq "SE"){
		$exonStart_0base = $exonStart_0base + 1;
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $exonStart_0base . ".." . $exonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE . "," . $exonStart_0base . ".." . $exonEnd;
		}
	}elsif($asType eq "MXE"){
		$firstExonStart_0base = $firstExonStart_0base + 1;
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



