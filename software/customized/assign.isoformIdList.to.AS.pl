#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--ASCatalogList\\\n" .
		"--ASTypeList \\\n" .
		"--ensemblGtf \\\n" . 
		"--improvementGtf \\\n" .
		"--outputTagFile \n";
        exit;
}

my ($ASCatalogList, $ASTypeList, $ensemblGtf, $improvementGtf, $outputTagFile);
GetOptions(
	'ASTypeList=s'=>\$ASTypeList,
	'ASCatalogList=s'=>\$ASCatalogList,
	'ensemblGtf=s'=>\$ensemblGtf,
	'improvementGtf=s'=>\$improvementGtf,
	'outputTagFile=s'=>\$outputTagFile,
);

my ($featureLine, @cols, $geneId, @isoformId, $isoformId, $isoformSymbol);
my (%ensembl, $ensemblHref);
$ensemblHref = \%ensembl;

open WW, ">$outputTagFile";
print WW join("\t", "ASID", "chr", "strand", "longAltExonSeries", "longAltEnsemblTrsptIdList", "longAltImprovedTrsptIdList", "shortAltExonSeries", "shortAltEnsemblTrsptIdList", "shortAltImprovedTrsptIdList") . "\n";

# 将ensembl中的trspt对应的外显子串读入到hash中
open FF, "<$ensemblGtf";
while($featureLine = <FF>){
        chomp($featureLine);
        next if($featureLine=~/#/);
        @cols = split(/\t/, $featureLine);

        $isoformId = &getIsoformIdInAttrs($cols[8]);
        $geneId = &getGeneIdInAttrs($cols[8]);
        if($cols[2] eq "transcript"){
                $ensemblHref->{$geneId}->{$isoformId}->{"chr"} = $cols[0];
                $ensemblHref->{$geneId}->{$isoformId}->{"strand"} = $cols[6];
        }elsif($cols[2] eq "exon"){
                $ensemblHref->{$geneId}->{$isoformId}->{"exonSeries"} .= $cols[3] . ".." . $cols[4] . ",";
        }
}
close FF;


# 将improvement对应外显子串读到hash中
my (%improvement, $improvementHref);
$improvementHref = \%improvement;
open FF, "<$improvementGtf";
while($featureLine = <FF>){
        chomp($featureLine);
        next if($featureLine=~/#/);
        @cols = split(/\t/, $featureLine);

        $isoformId = &getIsoformIdInAttrs($cols[8]);
        $geneId = &getGeneIdInAttrs($cols[8]);
        
	if($cols[2] eq "transcript"){
                $improvementHref->{$geneId}->{$isoformId}->{"chr"} = $cols[0];
                $improvementHref->{$geneId}->{$isoformId}->{"strand"} = $cols[6];
        }elsif($cols[2] eq "exon"){
		$improvementHref->{$geneId}->{$isoformId}->{"exonSeries"} .= $cols[3] . ".." . $cols[4] . ",";
        }
}
close FF;

# 从逐个打开AS文件，将每个AS读入到hash中进行检索
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
my (%as, @asCatalogFile, $asCatalogFile, @asType, $asType, $i, @gg);
my ($asLine, @valueField, @titleField, $fieldId, $asId);
my ($longAltEnsemblTrsptIdList, $longAltImprovedTrsptIdList, $shortAltEnsemblTrsptIdList, $shortAltImprovedTrsptIdList);
@asCatalogFile = split(/,/, $ASCatalogList);
@asType = split(/,/, $ASTypeList);
for($i=0; $i<=$#asCatalogFile; $i++){
	$asCatalogFile = $asCatalogFile[$i];
	$asType = $asType[$i];

	# 打开一个AS文件
	open FF, "<$asCatalogFile";
	# 获得列名进标题数组
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
		
		# 获得包含longAlt的trsptId列表。对于mxe而言longAlt指的是包含1st备用外显子
		# 获得AS的longAlt外显子串
		${$as{$asId}}{"longAltExonSeries"} = &generateLongAltAsExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});
		# 获得AS的longAlt外显子串对应的ensembl和improved中的trsptIdList
		$longAltEnsemblTrsptIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"longAltExonSeries"}, $ensemblHref);
		$longAltImprovedTrsptIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"longAltExonSeries"}, $improvementHref);
		# improvement中包含了ensembl，所以需要如下处理，才能得到真正improvedTrsptIdList：
		# 在longAltImprovementIsoformIdList中减去longAltEnsembIsoformIdList
		$longAltImprovedTrsptIdList = &removeEnsemblFromImprovement($longAltEnsemblTrsptIdList, $longAltImprovedTrsptIdList);


#############

		# 获得包含shortAltExon的外显子串。对于mxe而言shortAlt指的是包含2nd备用外显子
		${$as{$asId}}{"shortAltExonSeries"} = &generateShortAltAsExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});
		# 获得AS的shortAlt外显子串对应的ensembl和improved中的trsptIdList
		$shortAltEnsemblTrsptIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"shortAltExonSeries"}, $ensemblHref);
		$shortAltImprovedTrsptIdList = &getIsoformIdList(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"shortAltExonSeries"}, $improvementHref);
		# improvement中包含了ensembl，所以需要如下处理，才能得到真正improvedTrsptIdList：
		# 在longAltImprovementIsoformIdList中减去longAltEnsembIsoformIdList
		$shortAltImprovedTrsptIdList = &removeEnsemblFromImprovement($shortAltEnsemblTrsptIdList, $shortAltImprovedTrsptIdList);

		print WW join("\t", $asId, ${$as{$asId}}{"chr"}, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"longAltExonSeries"}, $longAltEnsemblTrsptIdList, $longAltImprovedTrsptIdList, ${$as{$asId}}{"shortAltExonSeries"}, $shortAltEnsemblTrsptIdList, $shortAltImprovedTrsptIdList) . "\n";
	}
	close FF;
}
close WW;


# 从improvementIsoformIdList中减去ensemblIsoformIList
sub removeEnsemblFromImprovement{
	my ($ensemblIdList, $improvementIdList) =@_;
	my ($ensemblId, @ensemblId, $improvementId, @improvementId, %ensemblId, $returnId);
	
	# ensemblIdList 为 “-”，
	# 	那么直接返回improvementIdList。如果improvementIdList为"-"，那么返回的正确；如果improvementIdList不为"-"，那么也是正确的。
	if($ensemblIdList eq "-"){
		return $improvementIdList;
	}
	# 如果improvementIdList = ensemblIdList，那么返回"-"
	if($ensemblIdList eq $improvementIdList){
		return "-";
	}
	
	# improvementIdList != ensemblIdList 或者 improvementIdList= ensemblIdList但是次序不一致时，
	$returnId = "";
	@ensemblId = split(/,/, $ensemblIdList);
	@improvementId = split(/,/, $improvementIdList);
	foreach $ensemblId(@ensemblId){
		$ensemblId{$ensemblId}=1;
	}
	foreach $improvementId(@improvementId){
		$returnId.= $improvementId . "," if(not exists($ensemblId{$improvementId}));
	}
	if($returnId eq ""){# improvementIdList != ensemblIdList，只是id次序不一致
		return "-";
	}else{		    # improvementIdList != ensemblIdList
		return substr($returnId, 0, length($returnId) - 1 );
	}
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


# 检查AS的某条链是否出现在指定的gtf中
sub getIsoformIdList{
	my ($geneId, $asExonSeries, $gtf) = @_;
	my (@isoformId, $isoformId, $isoformIdList);
	$isoformIdList = "";
	return "-" if(not exists($gtf->{$geneId}));

	@isoformId = keys(%{$gtf->{$geneId}});
	foreach $isoformId(@isoformId){
		if(index($gtf->{$geneId}->{$isoformId}->{"exonSeries"}, $asExonSeries)>=0){
			$isoformIdList.=$isoformId . ",";
		}
	}
	if($isoformIdList ne ""){
		return substr($isoformIdList, 0, length($isoformIdList) - 1);
	}else{
		return "-";
	}
}


# 对于mxe类型，longAlt指的是包含1st备用外显子
# 正链时从小到大，负链时从大到小
sub generateLongAltAsExonSeries{
	my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$flankingES = $flankingES + 1;		
		return $longExonStart_0base . ".." . $longExonEnd . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$flankingES = $flankingES + 1;
		return $flankingES . ".." . $flankingEE . "," . $longExonStart_0base . ".." . $longExonEnd;	
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
			return $downstreamES . ".." .  $downstreamEE . "," . $exonStart_0base . ".." . $exonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}elsif($asType eq "MXE"){
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

# 对于mxe类型AS，shortAlt为包含2nd备用外显子
sub generateShortAltAsExonSeries{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){

		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $shortES . ".." . $shortEE . "," . $flankingES . ".." . $flankingEE;		

	}elsif($asType eq "A3SS"){

		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;	
		return $flankingES . ".." . $flankingEE . "," . $shortES . ".." . $shortEE;

	}elsif($asType eq "RI"){
		$upstreamES = $upstreamES + 1;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE;
		}

	}elsif($asType eq "SE"){
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." . $upstreamEE;
		}
	}elsif($asType eq "MXE"){
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

