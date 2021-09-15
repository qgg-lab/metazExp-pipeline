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

# 将ensembl对应的cDNA坐标读入到hash
my ($featureLine, @cols, $geneId, @isoformId, $isoformId, $isoformSymbol);
my (%ensembl, $ensemblHref);
$ensemblHref = \%ensembl;
open WW, ">$outputTagFile";
print WW join("\t", "ASID", "ensemblGtf", "improvedGtf", "RNAseqMapping") . "\n";
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
print $ensemblHref->{"AT1G01010 "}->{"AT1G01010.1"}->{"exonSeries"};

# 将improvement对应cDNA坐标读入到hash中
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

#print "Finish load improvement\n";
# 依次将AS的坐标读入到hash中
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
my (%as, @asCatalogFile, $asCatalogFile, @asType, $asType, $i, @gg);
my ($asLine, @valueField, @titleField, $fieldId, $asId);
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
		
		# 获得as的longAlt对应的外显子串，mxe类型的longAlt包含1st备用外显子
		${$as{$asId}}{"longAltExonSeries"} = &generateLongAltAsExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});
		# 得as的shortAlt对应的外显子串, mxe类型的shortAlt包含2nd备用外显子
		${$as{$asId}}{"shortAltExonSeries"} = &generateShortAltAsExonSeries($asType, ${$as{$asId}}{"strand"}, ${$as{$asId}}{"1stExonEnd"}, ${$as{$asId}}{"1stExonStart_0base"}, ${$as{$asId}}{"2ndExonEnd"}, ${$as{$asId}}{"2ndExonStart_0base"}, ${$as{$asId}}{"downstreamEE"}, ${$as{$asId}}{"downstreamES"}, ${$as{$asId}}{"exonEnd"}, ${$as{$asId}}{"exonStart_0base"}, ${$as{$asId}}{"flankingEE"}, ${$as{$asId}}{"flankingES"}, ${$as{$asId}}{"longExonEnd"}, ${$as{$asId}}{"longExonStart_0base"}, ${$as{$asId}}{"riExonEnd"}, ${$as{$asId}}{"riExonStart_0base"}, ${$as{$asId}}{"shortEE"}, ${$as{$asId}}{"shortES"}, ${$as{$asId}}{"upstreamEE"}, ${$as{$asId}}{"upstreamES"});

		# 检查该AS的两个外显子串是否同时在ensembl注释的isoform坐标中出现
		if(&checkASInGtf(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"longAltExonSeries"}, ${$as{$asId}}{"shortAltExonSeries"}, $ensemblHref)==1){
			${$as{$asId}}{"ensembl"} = 1;
		}else{
			${$as{$asId}}{"ensembl"} = 0;
		}

		# 检查该AS的两个外显子串是否同时在improvement注释的isoform坐标中出现
		if(${$as{$asId}}{"ensembl"} ==0 and &checkASInGtf(${$as{$asId}}{"GeneID"}, ${$as{$asId}}{"longAltExonSeries"}, ${$as{$asId}}{"shortAltExonSeries"}, $improvementHref)==1){
			${$as{$asId}}{"improvement"} = 1;
		}else{
			${$as{$asId}}{"improvement"} = 0;
		}
	
		# 不在ensembl和improvement中出现，则将novel标注为1，否则标注为0
		if(${$as{$asId}}{"ensembl"} == 0 and ${$as{$asId}}{"improvement"} == 0){
			${$as{$asId}}{"novel"} = 1;
		}else{
			${$as{$asId}}{"novel"} = 0;
		}
		print WW join("\t", $asId, ${$as{$asId}}{"ensembl"}, ${$as{$asId}}{"improvement"}, ${$as{$asId}}{"novel"}) . "\n";
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


# 检查AS的两条链是否都出现在指定的gtf中
sub checkASInGtf{
	my ($geneId, $firstExonSeries, $secondExonSeries, $gtf) = @_;
	my (@isoformId, $isoformId, $firstFlag, $secondFlag);

	$firstFlag = 0;
	$secondFlag = 0;
	if(not exists($gtf->{$geneId})){
		goto RRTT;
	}
	@isoformId = keys(%{$gtf->{$geneId}});
	# 首先检查第一个外显子串
	foreach $isoformId(@isoformId){
		if(index($gtf->{$geneId}->{$isoformId}->{"exonSeries"}, $firstExonSeries)>=0){
			$firstFlag = 1;
			last;
		}
	}
	# 检查第二个外显子串
	foreach $isoformId(@isoformId){
		if(index($gtf->{$geneId}->{$isoformId}->{"exonSeries"}, $secondExonSeries)>=0){
			$secondFlag = 1;
			last;
		}
	}
RRTT:
	return $firstFlag * $secondFlag;
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


sub getAsStartAndEnd{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $start, $end)=@_;
        if($asType eq "A5SS"){
                if($strand eq "+"){
                        $$start = $longExonStart_0base+1;
                        $$end = $flankingEE;
                }else{
                        $$start = $flankingES+1;
                        $$end = $longExonEnd
                }
        }elsif($asType eq "A3SS"){
                if($strand eq "+"){
                        $$start = $flankingES+1;
                        $$end = $longExonEnd;
                }else{
                        $$start = $longExonStart_0base + 1;
                        $$end = $flankingEE;
                }
        }
        if($asType eq "RI"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
        if($asType eq "SE"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
        if($asType eq "MXE"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
}
