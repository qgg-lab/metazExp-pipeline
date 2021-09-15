#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf \\\n" .
                "--asCatalogFileList \\\n" .
		"--asTypeList \\\n" .
		"--pfamInOrigCdna \\\n" .
		"--goUniqListInOrigCdna \\\n" .
		"--KIdToEnzymeId \\\n" .
		"--geneToKId \\\n" .
		"--geneToPathwayId \\\n" .
		"--tissueSpecificRltTsv \\\n" .
		"--speciesAbbr ATHA \\\n" .
		"--symbolAnno \\\n" .
		"--geneMysqlTsv \n";
	exit;
}

my (
$gtf,
$asCatalogFileList, 
$asTypeList,
$pfamInOrigCdna, 
$goUniqListInOrigCdna, 
$geneToKId,
$KIdToEnzymeId,
$geneToPathwayId,
$tissueSpecificRltTsv,
$speciesAbbr,
$symbolAnno,
$geneMysqlTsv);

GetOptions(
        'gtf=s'=>\$gtf,
        'asCatalogFileList=s'=>\$asCatalogFileList,
        'asTypeList=s'=>\$asTypeList,
        'pfamInOrigCdna=s'=>\$pfamInOrigCdna,
	'goUniqListInOrigCdna=s'=>\$goUniqListInOrigCdna,
	'geneToKId=s'=>\$geneToKId,
	'KIdToEnzymeId=s'=>\$KIdToEnzymeId,
	'geneToPathwayId=s'=>\$geneToPathwayId,
	'tissueSpecificRltTsv=s'=>\$tissueSpecificRltTsv,
	'speciesAbbr=s'=>\$speciesAbbr,
	'symbolAnno=s'=>\$symbolAnno,
	'geneMysqlTsv=s'=>\$geneMysqlTsv,
);

my (%trsptToGene, %gene, $geneHref);
$geneHref=\%gene;

# 读取gtf中source
my (@line, $trsptId, $geneId, $geneSymbol, $trsptSymbol, $line, @field);
open FF, "<$gtf";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");

	($geneId, $trsptId, $geneSymbol, $trsptSymbol) = ("", "", "NA", "NA");
	&getGeneIdAndTrspt($field[8], \$geneId, \$trsptId, \$geneSymbol, \$trsptSymbol);

	# 登记trspt到gene之间映射
	$trsptToGene{$trsptId} = $geneId;

	# 登记gene的geneSymbol信息
	if(not exists($geneHref->{$geneId}->{"geneSymbol"})){
		$geneHref->{$geneId}->{"geneSymbol"} = $geneSymbol;
		$geneHref->{$geneId}->{"chr"} = $field[0];
		$geneHref->{$geneId}->{"strand"} = $field[6];
		$geneHref->{$geneId}->{"start"} = $field[3];
		$geneHref->{$geneId}->{"stop"} = $field[4];
	}else{
		if($field[3]< $geneHref->{$geneId}->{"start"}){
			$geneHref->{$geneId}->{"start"} = $field[3];
		}
		if($field[4] > $geneHref->{$geneId}->{"stop"}){
			$geneHref->{$geneId}->{"stop"} = $field[4]
		}
	}
	# 转录本列表登记到gene中
	if(not exists($geneHref->{$geneId}->{"trsptList"})){
		$geneHref->{$geneId}->{"trsptList"} = $trsptId;
	}else{
		$geneHref->{$geneId}->{"trsptList"}.= "," . $trsptId;
	}

}

# 读取gene的Symbol注释
my (@symbolField);
open FF, "<$symbolAnno";
# geneId          annoSymbol trId    trSymbol trFullName spId    spSymbol  spFullName
# LOC108854717    NA         NA      NA       NA         NA      NA        NA
<FF>;
while($line=<FF>){
	chomp($line);
	@symbolField = ();
	@symbolField = split(/\t/, $line);
	$geneId = shift(@symbolField);
	$geneHref->{$geneId}->{"annoSymbol"} = shift(@symbolField);
	$geneHref->{$geneId}->{"trId"} = shift(@symbolField);
	$geneHref->{$geneId}->{"trSymbol"} = shift(@symbolField);
	$geneHref->{$geneId}->{"trFullName"} = shift(@symbolField);
	$geneHref->{$geneId}->{"spId"} = shift(@symbolField);
	$geneHref->{$geneId}->{"spSymbol"} = shift(@symbolField);
	$geneHref->{$geneId}->{"spFullName"} = shift(@symbolField);
}
close FF;


# 将可变剪接登记到gene中
# ASID    GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE
# ATHASE0000010226        "AT1G01010"     "NAC001"        1       +       4485    4605    3995    4276    4705    5095
my (@asFile, $asFile, @asType, $asType, $i, @field);
@asFile = split(/,/, $asCatalogFileList);
@asType = split(/,/, $asTypeList);
for($i=0; $i<=$#asFile; $i++){
	$asFile = $asFile[$i];
	$asType = $asType[$i];
	open FF, "<$asFile";
	<FF>;
	while($line=<FF>){
		@field = split(/\t/, $line);
		if($field[1]=~/"(.*)"/){
			$geneId = $1;
			if(not exists($geneHref->{$geneId}->{$asType . "List"})){
				$geneHref->{$geneId}->{$asType . "List"} = $field[0];
			}else{
				$geneHref->{$geneId}->{$asType . "List"} .= "," . $field[0];
			}
		}
	}
	close FF;
}

my ($pfamPosList, @pfamPos, $pfamPos, $pfamId);
# 将pfamInCdna读入
open FF, "<$pfamInOrigCdna";
# AT3G11450.1     PF00226[547,786]|PF00249[2020,2160]
# AT5G58690.3     PF09279[226,450]|PF00388[487,915]|PF00387[1171,1431]|PF00168[1495,1803]
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamPosList) = split(/\t/, $line);
	$geneId = $trsptToGene{$trsptId};
	@pfamPos = ();
	@pfamPos = split(/\|/, $pfamPosList);
	foreach $pfamPos(@pfamPos){
		if($pfamPos=~/(PF\d+)\[\d+,\d+\]/){
			$pfamId = $1;
			if(not exists($geneHref->{$geneId}->{"origPfamList"})){
				$geneHref->{$geneId}->{"origPfamList"} = $pfamId;
			}elsif(index($geneHref->{$geneId}->{"origPfamList"}, $pfamId )<0){
				$geneHref->{$geneId}->{"origPfamList"} .= "," . $pfamId;
			}
		}
	}
}
close FF;

# 将as作用后的trspt上的pfam读入(这部分删除不需要了！)
my ($trsptIdAsId, $asId);
=beg
open FF, "<$pfamInAlteredCdna";
# AT2G01690.1XXXATHASE0000028750  PF12755[499,678]|PF11916[1471,2010]
while($line=<FF>){
	chomp($line);
	($trsptIdAsId, $pfamPosList) = split(/\t/, $line);
	($trsptId, $asId) = split(/XXX/, $trsptIdAsId);
	$geneId = $trsptToGene{$trsptId};
	@pfamPos = ();
	@pfamPos = split(/\|/, $pfamPosList);
	foreach $pfamPos(@pfamPos){
		if($pfamPos=~/(PF\d+)\[\d+,\d+\]/){
			$pfamId = $1;
			if(not exists($geneHref->{$geneId}->{"alteredPfamList"})){
				$geneHref->{$geneId}->{"alteredPfamList"} = $pfamId;
			}elsif(index($geneHref->{$geneId}->{"alteredPfamList"}, $pfamId)<0){
				$geneHref->{$geneId}->{"alteredPfamList"} .= "," . $pfamId;
			}
		}
	}
}
close FF;
=cut;

# 将go注释读入
my ($goList, @go, $go);
open FF, "<$goUniqListInOrigCdna";
# AT3G08840.13    GO:0005524|GO:0008716|GO:0046872
# AT1G49720.2     GO:0003700|GO:0006355 
while($line=<FF>){
	chomp($line);
	($trsptId, $goList) = split(/\t/, $line);
	$geneId = $trsptToGene{$trsptId};
	@go = ();
	@go = split(/\|/, $goList);
	foreach $go(@go){
		if(not exists($geneHref->{$geneId}->{"origGoList"})){
			$geneHref->{$geneId}->{"origGoList"} = $go;
		}elsif(index($geneHref->{$geneId}->{"origGoList"}, $go)<0){
			$geneHref->{$geneId}->{"origGoList"} .= "," . $go;
		}
	}
}
close FF;

# 将as作用后的go注释读入(这部分删除不需要了！)
my ($goList, @go, $go);
=beg
open FF, "<$goUniqListInAlteredCdna";
# AT3G08840.13    GO:0005524|GO:0008716|GO:0046872
# AT1G49720.2     GO:0003700|GO:0006355 
while($line=<FF>){
	chomp($line);
	($trsptIdAsId, $goList) = split(/\t/, $line);
	($trsptId, $asId) = split(/XXX/, $trsptIdAsId);
	$geneId = $trsptToGene{$trsptId};
	@go = ();
	@go = split(/\|/, $goList);
	foreach $go(@go){
		if(not exists($geneHref->{$geneId}->{"alteredGoList"})){
			$geneHref->{$geneId}->{"alteredGoList"} = $go;
		}elsif(index($geneHref->{$geneId}->{"alteredGoList"}, $go)<0){
			$geneHref->{$geneId}->{"alteredGoList"} .= "," . $go;
		}
	}
}
close FF;
=cut;


# 组织特异表达分析结果读入
# pvalue0.01Rlt   pvalue0.03Rlt   pvalue0.05Rlt   geneId  regulation      tissue/treatment        avgExpTpm       experimentNum   1stTissue/treatment|avgTpm|exptNum|pvalue       2ndTissue/treatment|avgTpm|exptNum|pvalue       3rdTissue/Treatment|avgTpm|exptNum|pvalue       4thTissue/Treatment|avgTpm|exptNum|pvalue       ...
# N       N       N       AT5G13565       High    meristem        77.655  6       ovule|57.702|4|0.057    receptacle|0.000|3|0.994        shoot|3.961|39|1.000    pistil|0.010|10|1.000   embryo|0.000|6|0.999    seed|5.972|48|0.841     silique|7.902|15|0.919  seedling|13.657|238|0.785       root|12.596|257|0.947   stem|15.885|9|0.072     cotyledon|0.021|9|1.000 leaf|2.043|112|0.996    endosperm|16.051|16|0.290       flower|26.797|19|0.564  hypocotyl|0.016|9|1.000 anther|10.275|4|0.832   inflorescence|16.935|18|0.742
my (@fieldName, @fieldValue, $valueList, %specificHash, $pvalueLevel, $tissueNum);
open FF, "<$tissueSpecificRltTsv";
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@fieldValue = split(/\t/, $line);
	next if($fieldValue[0] eq "N" and $fieldValue[1] eq "N" and $fieldValue[2] eq "N");
	if($fieldValue[2] eq "Y"){
		$pvalueLevel = "0.05";
	}
	if($fieldValue[1] eq "Y"){
		$pvalueLevel = "0.03";
	}
	if($fieldValue[0] eq "Y"){
		$pvalueLevel = "0.01";
	}

	$geneId = $fieldValue[3];
	
	$tissueNum = $#fieldValue + 1 - 6;
	if(uc($fieldValue[4]) eq "HIGH"){
		$geneHref->{$geneId}->{"specificHighTissue"} = $fieldValue[5];
		$geneHref->{$geneId}->{"specificHighTissuePvalue"} = $pvalueLevel;
		$geneHref->{$geneId}->{"specificHighTissueNum"} = $tissueNum;
	}elsif(uc($fieldValue[4]) eq "LOW"){
		$geneHref->{$geneId}->{"specificLowTissue"} = $fieldValue[5];
		$geneHref->{$geneId}->{"specificLowTissuePvalue"} = $pvalueLevel;
		$geneHref->{$geneId}->{"specificLowTissueNum"} = $tissueNum;
	}

}
close FF;


# 将基因的ORTHOLOGGROUP编号读入: 物种间用\t分隔，物种内的基因用", "分隔
# Orthogroup      AHAL    AHYP    ALYR    ATHA    BNAP    BRAP    BVUL    CANN    CARA    CMEL    CSAT    CSIN    DCAR
# OG0000000       AHAL_g00262_pep001, AHAL_g00524_pep001,
my (%tmpHash, @nameField, @valueField, $modifiedGeneIdList, @modifiedGeneId, $modifiedGeneId, $orthoGroupId);
=beg
open FF, "<$orthoroupTsv";
$line=<FF>,
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	# orthoGroupId
	$orthoGroupId = $valueField[0];

	# 将此orthGroup的ID读入到hash中
	for($i=0; $i<=$#nameField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	$modifiedGeneIdList = $tmpHash{$speciesAbbr};

	# 根据指定物种名缩写，从hash中取出指定物种下的修饰的geneId编号列表
	# 如果当前orthogroup下没有对应的基因，那么跳到下一行
	next if($modifiedGeneIdList eq "");
	@modifiedGeneId = ();
	@modifiedGeneId = split(/, /, $modifiedGeneIdList);
	foreach $modifiedGeneId(@modifiedGeneId){
		# AHAL_g00262_pep001
		if($modifiedGeneId=~/$speciesAbbr(_)(.*)(_pep.*)/){
			$geneId = $2;
			# 将orthoGroupId登记到总gene信息表中
			$geneHref->{$geneId}->{"orthogroupId"} = $orthoGroupId;
		}
	}	
}
close FF;
=cut;

# 读取koToEnzyme，获得每个koId对应的酶号
my (%koIdToEnzyme, @tt, $KId, $enzymeIdList);
open FF, "<$KIdToEnzymeId";
# K00028  EC:1.1.1.50 1.1.1.357 1.1.1.225
while($line=<FF>){
	chomp($line);
	($KId, $enzymeIdList)= split(/\t/, $line);
	if($enzymeIdList eq "-"){
		$enzymeIdList = "NA";
	}
	$koIdToEnzyme{$KId} = $enzymeIdList;
}
close FF;

# 将gene注释到的KId和对应的酶号添加到hash中
open FF, "<$geneToKId";
while(my $line=<FF>){
	if($line=~/^(.*)\t(.*)\n/){
		$geneId = $1;
		$KId = $2;
		$geneHref->{$geneId}->{"KId"} = $KId;
		if(exists($koIdToEnzyme{$KId})){
			$geneHref->{$geneId}->{"enzymeId"} = $koIdToEnzyme{$KId};
		}
	}
}
close FF;


# 读取geneToPathwayId，让geneId -> pathwayIdList
my ($nameList, @nameFile, $name, $valueList, @valueField, $value, %tmpHash, $pathwayId);
open FF, "<$geneToPathwayId";
# 1stClassId 1stClassName 2ndClassId 2ndClassName 3rdClassId 3rdClassName 3rdClassType pathwayId orthology enzymeAbbr enzymeName enzymeId geneId___- Metabolism - Global and overview maps - Metabolic pathways PATH ko01100 K00695 FKBP4_5 FK506-binding protein 4\/5 EC:5.2.1.8 LOC108847491
while(my $line=<FF>){
	chomp($line);
	($nameList, $valueList) = split(/___/, $line);
	@nameField = split(/\t/, $nameList);
	@valueField = split(/\t/, $valueList);
	for(my $i=0; $i<=$#nameField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$geneId = $tmpHash{"geneId"};
	$pathwayId = $tmpHash{"pathwayId"};
	$geneHref->{$geneId}->{"pathwayIdList"} .= $pathwayId . ",";
}
close FF;


# 将gene信息输出到mysqlTSV
my ($nameFieldString, $valueFieldString, @geneId, $mergedGoList, $mergedPfamList, $goImpact, $pfamImpact);
open WW, ">$geneMysqlTsv";
@geneId = keys(%gene);
foreach $geneId(@geneId){

############ nameFieldString ###############
	$nameFieldString = join("###", 
		"geneId", 
		"chr",
		"strand",
		"start",
		"stop",
		"geneSymbol", 
		"trsptList", 
		"A3ssList", 
		"A5ssList", 
		"MxeList", 
		"RiList", 
		"SeList", 
		"origGoList", 
		"origPfamList", 
		"specificHighTissue",
		"specificHighTissuePvalue",
		"specificHighTissueNum",
		"specificLowTissue",
		"specificLowTissuePvalue",
		"specificLowTissueNum",
		"koId",
		"enzymeId",
		"pathwayIdList",
		"annoSymbol",
		"trId",
		"trSymbol",
		"trFullName",
		"spId",
		"spSymbol",
		"spFullName"
);

	if(not exists($geneHref->{$geneId}->{"origGoList"})){
		$geneHref->{$geneId}->{"origGoList"} = "NA";
	}else{
		$geneHref->{$geneId}->{"origGoList"} = &sortTermList($geneHref->{$geneId}->{"origGoList"});
	}

	if(not exists($geneHref->{$geneId}->{"origPfamList"})){
		$geneHref->{$geneId}->{"origPfamList"} = "NA";
	}else{
		$geneHref->{$geneId}->{"origPfamList"} = &sortTermList($geneHref->{$geneId}->{"origPfamList"});
	}

	if(not exists($geneHref->{$geneId}->{"A3SSList"})){
		$geneHref->{$geneId}->{"A3SSList"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"A5SSList"})){
		$geneHref->{$geneId}->{"A5SSList"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"MXEList"})){
		$geneHref->{$geneId}->{"MXEList"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"RIList"})){
		$geneHref->{$geneId}->{"RIList"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"SEList"})){
		$geneHref->{$geneId}->{"SEList"} = "NA";
	}

	$geneHref->{$geneId}->{"geneSymbol"}=~tr/ /_/;

	# 是否为组织特异性
	if(not exists($geneHref->{$geneId}->{"specificHighTissue"})){
		$geneHref->{$geneId}->{"specificHighTissue"} = "NA";
		$geneHref->{$geneId}->{"specificHighTissueNum"} = -1;
		$geneHref->{$geneId}->{"specificHighTissuePvalue"} = "NA";
	}

	if(not exists($geneHref->{$geneId}->{"specificLowTissue"})){
		$geneHref->{$geneId}->{"specificLowTissue"} = "NA";
		$geneHref->{$geneId}->{"specificLowTissuePvalue"} = "NA";
		$geneHref->{$geneId}->{"specificLowTissueNum"} = -1;
	}

	# 是否具有KoId
	if(not exists($geneHref->{$geneId}->{"KId"})){
		$geneHref->{$geneId}->{"KId"} = "NA";
	}
	# 是否具有emzyme
	if(not exists($geneHref->{$geneId}->{"enzymeId"})){
		$geneHref->{$geneId}->{"enzymeId"} = "NA";
	}
	# 是否具有pathwayIdList
	if(not exists($geneHref->{$geneId}->{"pathwayIdList"})){
		$geneHref->{$geneId}->{"pathwayIdList"} = "NA";
	}else{
		$geneHref->{$geneId}->{"pathwayIdList"} = substr($geneHref->{$geneId}->{"pathwayIdList"}, 0, length($geneHref->{$geneId}->{"pathwayIdList"})-1);
	}

	# 是否具有swissprotDB的Symbol注释
	if(not exists($geneHref->{$geneId}->{"annoSymbol"})){
		$geneHref->{$geneId}->{"annoSymbol"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"trId"})){
		$geneHref->{$geneId}->{"trId"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"trSymbol"})){
		$geneHref->{$geneId}->{"trSymbol"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"trFullName"})){
		$geneHref->{$geneId}->{"trFullName"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"spId"})){
		$geneHref->{$geneId}->{"spId"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"spSymbol"})){
		$geneHref->{$geneId}->{"spSymbol"} = "NA";
	}
	if(not exists($geneHref->{$geneId}->{"spFullName"})){
		$geneHref->{$geneId}->{"spFullName"} = "NA";
	}

########## valuesFieldString #################
	$valueFieldString = join("###",
		$geneId, 
		$geneHref->{$geneId}->{"chr"},
		$geneHref->{$geneId}->{"strand"},
		$geneHref->{$geneId}->{"start"},
		$geneHref->{$geneId}->{"stop"},
		$geneHref->{$geneId}->{"geneSymbol"},
		$geneHref->{$geneId}->{"trsptList"}, 
		$geneHref->{$geneId}->{"A3SSList"}, 
		$geneHref->{$geneId}->{"A5SSList"}, 
		$geneHref->{$geneId}->{"MXEList"}, 
		$geneHref->{$geneId}->{"RIList"}, 
		$geneHref->{$geneId}->{"SEList"}, 
		$geneHref->{$geneId}->{"origGoList"}, 
		$geneHref->{$geneId}->{"origPfamList"}, 
		$geneHref->{$geneId}->{"specificHighTissue"},
		$geneHref->{$geneId}->{"specificHighTissuePvalue"},
		$geneHref->{$geneId}->{"specificHighTissueNum"},
		$geneHref->{$geneId}->{"specificLowTissue"},
		$geneHref->{$geneId}->{"specificLowTissuePvalue"},
		$geneHref->{$geneId}->{"specificLowTissueNum"},
		$geneHref->{$geneId}->{"KId"},
		$geneHref->{$geneId}->{"enzymeId"},
		$geneHref->{$geneId}->{"pathwayIdList"},
		$geneHref->{$geneId}->{"annoSymbol"},
		$geneHref->{$geneId}->{"trId"},
		$geneHref->{$geneId}->{"trSymbol"},
		$geneHref->{$geneId}->{"trFullName"},
		$geneHref->{$geneId}->{"spId"},
		$geneHref->{$geneId}->{"spSymbol"},
		$geneHref->{$geneId}->{"spFullName"}
);

	print WW $nameFieldString . "___" . $valueFieldString . "\n";
}
close WW;



# 排序termList
sub sortTermList{
	my ($termList) = @_;
	my (@term, @srtTerm, $term, $rltSrtTermList);
	@term = split(/,/, $termList);
	@srtTerm = sort(@term);
	foreach $term(@srtTerm){
		if($rltSrtTermList eq ""){
			$rltSrtTermList .= $term;
		}else{
			$rltSrtTermList .= "," . $term;
		}
	}
	return $rltSrtTermList;
}

# 合并两个termList，生成一个非冗余、有序的termlist
sub mergeTermList{
	my ($termList1, $termList2) = @_;
	my (%term, @term, $term, $termList);

	@term = ();
	@term = split(/,/, $termList1);
	foreach $term(@term){
		$term{$term} = 1;
	}

	@term = ();
	@term = split(/,/, $termList2);
	foreach $term(@term){
		$term{$term} = 1;
	}

	@term = ();
	@term = keys(%term);
	foreach $term(@term){
		$termList .= $term . ",";
	}

	$termList = substr($termList, 0, length($termList) - 1);
	return &sortTermList($termList);
}

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

sub removeAsIdFromAsIdList{
        my ($asIdList, $asId) = @_;
        my (@asId, $returnAsIdList, $tmpAsId);
        @asId = split(/,/, $asIdList);
        foreach $tmpAsId(@asId){
                $returnAsIdList .= $tmpAsId . "," if($tmpAsId ne $asId);
        }
        return substr($returnAsIdList, 0, length($returnAsIdList) - 1);
}

