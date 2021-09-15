#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputAsFileList \\\n" .
                "--taxonIdList \\\n" .
		"--criterion whole/inner \\\n" .
                "--asType \n" .
                "--outputOrthoAsStatisticFile \\\n" .
                "--outputOrthoAsListFile \n";
	exit;
}

my ($inputAsFileList, $taxonIdList, $asType, $criterion, $outputOrthoAsStatisticFile, $outputOrthoAsListFile);

GetOptions(
        'inputAsFileList=s'=>\$inputAsFileList,
        'taxonIdList=s'=>\$taxonIdList,
        'asType=s'=>\$asType,
	'criterion=s'=>\$criterion,
        'outputOrthoAsStatisticFile=s'=>\$outputOrthoAsStatisticFile,
        'outputOrthoAsListFile=s'=>\$outputOrthoAsListFile,
);

# 打开as文件将as读入hash：以全部/内部坐标组合为关键字，以taxonId为第2级关键字，以AS编号为列表为值
my ($i, $j, $line, @nameField, @valueField, %as, $asHref);
my (@inputAsFile, $inputAsFile, @taxonId, $taxonId);
my (%coordToTaxonToAsIdList, $coordToTaxonToAsIdListHref, $coordSeriesKey, %coordToTaxonToAsIdNum, $coordToTaxonToAsIdNumHref);
$coordToTaxonToAsIdListHref = \%coordToTaxonToAsIdList;
$coordToTaxonToAsIdNumHref = \%coordToTaxonToAsIdNum;
@inputAsFile = split(/,/, $inputAsFileList);
@taxonId = split(/,/, $taxonIdList);
for($i=0; $i<=$#taxonId; $i++){
	$taxonId = $taxonId[$i];
	$inputAsFile = $inputAsFile[$i];
	
	# 打开一个asFile
	open FF, "<$inputAsFile";
	$line = <FF>;
	chomp($line);
	@nameField = ();
	@nameField = split(/\t/, $line);
	while($line=<FF>){
		chomp($line);
		# 将1个as读入hash
		@valueField = ();
		@valueField = split(/\t/, $line);
		%as = ();
		$asHref = \%as;
		for($j=0; $j<=$#valueField; $j++){
			$asHref->{$nameField[$j]} = $valueField[$j];
		}

		# 根据as的类型和比较标准（全部坐标/内部坐标)生成coordToTaxonToAsIdList的关键字
		$coordSeriesKey = &getAsCoordSeriesKey($asHref, $asHref->{"strand"}, $asType, $criterion);
		my $checkRltt = 0;
		my @tt = ();
		@tt = split(/#/, $coordSeriesKey);
		if($criterion eq "all"){
			# 当该as的内部坐标都在orthoCdsAlignPos上时才登记到hash中
			# 检查AS是否存不在orthoCdsAlignPos上的内部坐标
			# AS至少需要4个坐标描述其位置
			# 如果内部坐标存在"-"，那么放弃该orthoAs
			# 也就是允许AS边缘的2个坐标为-

			for(my $j=1; $j<=$#tt-1; $j++){
				$checkRltt = 1 if($tt[$j] eq "-");
			}
		}else{
			for(my $j=0; $j<=$#tt; $j++){
				if($tt[$j] eq "-"){
					$checkRltt = 1;
				}
			}
		}
		if($checkRltt == 0){
			$coordToTaxonToAsIdNumHref->{$coordSeriesKey}->{$taxonId}++;	
			if(not exists($coordToTaxonToAsIdListHref->{$coordSeriesKey}->{$taxonId})){
				$coordToTaxonToAsIdListHref->{$coordSeriesKey}->{$taxonId} = $asHref->{"ASID"};
			}else{
				$coordToTaxonToAsIdListHref->{$coordSeriesKey}->{$taxonId} .= "," . $asHref->{"ASID"};
			}
		}
	}
	close FF;
	# 
}
close STAT;
close LIST;

# 对orthoAs统一编号
my ($orthoAsNum, @coordSeriesKey, $outputStat, $outputList);
open STAT, ">$outputOrthoAsStatisticFile";
print STAT join("\t", "orthAsId", @taxonId) . "\n";
open LIST, ">$outputOrthoAsListFile";
print LIST join("\t", "orthAsId", @taxonId) . "\n";
@coordSeriesKey = keys(%coordToTaxonToAsIdList);
foreach $coordSeriesKey(@coordSeriesKey){
	$orthoAsNum++;
	$outputStat = "orth" . $asType . sprintf("%08d", $orthoAsNum);
	$outputList = "orth" . $asType . sprintf("%08d", $orthoAsNum);
	foreach $taxonId(@taxonId){
		if(not exists($coordToTaxonToAsIdListHref->{$coordSeriesKey}->{$taxonId})){
			$outputStat .= "\t0";
			$outputList .= "\t-"; 
		}else{
			$outputStat .= "\t" . $coordToTaxonToAsIdNumHref->{$coordSeriesKey}->{$taxonId};
			$outputList .= "\t" . $coordToTaxonToAsIdListHref->{$coordSeriesKey}->{$taxonId};
		}
	}
	print STAT $outputStat . "\n";
	print LIST $outputList . "\n";
}
close STAT;
close LIST;

# 读取as的坐标，按照正义cDNA的次序生成坐标组合
# 即：按照cDNA的5'->3'次序
# longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
sub getAsCoordSeriesKey{
	my ($asHref, $strand, $asType, $criterion) = @_;
	if($asType eq "A3SS"){
		# longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
		if($strand eq "+"){
			if($criterion eq "all"){
				return join("#", $asHref->{"flankingES"}, $asHref->{"flankingEE"}, $asHref->{"longExonStart_0base"}, $asHref->{"shortES"}, $asHref->{"longExonEnd"});
			}else{
				return join("#", $asHref->{"flankingEE"}, $asHref->{"longExonStart_0base"}, $asHref->{"shortES"});
			}
		}else{
			if($criterion eq "all"){
				return join("#", $asHref->{"flankingEE"}, $asHref->{"flankingES"}, $asHref->{"longExonEnd"}, $asHref->{"shortEE"}, $asHref->{"longExonStart_0base"});
			}else{
				return join("#", $asHref->{"flankingES"}, $asHref->{"longExonEnd"}, $asHref->{"shortEE"});
			}
		}
	}elsif($asType eq "A5SS"){
		# longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
		if($strand eq "+"){
			if($criterion eq "all"){
				return join("#", $asHref->{"longExonStart_0base"}, $asHref->{"shortEE"}, $asHref->{"longExonEnd"}, $asHref->{"flankingES"}, $asHref->{"flankingEE"});
			}else{
				return join("#", $asHref->{"shortEE"}, $asHref->{"longExonEnd"}, $asHref->{"flankingES"});
			}
		}else{
			if($criterion eq "all"){
				return join("#", $asHref->{"longExonEnd"}, $asHref->{"shortES"}, $asHref->{"longExonStart_0base"}, $asHref->{"flankingEE"}, $asHref->{"flankingES"});
			}else{
				return join("#", $asHref->{"shortES"}, $asHref->{"longExonStart_0base"}, $asHref->{"flankingEE"});
			}
		}
	}elsif($asType eq "MXE"){
		# 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
		if($strand eq "+"){
			if($criterion eq "all"){
				return join("#", $asHref->{"upstreamES"}, $asHref->{"upstreamEE"}, $asHref->{"1stExonStart_0base"}, $asHref->{"1stExonEnd"}, $asHref->{"2ndExonStart_0base"}, $asHref->{"2ndExonEnd"}, $asHref->{"downstreamES"}, $asHref->{"downstreamEE"});
			}else{
				return join("#", $asHref->{"upstreamEE"}, $asHref->{"1stExonStart_0base"}, $asHref->{"1stExonEnd"}, $asHref->{"2ndExonStart_0base"}, $asHref->{"2ndExonEnd"}, $asHref->{"downstreamES"});
			}
		}else{
			if($criterion eq "all"){
				return join("#", $asHref->{"downstreamEE"}, $asHref->{"downstreamES"}, $asHref->{"2ndExonEnd"}, $asHref->{"2ndExonStart_0base"}, $asHref->{"1stExonEnd"}, $asHref->{"1stExonStart_0base"}, $asHref->{"upstreamEE"}, $asHref->{"upstreamES"});
			}else{
				return join("#", $asHref->{"downstreamES"}, $asHref->{"2ndExonEnd"}, $asHref->{"2ndExonStart_0base"}, $asHref->{"1stExonEnd"}, $asHref->{"1stExonStart_0base"}, $asHref->{"upstreamEE"});
			}
		}
	}elsif($asType eq "RI"){
		# riExonStart_0base       riExonEnd       upstreamES      upstreamEE      downstreamES    downstreamEE
		if($strand eq "+"){
			if($criterion eq "all"){
				return join("#", $asHref->{"upstreamES"}, $asHref->{"upstreamEE"}, $asHref->{"downstreamES"}, $asHref->{"downstreamEE"});
			}else{
				return join("#", $asHref->{"upstreamEE"}, $asHref->{"downstreamES"});
			}
		}else{
			if($criterion eq "all"){
				return join("#", $asHref->{"downstreamEE"}, $asHref->{"downstreamES"}, $asHref->{"upstreamEE"}, $asHref->{"upstreamES"});
			}else{
				return join("#", $asHref->{"downstreamES"}, $asHref->{"upstreamEE"});
			}
		}
	}elsif($asType eq "SE"){
		# exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE
		if($strand eq "+"){
			if($criterion eq "all"){
				return join("#", $asHref->{"upstreamES"}, $asHref->{"upstreamEE"}, $asHref->{"exonStart_0base"}, $asHref->{"exonEnd"}, $asHref->{"downstreamES"}, $asHref->{"downstreamEE"});
			}else{
				return join("#", $asHref->{"upstreamEE"}, $asHref->{"exonStart_0base"}, $asHref->{"exonEnd"}, $asHref->{"downstreamES"});
			}
		}else{
			if($criterion eq "all"){
				return join("#", $asHref->{"downstreamEE"}, $asHref->{"downstreamES"}, $asHref->{"exonEnd"}, $asHref->{"exonStart_0base"}, $asHref->{"upstreamEE"}, $asHref->{"upstreamES"});
			}else{
				return join("#", $asHref->{"downstreamES"}, $asHref->{"exonEnd"}, $asHref->{"exonStart_0base"}, $asHref->{"upstreamEE"});
			}
		} 
	}
}

