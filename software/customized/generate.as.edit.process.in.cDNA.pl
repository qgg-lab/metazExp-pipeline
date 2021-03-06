#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asMappingTrsptTsv \\\n" .
                "--genomeFasta \\\n" .
                "--asProcessInCdnaTsv \n";
	exit;
}
# 注意：
my ($asMappingTrsptTsv, $genomeFasta, $asProcessInCdnaTsv);

GetOptions(
        'asMappingTrsptTsv=s'=>\$asMappingTrsptTsv,
        'genomeFasta=s'=>\$genomeFasta,
        'asProcessInCdnaTsv=s'=>\$asProcessInCdnaTsv,
);

# 读取基因组序列如散列
my (%genomeSeq, $chr, $line, @tt);
open FF, "<$genomeFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$chr = $1;
		@tt = ();
		@tt = split(/ /, $chr);
		$chr = $tt[0];
	}else{
		$genomeSeq{$chr}.=uc($line);
	}
}
close FF;

open WW, ">$asProcessInCdnaTsv";
print WW join("\t", "chr", "strand", "asId", "asType", "trsptId", "residentType", "editSite", "cutSize", "insertSize", "insertSeq") . "\n";
# 读取AS和Trspt之间映射关系，将其读入到hash
my (@titleField, $titleField, @valueField, $valueField, $mappingLine, $i, @exon, @exonString, $exonString, @tt, $exonNum);
my (%asMapTrspt, $asMapTrsptHref, $editSite, $cutSize, $insertSize, $insertSeq, $editBegin, $editEnd, $cutSeq, $residentType);
$asMapTrsptHref = \%asMapTrspt;
open FF, "<$asMappingTrsptTsv";
$mappingLine=<FF>;
chomp($mappingLine);
@titleField = split(/\t/, $mappingLine);
while($mappingLine=<FF>){
	chomp($mappingLine);
	@valueField = ();
	@valueField = split(/\t/, $mappingLine);
	# 将这一行值读入到hash中
	for($i=0; $i<=$#titleField; $i++){
		$asMapTrsptHref->{$titleField[$i]} = $valueField[$i];
	}
	
	# 生成外显子坐标二维数组
	@exonString = ();
	@exonString = split(/,/, $asMapTrsptHref->{"trsptExonSeries"});
	$exonNum = 0;
#	print "ASId:" . $asMapTrsptHref->{"asId"} . "\n";
#	print "ASType:" . $asMapTrsptHref->{"asType"} . "\n";
#	print "trsptId：" . $asMapTrsptHref->{"trsptId"} . "\n";
#	print "trsptExonSeries:" . $asMapTrsptHref->{"trsptExonSeries"} . "\n";
	foreach $exonString(@exonString){
		@tt = ();
		@tt = split(/\.\./, $exonString);
		$exon[$exonNum][0] = $tt[0];
		$exon[$exonNum][1] = $tt[1];
#		print $exon[$exonNum][0] . "\t" . $exon[$exonNum][1] . "\n";
		$exonNum++;
	}
#	<STDIN>;
	# 计算AS对trspt的操作方式: ASID, ASType, trsptId, editSite, cutSize, insertSize, insertSeq
	# 注意：editSite在操作时并不涉及，涉及的只是它之后的碱基

	if($asMapTrsptHref->{"residentType"} eq "shortAlt"){

###### Trspt和AS的关系：Trspt包含了AS的shortAlt部分 ########
		$editSite = 0;
		$cutSize = 0;
		$insertSize = 0;
		$insertSeq = "";
		# A3SS
		if($asMapTrsptHref->{"asType"} eq "A3SS"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"flankingEE"});
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"longExonStart_0base"}, $asMapTrsptHref->{"shortES"} - ($asMapTrsptHref->{"longExonStart_0base"} + 1) + 1);
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"flankingES"}+1);
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"shortEE"}, $asMapTrsptHref->{"longExonEnd"} - ($asMapTrsptHref->{"shortEE"} +1) + 1);
				$insertSeq=reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
			}
		}

		# A5SS
		if($asMapTrsptHref->{"asType"} eq "A5SS"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"shortEE"});
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"shortEE"}, $asMapTrsptHref->{"longExonEnd"} - ($asMapTrsptHref->{"shortEE"} + 1) + 1);
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"shortES"}+1);
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"longExonStart_0base"}, $asMapTrsptHref->{"shortES"} - ($asMapTrsptHref->{"longExonStart_0base"} + 1) + 1);
				$insertSeq = reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
			}
		}

		# RI
		if($asMapTrsptHref->{"asType"} eq "RI"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"upstreamEE"}, $asMapTrsptHref->{"downstreamES"} - ($asMapTrsptHref->{"upstreamEE"} + 1) + 1);
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"upstreamEE"}, $asMapTrsptHref->{"downstreamES"} - ($asMapTrsptHref->{"upstreamEE"} + 1) + 1);
				$insertSeq = reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
			}
		}

		# SE
		if($asMapTrsptHref->{"asType"} eq "SE"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"exonStart_0base"}, $asMapTrsptHref->{"exonEnd"} - ($asMapTrsptHref->{"exonStart_0base"} + 1) + 1);
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"exonStart_0base"}, $asMapTrsptHref->{"exonEnd"} - ($asMapTrsptHref->{"exonStart_0base"} + 1) + 1);
				$insertSeq = reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
			}
		}

		# MXE
		if($asMapTrsptHref->{"asType"} eq "MXE"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"1stExonStart_0base"}, $asMapTrsptHref->{"1stExonEnd"} - ($asMapTrsptHref->{"1stExonStart_0base"} + 1) + 1);
				$cutSize = $asMapTrsptHref->{"2ndExonEnd"} - $asMapTrsptHref->{"2ndExonStart_0base"};
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"1stExonStart_0base"}, $asMapTrsptHref->{"1stExonEnd"} - ($asMapTrsptHref->{"1stExonStart_0base"} + 1) + 1);
				$insertSeq = reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
				$cutSize = $asMapTrsptHref->{"2ndExonEnd"} - $asMapTrsptHref->{"2ndExonStart_0base"};
			}
		}

		$insertSize = length($insertSeq);
		print WW join("\t", $asMapTrsptHref->{"chr"}, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"asId"}, $asMapTrsptHref->{"asType"}, $asMapTrsptHref->{"trsptId"}, $asMapTrsptHref->{"residentType"}, $editSite, $cutSize, $insertSize, $insertSeq) . "\n";

	}else{
#### Trspt和AS的关系：Trspt包含了AS的longAlt部分  #######
		$editSite = 0;
		$insertSize = 0;
		$insertSeq = "";
		$cutSize = 0;

		# A3SS
		if($asMapTrsptHref->{"asType"} eq "A3SS"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"flankingEE"});
				$cutSize = $asMapTrsptHref->{"shortES"} - ($asMapTrsptHref->{"longExonStart_0base"} + 1) + 1;
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"flankingES"} + 1);
				$cutSize = $asMapTrsptHref->{"longExonEnd"} - ($asMapTrsptHref->{"shortEE"} + 1) + 1;
			}
		}

		# A5SS
		if($asMapTrsptHref->{"asType"} eq "A5SS"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"shortEE"});
				$cutSize = $asMapTrsptHref->{"longExonEnd"} - ($asMapTrsptHref->{"shortEE"} + 1) + 1;
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"shortES"} + 1);
				$cutSize = $asMapTrsptHref->{"shortES"} - ($asMapTrsptHref->{"longExonStart_0base"} + 1) + 1;
			}
		}

		# RI
		if($asMapTrsptHref->{"asType"} eq "RI"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$cutSize = $asMapTrsptHref->{"downstreamES"} - ($asMapTrsptHref->{"upstreamEE"} + 1) + 1;	
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$cutSize = $asMapTrsptHref->{"downstreamES"} - ($asMapTrsptHref->{"upstreamEE"} + 1) + 1;
			}
		}	

		# SE
		if($asMapTrsptHref->{"asType"} eq "SE"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$cutSize = $asMapTrsptHref->{"exonEnd"} - ($asMapTrsptHref->{"exonStart_0base"} + 1) + 1;
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$cutSize = $asMapTrsptHref->{"exonEnd"} - ($asMapTrsptHref->{"exonStart_0base"} + 1) + 1;
			}
		}

		# MXE
		if($asMapTrsptHref->{"asType"} eq "MXE"){
			if($asMapTrsptHref->{"strand"} eq "+"){
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"upstreamEE"});
				$cutSize = $asMapTrsptHref->{"1stExonEnd"} - ($asMapTrsptHref->{"1stExonStart_0base"} + 1) + 1;
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"2ndExonStart_0base"}, $asMapTrsptHref->{"2ndExonEnd"} - ($asMapTrsptHref->{"2ndExonStart_0base"} + 1) + 1);
				$insertSize = length($insertSeq);
			}else{
				$editSite = &getCoordinateInCdna(\@exon, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"downstreamES"} + 1);
				$cutSize = $asMapTrsptHref->{"1stExonEnd"} - ($asMapTrsptHref->{"1stExonStart_0base"} + 1) + 1;
				$insertSeq = substr($genomeSeq{$asMapTrsptHref->{"chr"}}, $asMapTrsptHref->{"2ndExonStart_0base"}, $asMapTrsptHref->{"2ndExonEnd"} - ($asMapTrsptHref->{"2ndExonStart_0base"} + 1) + 1);
				$insertSeq = reverse($insertSeq);
				$insertSeq=~tr/ACGT/TGCA/;
				$insertSize = length($insertSeq);
			}
		}
	
		$insertSize = length($insertSeq);
		print WW join("\t", $asMapTrsptHref->{"chr"}, $asMapTrsptHref->{"strand"}, $asMapTrsptHref->{"asId"}, $asMapTrsptHref->{"asType"}, $asMapTrsptHref->{"trsptId"}, $asMapTrsptHref->{"residentType"}, $editSite, $cutSize, $insertSize, $insertSeq) . "\n";
	}
}

close FF;

# 转录本外显子在DNA中的坐标以及所在链，转换1个DNA坐标为转录本中的坐标
sub getCoordinateInCdna{
        my ($exonInDna, $strand, $siteInDna) = @_;
        my ($exonNum, $position, @tt, $i);

        $exonNum = $#$exonInDna;
        if($strand eq "+"){
                for($i=0; $i<=$exonNum; $i++){
                        if($siteInDna > $$exonInDna[$i][1]){
                                $position+=$$exonInDna[$i][1] - $$exonInDna[$i][0] + 1;
                        }else{
                                $position+= $siteInDna - $$exonInDna[$i][0] + 1;
#				print join("\t", "siteInDna:" . $siteInDna, "siteIncDNA" . $position) . "\n";
#				<STDIN>;
                                return $position;
                        }
                }
        }else{
                for($i=0; $i<=$exonNum; $i++){
                        if($siteInDna < $$exonInDna[$i][0]){
                                $position+=$$exonInDna[$i][1] - $$exonInDna[$i][0] + 1;
                        }else{
                                $position+= $$exonInDna[$i][1] - $siteInDna + 1;
#				print join("\t", "siteInDna:" . $siteInDna, "siteIncDNA" . $position) . "\n";
#				<STDIN>;
                                return $position;
                        }
                }
        }
}

