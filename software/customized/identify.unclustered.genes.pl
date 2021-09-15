#!/usr/bin/perl
use strict;
use List::Util qw/max min/;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtfFile original.trspt.anno.gtf \\\n" .
		"--refGtfFileList ensembl.gtf \\\n" .
		"--refGffFileList refseq.gff3 \\\n" .
		"--gffChrToGtfChrFile refseq.seqId.map.ensembl.seqId.tsv \\\n" .
		"--RnaseqAssemblyListFile RnaseqAssemblyList.tsv \\\n" .
                "--unClusteredGeneList unClusterGeneList.tsv\n";
	exit;
}

my ($gtfFile, $unClusteredGeneList, $refGtfFileList, $refGffFileList, $gffChrToGtfChrFile,  $RnaseqAssemblyListFile);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
	'refGtfFileList=s'=>\$refGtfFileList,
	'refGffFileList=s'=>\$refGffFileList,
        'unClusteredGeneList=s'=>\$unClusteredGeneList,
	'gffChrToGtfChrFile=s'=>\$gffChrToGtfChrFile,
	'RnaseqAssemblyListFile=s'=>\$RnaseqAssemblyListFile,
);

my ($line, @field, $field, @attr, $attr, $geneId, $trsptId);
my (@refGtfFile, $refGtfFile, @refGffFile, $refGffFile, %refGeneCoord, %gffChrToGtfChr, $chr, $start, $chain, $stop);

# 读取RNAseq组装进入散列，登记experimentID -> assemblyFile 之间的映射关系
my (%rnaseqAssembly, @rnaseqAssembly);
open FF, "<$RnaseqAssemblyListFile";
#print $RnaseqAssemblyListFile;
while($line=<FF>){
	# /mnt/home/liujind1/workAS/9913/004-combine-assemblies-and-annos/../002-assemble-trsptome-on-goodExps/psiOutputDir/SRX2160810/transcriptomeByStringtie.gtf
	chomp($line);
	@field = ();
	@field = split(/\//, $line);
#	print $field[$#field-1];
#	<STDIN>;
	$rnaseqAssembly{$field[$#field-1]} = $line;
}
close FF;

# 读取gff3文件中chromosome编号和gtf中chromosome编号之间的映射关系
open FF, "<$gffChrToGtfChrFile";
while($line=<FF>){
	# NW_020192282.1  NKLS02002198.1
	chomp($line);
	@field = split(/\t/, $line);
	$gffChrToGtfChr{$field[0]} = $field[1];
}
close FF;


# 将gtf中的gene位置登记入散列
@refGtfFile = split(/,/, $refGtfFileList);
foreach $refGtfFile(@refGtfFile){
	open FF, "<$refGtfFile";
	while($line=<FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		next if($field[2] ne "gene");
		
		$geneId = "";
		@attr = ();
		@attr= split(/;/, $field[8]);
		foreach $attr(@attr){
			$geneId = $1 if($attr=~/gene_id "(.*)"/);
		}

		${$refGeneCoord{$geneId}}{"chr"} = $field[0];
		${$refGeneCoord{$geneId}}{"chain"} = $field[6];
		${$refGeneCoord{$geneId}}{"start"} = $field[3];
		${$refGeneCoord{$geneId}}{"stop"} = $field[4];
	}
	close FF;
}

# 将gff中的gene位置登记在散列中
@refGffFile = split(/,/, $refGffFileList);
foreach $refGffFile(@refGffFile){
	open FF, "<$refGffFile";
	while($line=<FF>){
                chomp($line);
                @field = ();
                @field = split(/\t/, $line);
                next if($field[2] ne "gene");

		$geneId = $1 if($field[8]=~/ID=(.*?);/);
                ${$refGeneCoord{$geneId}}{"chr"} = $gffChrToGtfChr{$field[0]};
                ${$refGeneCoord{$geneId}}{"chain"} = $field[6];
                ${$refGeneCoord{$geneId}}{"start"} = $field[3];
                ${$refGeneCoord{$geneId}}{"stop"} = $field[4];
	}
	close FF;
}


# 读取原始gtf文件，建立：
#（1）geneId - > trsptIdList的映射关系
#（2）trspt -> trspt position的映射关系
my (%geneToTrspt, %trspt);
open FF, "<$gtfFile";
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");

	$geneId = "";
	$trsptId = "";

	&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);
# 由于refseq的gtf和ensemb不一样，因此改动此处. 用&getGeneIdAndTrsptId替换
#	@attr = ();
#	@attr= split(/;/, $field[8]);
#	foreach $attr(@attr){
#		$geneId = $1 if($attr=~/gene_id "(.*)"/);
#		$trsptId = $1 if($attr=~/transcript_id "(.*)"/);
#		print "geneId=$geneId  trsptId=$trsptId\n";
#		<STDIN>;
#	}

	$geneToTrspt{$geneId} .= $trsptId . "#";;
	${$trspt{$trsptId}}{"chr"} = $field[0];
	${$trspt{$trsptId}}{"chain"} = $field[6];
	${$trspt{$trsptId}}{"start"} = $field[3];
	${$trspt{$trsptId}}{"stop"} = $field[4];
}
close FF;


# 读取每个geneId，然后在gene内部对trspt进行聚类
open WW, ">$unClusteredGeneList";
my (@geneId, @trsptId, @clusterId, $clusterId, %clusterToTrspt, $clusterNum, $i,  $j, $delClusterNum, @tmp);
@geneId = keys(%geneToTrspt);

foreach $geneId(@geneId){
	@trsptId = ();
	@trsptId = split(/#/, $geneToTrspt{$geneId});
	
	# 初始化每个trspt为一个cluster
	%clusterToTrspt = ();
	$clusterNum = 0;

	foreach $trsptId(@trsptId){
		$clusterNum++;
		$clusterId = "cluster" . $clusterNum;
		${$clusterToTrspt{$clusterId}}{"trsptList"} = $trsptId;
		${$clusterToTrspt{$clusterId}}{"chr"} = ${$trspt{$trsptId}}{"chr"};
		${$clusterToTrspt{$clusterId}}{"chain"} = ${$trspt{$trsptId}}{"chain"};
		${$clusterToTrspt{$clusterId}}{"start"} = ${$trspt{$trsptId}}{"start"};
		${$clusterToTrspt{$clusterId}}{"stop"} = ${$trspt{$trsptId}}{"stop"};
	}

	# 重复合并cluster，直到cluster的数量不再减少
	while(1){
		@clusterId = ();
		@clusterId = keys(%clusterToTrspt);
		$delClusterNum = 0;

		# 列出所有两个cluster组合，进行比较是否能够合并
		# 将 $clusterToTrspt{$i}合并到 $clusterToTrspt{$j} 中，
		# 如果能够合并，那么：
		#     将$clusterToTrspt{$i}的trsptList添加到$clusterToTrspt{$j}trsptList
		#     更新$clusterToTrspt{$j}的start和stop
		#     将$clusterToTrspt{$i}删除，$delClusterNum++
		for($i=0; $i<=$#clusterId-1; $i++){
			for($j=$i+1; $j<=$#clusterId; $j++){

				# 不   分离情况1                                                                                    分离情况2
				if(not(${$clusterToTrspt{$clusterId[$i]}}{"stop"} < ${$clusterToTrspt{$clusterId[$j]}}{"start"} or ${$clusterToTrspt{$clusterId[$i]}}{"start"} > ${$clusterToTrspt{$clusterId[$j]}}{"stop"})){

					# 更新$clusterToTrspt{$clusterId[$j]}的坐标
					${$clusterToTrspt{$clusterId[$j]}}{"stop"} = max (${$clusterToTrspt{$clusterId[$i]}}{"start"}, ${$clusterToTrspt{$clusterId[$i]}}{"stop"}, ${$clusterToTrspt{$clusterId[$j]}}{"start"}, ${$clusterToTrspt{$clusterId[$j]}}{"stop"});
					${$clusterToTrspt{$clusterId[$j]}}{"start"} = min (${$clusterToTrspt{$clusterId[$i]}}{"start"}, ${$clusterToTrspt{$clusterId[$i]}}{"stop"}, ${$clusterToTrspt{$clusterId[$j]}}{"start"}, ${$clusterToTrspt{$clusterId[$j]}}{"stop"});

					# 更新$clusterToTrspt{$clusterId[$j]}的trsptId
					${$clusterToTrspt{$clusterId[$j]}}{"trsptList"} .= "#" . ${$clusterToTrspt{$clusterId[$i]}}{"trsptList"};

					# 删除$clusterToTrspt{$clusterId[$i]}
					delete($clusterToTrspt{$clusterId[$i]});
					$delClusterNum++;
					last;

				}
			}
		}

		last if($delClusterNum == 0);
	}

	# 输出当前gene中所有的cluster：
	# gene1 \t cluster1:chr1#+#1222000#1345000,trspt1:chr#+#12222000#134000,trspt2:chr1#+#123232-124234 \t cluster2:chr1#+#123434#1324234,trspt4:chr1#+#234234#213423,trspt5:chr1#+#345345#8734
	
	@clusterId = ();
	@clusterId = keys(%clusterToTrspt);

	if($#clusterId >=1){
		if(exists($refGeneCoord{$geneId})){
			$chr = ${$refGeneCoord{$geneId}}{"chr"};
			$chain = ${$refGeneCoord{$geneId}}{"chain"};
			$start = ${$refGeneCoord{$geneId}}{"start"};
			$stop = ${$refGeneCoord{$geneId}}{"stop"};
		}elsif($geneId=~/(.*)\.(\d+)/){
			my $expId = $1;
			$chr = "";
			$chain = "";
			$start = 0;
			$stop = 0;
			&getRnaseqAsseblyCoord($rnaseqAssembly{$expId}, $geneId, \$chr, \$chain, \$start, \$stop);
		}

		print WW $geneId . ":" . $chr . "#" . $chain . "#" . $start . "#" . $stop;

		foreach $clusterId(@clusterId){
			print WW "\t";
			print WW $clusterId . ":" . ${$clusterToTrspt{$clusterId}}{"chr"} . "#" . ${$clusterToTrspt{$clusterId}}{"chain"} . "#" . ${$clusterToTrspt{$clusterId}}{"start"} . "#" . ${$clusterToTrspt{$clusterId}}{"stop"};

			@trsptId = ();
			@trsptId = split(/#/, ${$clusterToTrspt{$clusterId}}{"trsptList"});
			foreach $trsptId(@trsptId){
				print WW ",";
				print WW $trsptId . ":" . ${$trspt{$trsptId}}{"chr"} . "#" . ${$trspt{$trsptId}}{"chain"} . "#" . ${$trspt{$trsptId}}{"start"} . "#" . ${$trspt{$trsptId}}{"stop"};
			}
		}
		print WW "\n";
	}

}
close WW;

sub getRnaseqAsseblyCoord{
	my ($assemblyFile, $geneId, $chr, $chain, $start, $stop) = @_;
	my ($line, @field, $field, @attr, $attr);

	open FFF, "<$assemblyFile";
	while($line=<FFF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);

		next if($field[2] ne "transcript");
		if($field[8]=~/gene_id "(.*?)";/){
			next if($1 ne $geneId);
		}

		if($$start == 0){
			$$chr = $field[0];
			$$chain = $field[6];
			$$start = $field[3];
			$$stop = $field[4];
		}else{
			$$start = $field[3] if($field[3] < $$start);
			$$stop = $field[4] if($field[4] > $$stop);
		}
	}
	close FFF;
}

sub getGeneIdAndTrsptId{
	my ($attrString, $geneId, $trsptId) = @_;
	my (@attr, $attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}
	}
}
