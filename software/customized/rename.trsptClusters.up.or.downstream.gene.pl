#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--unClusteredGeneList unClusterGeneList.tsv \\\n" .
		"--originalGtfFile origin.trspt.gtf \\\n" .
		"--outputNewGtfFile new.trspt.gtf \\\n" .
                "--outputRenameTrspt renamed.trsptId.list \n";
	exit;
}

my ($unClusteredGeneList, $outputRenameTrspt, $originalGtfFile, $outputNewGtfFile);
my ($line);
my ($geneId, $geneChain, $geneStart, $geneStop, $geneChr);
my (@cluster, $cluster, $clusterId, $clusterChr, $clusterChain, $clusterStart, $clusterStop);
my ($genePosition,  $clusterPosition);
my (@trspt, $trspt, $trsptId, $newGeneId);
my (%trsptIdToGeneId);

GetOptions(
        'unClusteredGeneList=s'=>\$unClusteredGeneList,
	'originalGtfFile=s'=>\$originalGtfFile,
	'outputNewGtfFile=s'=>\$outputNewGtfFile,
        'outputRenameTrspt=s'=>\$outputRenameTrspt,
);

open WW, ">$outputRenameTrspt";
open FF, "<$unClusteredGeneList";
# ENSBTAG00000017895:19#-#26833438#26837327       cluster1:19#-#26825722#26828000,SRX1354441.16485.6:19#-#26825722#26828000,SRX1968515.28662.7:19#-#26825722#26828000,SRX1354421.16327.5:19#-#26825722#26828000   cluster7:19#-#26828128#26837370,ERX1545694.40145.1:19#-#26833450#26834992,ENSBTAT00000054378:19#-#26833438#26837327,rna54636:19#-#26828128#26837368,rna54635:19#-#26828128#26837370,ERX1545690.34785.2:19#-#26835316#26837327
while($line=<FF>){
	chomp($line);
	@cluster = ();
	@cluster = split(/\t/, $line);

	# 提取整个基因的位置范围
	$genePosition = shift(@cluster);
	if($genePosition=~/(.*):(.*)#(\+|\-)#(\d+)#(\d+)/){
		$geneId = $1;
		$geneChr = $2;
		$geneChain = $3;
		$geneStart = $4;
		$geneStop = $5;
	}
	
	# 对每个cluster进行分析
	foreach $cluster(@cluster){
		@trspt = ();
		@trspt = split(/,/, $cluster);

		# 提取该cluster的位置范围
		$clusterPosition = shift(@trspt);

		if($clusterPosition=~/(.*):(.*)#(\+|\-)#(\d+)#(\d+)/){
	                $clusterId = $1;
        	        $clusterChr = $2;
	                $clusterChain = $3;
        	        $clusterStart = $4;
	                $clusterStop = $5;
		}

		# 如果重叠，那么跳过去不处理该cluster中的trspt的geneId，仍然用原来的gene ID; 
		# 否则重新命名该cluster中的转录本geneId，格式：临近geneId.up/down.clusterId
		next if(not($geneStart > $clusterStop or $geneStop < $clusterStart));

		if($geneStart > $clusterStop){
			$newGeneId = $geneId . ".up." . $clusterId;
		}elsif($geneStop < $clusterStart){
			$newGeneId = $geneId . ".down." . $clusterId;
		}

		# 抽出trsptId
		foreach $trspt(@trspt){
			if($trspt=~/(.*):(.*)#(\+|\-)#(\d+)#(\d+)/){
				$trsptId = $1;
				print WW $trsptId . "\t" . $newGeneId . "\n";
				$trsptIdToGeneId{$trsptId} = $newGeneId;
			}
		}
	}
}
close FF;
close WW;
#$originalGtfFile, $outputNewGtfFile
my ($geneId, $trsptId, @field, $attrString);
open FF, "<$originalGtfFile";
open WW, ">$outputNewGtfFile";
while($line=<FF>){
	# 1       StringTie       transcript      195058  207560  1000    +       .       gene_id "SRX1481222.1"; transcript_id "SRX1481222.1.1"; cov "4.421379"; FPKM "0.230956"; TPM "0.541498";
	($geneId, $trsptId) = ("", "");
	@field = ();
	@field = split(/\t/, $line);
	$attrString = $field[8];
	&getGeneIdAndTrsptId($attrString, \$geneId, \$trsptId);
	
	
	if(exists($trsptIdToGeneId{$trsptId})){
		# 用新的geneId替换掉原来的geneId 
		pop(@field);
		$geneId = $trsptIdToGeneId{$trsptId};
		if($attrString=~/gene_id "$geneId"; transcript_id "$trsptId"; (.*)\n/){
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\"; " . $1 . "\n");

		}elsif($attrString=~/gene_id "$geneId"; transcript_id "$trsptId"\n/){
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";\n");
		}elsif($attrString=~/gene_id "$geneId"; transcript_id "$trsptId";\n/){
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";\n");
		}
	}else{
		if($attrString=~/gene_id "$geneId"; transcript_id "$trsptId"\n/){
			pop(@field);
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";\n") ;
		}else{
			print WW $line;
		}

	}
}
close FF;
close WW;

sub getGeneIdAndTrsptId{
        my ($attrString, $geneId, $trsptId) = @_;
        my (@attr, $attr);
	chomp($attrString);
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

