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

# SRX2160810.81987:NKLS02000075.1#+#7041#33013    cluster2:NKLS02000075.1#+#15608#32416,SRX2160810.81987.8:NKLS02000075.1#+#15608#32416   cluster1:NKLS02000075.1#+#7620#15493,SRX2160809.1092628.26:NKLS02000075.1#+#7620#15493
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
	
	# 对每个cluster进行，检查是否都是和基因位置有重叠
	my $allOverlap = 1;
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
		$allOverlap = 0 if($geneStart > $clusterStop or $geneStop < $clusterStart);
	}

	# 如果所有的转录本簇不是都和基因重叠，那么放弃重命名基因id
	next if($allOverlap == 0);

	# 如果所有的转录本簇都和基因重叠，那么重新命名这些簇内转录本的基因ID
	# 格式是:geneId.inner.cluster。
	# 新命名的基因ID不在采用原来关联的ID了.
	
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

		$newGeneId = $geneId . ".inner." . $clusterId;

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
my ($geneId, $trsptId, $attrString, @field);
open FF, "<$originalGtfFile";
open WW, ">$outputNewGtfFile";
while($line=<FF>){
        @field = ();
        @field = split(/\t/, $line);
        $attrString = $field[8];
        &getGeneIdAndTrsptId($attrString, \$geneId, \$trsptId);
        if(exists($trsptIdToGeneId{$trsptId})){
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

