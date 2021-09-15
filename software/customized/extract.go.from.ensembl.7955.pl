#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputGtf \\\n" .
                "--inputEmblAnnoData \\\n" .
                "--outputGeneGoAnn \\\n" .
		"--outputTrsptGoAnn \n";
	exit;
}

my ($inputGtf, $inputEmblAnnoData, $outputGeneGoAnn, $outputTrsptGoAnn);

GetOptions(
        'inputGtf=s'=>\$inputGtf,
        'inputEmblAnnoData=s'=>\$inputEmblAnnoData,
        'outputGeneGoAnn=s'=>\$outputGeneGoAnn,
        'outputTrsptGoAnn=s'=>\$outputTrsptGoAnn,
);

my (%trsptTogene, %geneToTrsptIdList, @geneId, $geneId, @trsptId, $trsptId, @field, $field, @attr, $attr, $line, $trsptVersion);
# 建立gene和transcript之间的映射关系
open FF, "<$inputGtf";
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	next if($#field!=8);
	($geneId, $trsptId) = ("", "");
	@attr = ();
	@attr = split(/;/, $field[8]);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*?)"/){
			$geneId = $1;
		}
		if($attr=~/transcript_id "(.*?)"/){
			$trsptId = $1;
		}
		if($attr=~/transcript_version "(\d+)"/){
			$trsptVersion = $1;
		}
	}
	if($geneId ne "" and $trsptId ne ""){
		#$trsptTogene{$trsptId . "." . $trsptVersion} = $geneId;
		$trsptTogene{$trsptId}= $geneId;
		#print join("\t", $trsptId, $geneId);
	}
}
close FF;

@trsptId = keys(%trsptTogene);
foreach $trsptId(@trsptId){
	$geneToTrsptIdList{$trsptTogene{$trsptId}} .= $trsptId . ",";
}

# 将embl注释整体读入文本变量
my $emblTxt = "";
open FF, "<$inputEmblAnnoData";
while($line=<FF>){
	$emblTxt .= $line;
}
close FF;


# 以转录本为单位对文本进行分隔，使之成为数组
my (@trsptTxt, $trsptTxt, $go, $goList, @line,%trsptToGoList);
@trsptTxt = split(/FT   mRNA            /, $emblTxt);
foreach $trsptTxt(@trsptTxt){
	($trsptId, $goList) = ("", "");	

	@line = ();
	@line = split(/\n/, $trsptTxt);

	foreach $line(@line){
		if($line=~/FT +\/standard_name="(.*)"/){
			$trsptId = $1;
		}
		if($line=~/FT +\/db_xref="(GO:\d+)"/){
			$go = $1;
			if(index($goList, $go)<0){
				$goList .= $go . ",";
			}
		}
	}

	next if($trsptId eq "" or $goList eq "");
	$trsptToGoList{$trsptId} = substr($goList, 0, length($goList) - 1);
	#print $trsptId . "\t" . $trsptToGoList{$trsptId};
	#<STDIN>;
}

# 输出transcript的go注释
open WW, ">$outputTrsptGoAnn";
@trsptId = keys(%trsptToGoList);
foreach $trsptId(@trsptId){
	print WW join("\t", $trsptId, $trsptToGoList{$trsptId}) . "\n";
}
close WW;


# 输出gene的go注释
my (@go);
open WW, ">$outputGeneGoAnn";
@geneId = keys(%geneToTrsptIdList);
foreach $geneId(@geneId){
	$goList = "";
	@trsptId = ();
	@trsptId = split(/,/, $geneToTrsptIdList{$geneId});
	foreach $trsptId(@trsptId){
		next if(not exists($trsptToGoList{$trsptId}));
		@go = ();
		@go = split(/,/, $trsptToGoList{$trsptId});
		foreach $go(@go){
			if(index($goList, $go)<0){
				$goList .= $go . ",";
			}
		}
	}
	if($goList ne ""){
		print WW join("\t", $geneId, substr($goList, 0, length($goList) - 1)) . "\n";
	}
}
close WW;
