#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtfFile \\\n" .
                "--beforeAsPfamListFile \\\n" .
                "--afterAsPfamListFile \\\n" .
                "--beforeAsGoListFile \\\n" .
                "--afterAsGoListFile \\\n" .
		"--outputGenePfamGoAnno \n";
	exit;
}

my ($gtfFile, $beforeAsPfamListFile, $afterAsPfamListFile, $beforeAsGoListFile, $afterAsGoListFile, $outputGenePfamGoAnno);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
        'beforeAsPfamListFile=s'=>\$beforeAsPfamListFile,
        'afterAsPfamListFile=s'=>\$afterAsPfamListFile,
        'beforeAsGoListFile=s'=>\$beforeAsGoListFile,
        'afterAsGoListFile=s'=>\$afterAsGoListFile,
        'outputGenePfamGoAnno=s'=>\$outputGenePfamGoAnno,
);
my ($line, @field, $geneId, $trsptId);
my (%trsptIdToGeneId);

# 读取gtf，获得gene->trspt的映射关系
open FF, "<$gtfFile";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");
	($geneId, $trsptId) = ("", "");
	&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);
	$trsptIdToGeneId{$trsptId} = $geneId;
}
close FF;


my (%geneAnno, $geneAnnoHref);
$geneAnnoHref = \%geneAnno;
my ($pfamPosList, @pfamPos, $pfamPos, $pfamId, $pos);
# 读取AS作用前transcript的pfam注释
open FF, "<$beforeAsPfamListFile";
# AT2G47000.1     PF00664[449,1267]|PF00005[1472,1915]|PF00664[2420,3241]|PF00005[3443,3892]
while($line=<FF>){
	chomp($line);
	($geneId, $trsptId, $pfamPosList) = ("", "", "");
	($trsptId, $pfamPosList) = split(/\t/, $line);
	$geneId = $trsptIdToGeneId{$trsptId};

	@pfamPos = ();
	@pfamPos = split(/\|/, $pfamPosList);
	foreach $pfamPos(@pfamPos){
		($pfamId, $pos) = ("", "");
		($pfamId, $pos) = split(/\[/, $pfamPos);
		if(not(exists($geneAnnoHref->{$geneId}->{"pfamList"}))){
			$geneAnnoHref->{$geneId}->{"pfamList"} = $pfamId;
		}elsif(index($geneAnnoHref->{$geneId}->{"pfamList"}, $pfamId)<0){
			$geneAnnoHref->{$geneId}->{"pfamList"} .= "," . $pfamId;
		}
	}
}
close FF;

# 读取AS作用后transcript的pfam注释
my ($asTrsptId, $asId);
open FF, "<$afterAsPfamListFile";
#AT5G06600.2XXXATHASE0000009673  PF00917[289,651]|PF00443[703,1383]|PF12436[1906,2658]|PF14533[2686,3318]
while($line=<FF>){
	chomp($line);
	($geneId, $trsptId, $pfamPosList) = ("", "", "");
	($asTrsptId, $pfamPosList) = split(/\t/, $line);
	($trsptId, $asId) = ("", "");
	($trsptId, $asId) = split(/XXX/, $asTrsptId);
	$geneId = $trsptIdToGeneId{$trsptId};

	@pfamPos = ();
	@pfamPos = split(/\|/, $pfamPosList);
	foreach $pfamPos(@pfamPos){
		($pfamId, $pos) = ("", "");
		($pfamId, $pos) = split(/\[/, $pfamPos);
		if(not(exists($geneAnnoHref->{$geneId}->{"pfamList"}))){
			$geneAnnoHref->{$geneId}->{"pfamList"} = $pfamId;
		}elsif(index($geneAnnoHref->{$geneId}->{"pfamList"}, $pfamId)<0){
			$geneAnnoHref->{$geneId}->{"pfamList"} .= "," . $pfamId;
		}
	}
}
close FF;


# 处理GO
my ($goList, @go, $go);
# 读取AS作用前transcript的go注释
open FF, "<$beforeAsGoListFile";
# AT5G20320.1     GO:0003676|GO:0004525|GO:0005515|GO:0005524|GO:0006396|GO:0016891
while($line=<FF>){
	chomp($line);
	($geneId, $trsptId, $goList) = ("", "", "");
	($trsptId, $goList) = split(/\t/, $line);
	$geneId = $trsptIdToGeneId{$trsptId};

	@go = ();
	@go = split(/\|/, $goList);
	foreach $go(@go){
		if(not(exists($geneAnnoHref->{$geneId}->{"goList"}))){
			$geneAnnoHref->{$geneId}->{"goList"} = $go;
		}elsif(index($geneAnnoHref->{$geneId}->{"goList"}, $go)<0){
			$geneAnnoHref->{$geneId}->{"goList"} .= "," . $go;
		}
	}
}
close FF;

# 读取AS作用后transcript的go注释
my ($asTrsptId, $asId);
open FF, "<$afterAsGoListFile";
# AT1G16710.7XXXATHARI0000004237  GO:0003712|GO:0004402|GO:0005634|GO:0006355|GO:0008270|GO:0016573
while($line=<FF>){
	chomp($line);
	($geneId, $trsptId, $goList) = ("", "", "");
	($asTrsptId, $goList) = split(/\t/, $line);
	($trsptId, $asId) = ("", "");
	($trsptId, $asId) = split(/XXX/, $asTrsptId);
	$geneId = $trsptIdToGeneId{$trsptId};

	@go = ();
	@go = split(/\|/, $goList);
	foreach $go(@go){
		if(not(exists($geneAnnoHref->{$geneId}->{"goList"}))){
			$geneAnnoHref->{$geneId}->{"goList"} = $go;
		}elsif(index($geneAnnoHref->{$geneId}->{"goList"}, $go)<0){
			$geneAnnoHref->{$geneId}->{"goList"} .= "," . $go;
		}
	}
}
close FF;


# 将基因的GO和pfam注释结果输出
my @geneId = keys(%geneAnno);
open WW, ">$outputGenePfamGoAnno";
print WW join("\t", "geneId", "pfamList", "GOList") . "\n";
foreach $geneId(@geneId){
	print WW join("\t", $geneId, $geneAnnoHref->{$geneId}->{"pfamList"}, $geneAnnoHref->{$geneId}->{"goList"}) . "\n";
}
close WW;

# 获得geneId和transcriptId
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
