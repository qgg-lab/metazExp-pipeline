#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--duplicationGeneGroup \\\n" .
		"--conservedAsGroup \\\n" .
		"--asFileList \\\n" .
		"--dpgPrefix \\\n" .
                "--outputDPG \n";
	exit;
}

my ($duplicationGeneGroup, $conservedAsGroup, $asFileList, $dpgPrefix, $outputDPG);

GetOptions(
        'duplicationGeneGroup=s'=>\$duplicationGeneGroup,
        'conservedAsGroup=s'=>\$conservedAsGroup,
	'asFileList=s'=>\$asFileList,
	'dpgPrefix=s'=>\$dpgPrefix,
	'outputDPG=s'=>\$outputDPG,
);

my (%asToGene, %geneToDpgg);

# 从asFile读入as，获得as->gene的映射关系
my ($line, @asFile, $asFile, @field);
@asFile = split(/,/, $asFileList);
foreach $asFile(@asFile){
	open FF, "<$asFile";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		$asToGene{$field[0]} = substr($field[1], 1, length($field[1]) - 2 );
	}	
	close FF;
}

# 从duplication gene group中获得gene->duplicationGeneGroup
my ($geneIdList, @geneId, $geneId,  $dpggId, $num1, $num2);
open FF, "<$duplicationGeneGroup";
# AthaDpgg00030   AT1G59406,AT1G58525,AT3G43550,AT1G59030,AT1G58725,AT3G43570     6       6
<FF>;
while($line=<FF>){
	chomp($line);
	($dpggId, $geneIdList, $num1, $num2) = split(/\t/, $line);
	@geneId = ();
	@geneId = split(/,/, $geneIdList);
	foreach $geneId(@geneId){
		$geneToDpgg{$geneId} = $dpggId;
	}
}
close FF;


open WW, ">$outputDPG";
open FF, "<$conservedAsGroup";
my ($oldConservedAsGroup, $newConservedAsGroup, $asIdList, @asId, $asId);
# 读取保守AS的group信息
# orthA3SS00000204        ATHAA3SS0000000642,ATHAA3SS0000003172
# 输出：
# dupGeneGroupId        dupGeneAsGroupId      AsIdList
# AthaDpgg00030         AthaDpggA3SS00000204    ATHAA3SS0000000642,ATHAA3SS0000003172
print WW join("\t", "dupGeneGroupId", "dupGeneAsGroupId", "AsIdList") . "\n";
while($line=<FF>){
#	print $line;
#	<STDIN>;
	chomp($line);
	($oldConservedAsGroup, $asIdList) = split(/\t/, $line);

	@asId = ();
	@asId = split(/,/, $asIdList);
	$asId = $asId[0];
#	print "asId:" . $asId;
#	<STDIN>;
	$geneId = $asToGene{$asId};
#	print "geneId:" . $geneId;
#	<STDIN>;
	$dpggId = $geneToDpgg{$geneId};
#	print "dpggId:" . $dpggId;
#	<STDIN>;
	# 修改asGroupId的格式
	$newConservedAsGroup = $dpgPrefix . substr($oldConservedAsGroup, 4);

	# 输出asGroup
	print WW join("\t", $dpggId, $newConservedAsGroup, $asIdList) . "\n";
}
close FF;
close WW;
