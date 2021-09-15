#!/usr/bin/perl
use strict;
use Getopt::Long; 
if($#ARGV<0){
	print "This script is used to draw venn to describe final annotation origin,\n" .
		"and then the venn will be put on website\n\n";
	print "perl $0 \\\n" .
		"--Rscript /opt/software/R/3.5.1-foss-2018a-X11-20180131/bin/Rscript \\\n" .
		"--xvfbRun  /usr/bin/xvfb-run \\\n" .
		"--drawVennRscript ~/software/customized/draw.vennDiagram.for.three.sets.R \\\n" .
		"--set1Label Ensembl \\\n" .
		"--set2Label Refseq \\\n" .
		"--set3Label RNAseq \\\n" .
		"--set1Size 23585 \\\n" .
		"--set2Size 51071 \\\n" .
		"--set3Size 65570 \\\n" .
		"--set12Size 20362 \\\n" .
		"--set13Size 0 \\\n" .
		"--set23Size 0 \\\n" .
		"--set123Size 0 \\\n" .
		"--tmpDir ./tmp \\\n" .
		"--outputPng ./btau.anno.venn.png \n\n";
		exit;
}
my ($Rscript, $xvfbRun, $drawVennRscript, $set1Label, $set2Label, $set3Label );
my ($set1Size, $set2Size, $set3Size, $set12Size, $set13Size, $set23Size, $set123Size, $outputPng, $tmpDir);

GetOptions(
	'Rscript=s'=>\$Rscript,
	'xvfbRun=s'=>\$xvfbRun, 
	'drawVennRscript=s'=>\$drawVennRscript,
	'set1Label=s'=>\$set1Label, 
	'set2Label=s'=>\$set2Label,
	'set3Label=s'=>\$set3Label,
	'set1Size=i'=>\$set1Size, 
	'set2Size=i'=>\$set2Size, 
	'set3Size=i'=>\$set3Size, 
	'set12Size=i'=>\$set12Size,
	'set13Size=i'=>\$set13Size,
	'set23Size=i'=>\$set23Size,
	'set123Size=i'=>\$set123Size,
	'tmpDir=s'=>\$tmpDir,
	'outputPng=s'=>\$outputPng,
);

system("mkdir $tmpDir");
open ONE, ">$tmpDir/one.txt";
open TWO, ">$tmpDir/two.txt";
open THREE, ">$tmpDir/three.txt";

my $i;
# output set1
for($i = 1; $i <= $set1Size; $i++){
	print ONE "1_" . $i . "\n";
}
for($i = 1; $i <= $set12Size; $i++){
	print ONE "12_" . $i . "\n";
}
for($i = 1; $i <= $set13Size; $i++){
	print ONE "13_" . $i . "\n";
}
for($i = 1; $i <= $set123Size; $i++){
	print ONE "123_" . $i . "\n";
}

# output set2
for($i = 1; $i <= $set2Size; $i++){
	print TWO "2_" . $i . "\n";
}
for($i = 1; $i <= $set12Size; $i++){
	print TWO "12_" . $i . "\n";
}
for($i = 1; $i <= $set23Size; $i++){
	print TWO "23_" . $i . "\n";
}
for($i = 1; $i <= $set123Size; $i++){
	print TWO "123_" . $i . "\n";
}

# output set3
for($i = 1; $i <= $set3Size; $i++){
	print THREE "3_" . $i . "\n";
}
for($i = 1; $i <= $set13Size; $i++){
	print THREE "13_" . $i . "\n";
}
for($i = 1; $i <= $set23Size; $i++){
	print THREE "23_" . $i . "\n";
}
for($i = 1; $i <= $set123Size; $i++){
	print THREE "123_" . $i . "\n";
}


close ONE;
close TWO;
close THREE;

my $cmd = $xvfbRun . " " . $Rscript . " " . $drawVennRscript . " " .
        " -A $tmpDir/one.txt" .
        " -a $set1Label" .
        " -B $tmpDir/two.txt" .
        " -b $set2Label" .
        " -C $tmpDir/three.txt" .
        " -c $set3Label" .
        " -i $outputPng" .
        " -t png";
system($cmd);

