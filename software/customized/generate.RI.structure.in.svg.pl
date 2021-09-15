#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV<0){
	print "\nperl $0 \\\n" .
		"--strand + \\\n" .
		"--riExonStart_0base 6613000 \\\n" .
		"--riExonEnd 6614211  \\\n" .
		"--upstreamES 6613000 \\\n" .
		"--upstreamEE 6613085 \\\n" .
		"--downstreamES 6613400 \\\n" .
		"--downstreamEE 6614211  \\\n" .
		"--containerWidth 1000 \\\n" .
		"--containerHeight 200 \\\n" .
		"--middleLineWidth 4 \\\n" .
		"--middleLineColor 000 \\\n" .
		"--upperLineWidth 0.5 \\\n" .
		"--upperLineColor 00F \\\n" .
		"--belowLineWidth 0.5 \\\n" .
		"--belowLineColor 906 \\\n" .
		"--firstExonFillColor  06C \\\n" .
		"--secondExonFillColor  F90 \\\n" .
		"--thirdExonFillColor  06C \\\n" .
		"--exonBorderColor 000 \\\n" .
		"--exonBorderWidth 0.5 \\\n" .
		"--inputVarFile var.tsv \\\n" .
		"--varOvalBorderColor 000 \\\n" .
		"--varOvalBorderWidth 0.5 \\\n" .
		"--varOvalXLen 3 \\\n" .
		"--varOvalYLen 10 \\\n" .
		"--varSnvColor F00 \\\n" .
		"--varInsertColor 090 \\\n" .
		"--varDeleteColor 000 \\\n" .
		"--varSubstitColor FF0 \\\n" .
		"--outputSvgFile BtusSE000012345.svg\n\n";
	exit;
}

my ($strand, $riExonStart_0base, $riExonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE, $outputSvgFile);
my ($containerWidth, $containerHeight);
my ($firstExonFillColor, $secondExonFillColor, $thirdExonFillColor);
my ($exonBorderColor, $exonBorderWidth);
my ($middleLineWidth, $upperLineWidth, $belowLineWidth);
my ($middleLineColor, $upperLineColor, $belowLineColor);
my ($inputVarFile, $varOvalBorderColor, $varOvalBorderWidth, $varOvalXLen, $varOvalYLen, $varSnvColor,
 $varInsertColor, $varDeleteColor, $varSubstitColor);

GetOptions(
	'strand=s'=>\$strand,
	'riExonStart_0base=s'=>\$riExonStart_0base,
	'riExonEnd=s'=>\$riExonEnd,
	'upstreamES=s'=>\$upstreamES,
	'upstreamEE=s'=>\$upstreamEE,
	'downstreamES=s'=>\$downstreamES,
	'downstreamEE=s'=>\$downstreamEE,
	'containerWidth=s'=>\$containerWidth,
	'containerHeight=s'=>\$containerHeight,
	'middleLineWidth=s'=>\$middleLineWidth,
	'middleLineColor=s'=>\$middleLineColor,
	'upperLineWidth=s'=>\$upperLineWidth,
	'upperLineColor=s'=>\$upperLineColor,
	'belowLineWidth=s'=>\$belowLineWidth,
	'belowLineColor=s'=>\$belowLineColor,
	'firstExonFillColor=s'=>\$firstExonFillColor,
	'secondExonFillColor=s'=>\$secondExonFillColor,
	'thirdExonFillColor=s'=>\$thirdExonFillColor,
	'exonBorderColor=s'=>\$exonBorderColor,
	'exonBorderWidth=s'=>\$exonBorderWidth,
	'inputVarFile=s'=>\$inputVarFile,
	'varOvalBorderColor=s'=>\$varOvalBorderColor,
	'varOvalBorderWidth=s'=>\$varOvalBorderWidth,
	'varOvalXLen=s'=>\$varOvalXLen,
	'varOvalYLen=s'=>\$varOvalYLen,
	'varSnvColor=s'=>\$varSnvColor,
	'varInsertColor=s'=>\$varInsertColor,
	'varDeleteColor=s'=>\$varDeleteColor,
	'varSubstitColor=s'=>\$varSubstitColor,
	'outputSvgFile=s'=>\$outputSvgFile,
);

# declare

my ($contextWidth, $contextHeight, $paddingLeft, $paddingTop);
my ($contextHeightUnit, $intronContextY, $flagX, $flagY);
$contextWidth = 0.9 * $containerWidth;
$contextHeight = 0.9 * $containerHeight;
$paddingLeft = 0.05 * $containerWidth;
$paddingTop = 0.05 * $containerHeight;
$contextHeightUnit = $contextHeight/10;
$intronContextY = $contextHeightUnit * 5 + $paddingTop;

open WW, ">$outputSvgFile";
# calculate
print WW "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" height=\"" . $containerHeight . "\">\n";
print WW "<rect x=\'0\' y=\'0\' width=\'200\' height=\'25\' style=\"fill:#FFF; stroke:#FFF;stroke-width:0.5\"/>\n";
print WW "<text x=\"15\" y=\"20\" fill=\"red\">: SNV</text>\n";
print WW "<ellipse cx=\"10\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "5" . "\" style=\"fill:#F00; stroke:#F00; stroke-width:0.5\"/>\n";
print WW "<text x=\"75\" y=\"20\" fill=\"black\">: deletion</text>\n";
print WW "<ellipse cx=\"70\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "10" . "\" style=\"fill:#000; stroke:#000; stroke-width:0.5\"/>\n";
print WW "<text x=\"165\" y=\"20\" fill=\"green\">: insertion</text>\n";
print WW "<ellipse cx=\"160\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "20" . "\" style=\"fill:#090; stroke:#090; stroke-width:0.5\"/>\n";

my ($asGenomeSpan);
if($strand eq "+"){
	$asGenomeSpan = $riExonEnd - $riExonStart_0base;
	$flagX = $paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">5\'</text>\n";
	$flagX = $contextWidth + 3*$paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">3\'</text>\n";
}else{
	$asGenomeSpan = $riExonEnd - $riExonStart_0base;
	$flagX = $paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">3\'</text>\n";
	$flagX = $contextWidth + 3*$paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">5\'</text>\n";
}


#############################
# draw left external Intron
############################
print WW "<line x1=\"0\" y1=\"" . $intronContextY . "\" x2=\"" . $paddingLeft . "\" y2=\"" . $intronContextY . "\" style=\"stroke-dasharray:5;stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . ";\"/>\n";

#############################
# draw first exon
#############################
my ($firstExonGenomeSpan, $firstExonContextWidth, $firstExonHeight, $firstExonStartX, $firstExonStartY);
$firstExonGenomeSpan = $upstreamEE - $upstreamES;
$firstExonContextWidth = $contextWidth * ($firstExonGenomeSpan/$asGenomeSpan);
$firstExonHeight = $contextHeightUnit * 2;
$firstExonStartX = $paddingLeft;
$firstExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $firstExonStartX . "\' y=\'" . $firstExonStartY . "\' width=\'" . $firstExonContextWidth . "\' height=\'" . $firstExonHeight . "\' style=\"fill:#" . $firstExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#############################
# draw second exon
#############################
my ($secondExonGenomeSpan, $secondExonContextWidth, $secondExonHeight, $secondExonStartX, $secondExonStartY);
$secondExonGenomeSpan = $downstreamES - ($upstreamEE + 1) + 1;
$secondExonContextWidth = $contextWidth * ($secondExonGenomeSpan/$asGenomeSpan);
$secondExonHeight = $contextHeightUnit * 2;
$secondExonStartX = $paddingLeft + $firstExonContextWidth;
$secondExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $secondExonStartX . "\' y=\'" . $secondExonStartY . "\' width=\'" . $secondExonContextWidth . "\' height=\'" . $secondExonHeight . "\' style=\"fill:#" . $secondExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#############################
# draw third exon
#############################
my ($thirdExonGenomeSpan, $thirdExonContextWidth, $thirdExonHeight, $thirdExonStartX, $thirdExonStartY);
$thirdExonGenomeSpan = $downstreamEE - $downstreamES;
$thirdExonContextWidth = $contextWidth * ($thirdExonGenomeSpan/$asGenomeSpan);
$thirdExonHeight = $contextHeightUnit * 2;
$thirdExonStartX = $paddingLeft + $firstExonContextWidth + $secondExonContextWidth;
$thirdExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $thirdExonStartX . "\' y=\'" . $thirdExonStartY . "\' width=\'" . $thirdExonContextWidth . "\' height=\'" . $thirdExonHeight . "\' style=\"fill:#" . $thirdExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#############################
# draw right external intron
############################
my $externalIntronStartX = $thirdExonStartX + $thirdExonContextWidth;
my $externalIntronEndX = $containerWidth;
print WW "<line x1=\"" . $externalIntronStartX . "\" y1=\"" . $intronContextY . "\" x2=\"" . $externalIntronEndX . "\" y2=\"" . $intronContextY . "\" style=\"stroke-dasharray:5;stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . ";\"/>\n";

###########################
# draw upperLine1
###########################
my ($upperLine1StartX, $upperLine1StartY, $upperLine1EndX, $upperLine1EndY);
$upperLine1StartX = $secondExonStartX;
$upperLine1StartY = $contextHeightUnit * 4 + $paddingTop;
$upperLine1EndX = $upperLine1StartX + ($secondExonContextWidth)/2;
$upperLine1EndY = $paddingTop;
print WW "<line x1=\"" . $upperLine1StartX . "\" y1=\"" . $upperLine1StartY . "\" x2=\"" . $upperLine1EndX . "\" y2=\"" . $upperLine1EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";

##########################
# draw upperLine2
##########################
my ($upperLine2StartX, $upperLine2StartY, $upperLine2EndX, $upperLine2EndY);
$upperLine2StartX = $upperLine1EndX;
$upperLine2StartY = $upperLine1EndY;
$upperLine2EndX = $thirdExonStartX;
$upperLine2EndY = $contextHeightUnit * 4 + $paddingTop;
print WW "<line x1=\"" . $upperLine2StartX . "\" y1=\"" . $upperLine2StartY . "\" x2=\"" . $upperLine2EndX . "\" y2=\"" . $upperLine2EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";


#######################
#
# draw variant position
#
#######################

my ($line, @fields, $varPositionInAS, $varPositionInSvgX, $varPositionInSvgY, $varFillColor, $genomeBeg);
open FF, "<$inputVarFile";

if($strand eq "+"){
	$genomeBeg = $riExonStart_0base;
}else{
	$genomeBeg = $riExonStart_0base;
}
while($line=<FF>){
	@fields = ();
	@fields = split(/\t/, $line);
	$varPositionInAS = $fields[3] - $genomeBeg;
	$varPositionInSvgX = $paddingLeft + $contextWidth * $varPositionInAS/$asGenomeSpan;
	$varPositionInSvgY = $paddingTop + $contextHeightUnit * 5;
	if($fields[6]=~/SNV/){
		$varFillColor = $varSnvColor;
		print WW "<ellipse cx=\"" . $varPositionInSvgX . "\" cy=\"" . $varPositionInSvgY . "\" rx=\"" . "0.5" . "\" ry=\"" . "5" . "\" style=\"fill:#" . $varFillColor . "; stroke:#" . $varOvalBorderColor . "; stroke-width:" . $varOvalBorderWidth . "\"/>\n";
	}elsif($fields[6]=~/deletion/){
		$varFillColor = $varDeleteColor;
		print WW "<ellipse cx=\"" . $varPositionInSvgX . "\" cy=\"" . $varPositionInSvgY . "\" rx=\"" . "0.5" . "\" ry=\"" . "10" . "\" style=\"fill:#" . $varFillColor . "; stroke:#" . $varOvalBorderColor . "; stroke-width:" . $varOvalBorderWidth . "\"/>\n";
	}elsif($fields[6]=~/insertion/){
		$varFillColor = $varInsertColor;
		print WW "<ellipse cx=\"" . $varPositionInSvgX . "\" cy=\"" . $varPositionInSvgY . "\" rx=\"" . "0.5" . "\" ry=\"" . "20" . "\" style=\"fill:#" . $varFillColor . "; stroke:#" . $varOvalBorderColor . "; stroke-width:" . $varOvalBorderWidth . "\"/>\n";
	}elsif($fields[6]=~/substitution/){
		$varFillColor = $varSubstitColor;
	}
	# print WW "<ellipse cx=\"" . $varPositionInSvgX . "\" cy=\"" . $varPositionInSvgY . "\" rx=\"" . $varOvalXLen . "\" ry=\"" . $varOvalYLen . "\" style=\"fill:#" . $varFillColor . "; stroke:#" . $varOvalBorderColor . "; stroke-width:" . $varOvalBorderWidth . "\"/>\n";
}
close FF;
print WW "</svg>\n";
close WW;