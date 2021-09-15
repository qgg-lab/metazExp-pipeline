#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV<0){
	print "\nperl $0 \\\n" .
		"--strand + \\\n" .
		"--1stExonStart_0base 52389318 \\\n" .
		"--1stExonEnd 52389342  \\\n" .
		"--2ndExonStart_0base 52391234 \\\n" .
		"--2ndExonEnd 52391252 \\\n" .
		"--upstreamES 52380606 \\\n" .
		"--upstreamEE 52380780  \\\n" .
		"--downstreamES 52393370 \\\n" .
		"--downstreamEE 52393472 \\\n" .
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
		"--thirdExonFillColor  F90 \\\n" .
		"--fourthExonFillColor  06C \\\n" .
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

my ($strand, $FstExonStart_0base, $FstExonEnd, $ScndExonStart_0base, $ScndExonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE, $outputSvgFile);
my ($containerWidth, $containerHeight);
my ($firstExonFillColor, $secondExonFillColor, $thirdExonFillColor, $fourthExonFillColor);
my ($exonBorderColor, $exonBorderWidth);
my ($middleLineWidth, $upperLineWidth, $belowLineWidth);
my ($middleLineColor, $upperLineColor, $belowLineColor);
my ($inputVarFile, $varOvalBorderColor, $varOvalBorderWidth, $varOvalXLen, $varOvalYLen, $varSnvColor,
 $varInsertColor, $varDeleteColor, $varSubstitColor);

GetOptions(
	'strand=s'=>\$strand,
	'1stExonStart_0base=s'=>\$FstExonStart_0base,
	'1stExonEnd=s'=>\$FstExonEnd,
	'2ndExonStart_0base=s'=>\$ScndExonStart_0base,
	'2ndExonEnd=s'=>\$ScndExonEnd,
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
	'fourthExonFillColor=s'=>\$fourthExonFillColor,
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
my ($asGenomeSpan);

print WW "<rect x=\'0\' y=\'0\' width=\'200\' height=\'25\' style=\"fill:#FFF; stroke:#FFF;stroke-width:0.5\"/>\n";
print WW "<text x=\"15\" y=\"20\" fill=\"red\">: SNV</text>\n";
print WW "<ellipse cx=\"10\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "5" . "\" style=\"fill:#F00; stroke:#F00; stroke-width:0.5\"/>\n";
print WW "<text x=\"75\" y=\"20\" fill=\"black\">: deletion</text>\n";
print WW "<ellipse cx=\"70\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "10" . "\" style=\"fill:#000; stroke:#000; stroke-width:0.5\"/>\n";
print WW "<text x=\"165\" y=\"20\" fill=\"green\">: insertion</text>\n";
print WW "<ellipse cx=\"160\" cy=\"15\" rx=\"" . "0.5" . "\" ry=\"" . "20" . "\" style=\"fill:#090; stroke:#090; stroke-width:0.5\"/>\n";



if($strand eq "+"){
	$asGenomeSpan = $downstreamEE - $upstreamES;
	$flagX = $paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">5\'</text>\n";
	$flagX = $contextWidth + 3*$paddingLeft/2;
	$flagY = $intronContextY - 10;
	print WW "<text x=\"" . $flagX . "\" y=\"" . $flagY . "\" fill=\"black\">3\'</text>\n";
}else{
	$asGenomeSpan = $downstreamEE - $upstreamES;
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
#print WW $contextWidth . "\t" . $firstExonGenomeSpan . "\t" . $asGenomeSpan . "\n";
print WW "<rect x=\'". $firstExonStartX . "\' y=\'" . $firstExonStartY . "\' width=\'" . $firstExonContextWidth . "\' height=\'" . $firstExonHeight . "\' style=\"fill:#" . $firstExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#########################
# draw first intron
#########################
my ($firstIntronGenomeSpan, $firstIntronContextWidth, $firstIntronStartX, $firstIntronStartY, $firstIntronEndX, $firstIntronEndY);
$firstIntronGenomeSpan = $FstExonStart_0base - ($upstreamEE + 1) + 1;
$firstIntronContextWidth = $contextWidth * ($firstIntronGenomeSpan/$asGenomeSpan);
$firstIntronStartX = $paddingLeft + $firstExonContextWidth;
$firstIntronStartY = $intronContextY;
$firstIntronEndX = $firstIntronStartX + $firstIntronContextWidth;
$firstIntronEndY = $intronContextY;;
print WW "<line x1=\"" . $firstIntronStartX . "\" y1=\"" . $firstIntronStartY . "\" x2=\"" . $firstIntronEndX . "\" y2=\"" . $firstIntronEndY . "\" style=\"stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . "\"/>\n";

#############################
# draw second exon
#############################
my ($secondExonGenomeSpan, $secondExonContextWidth, $secondExonHeight, $secondExonStartX, $secondExonStartY);
$secondExonGenomeSpan = $FstExonEnd - $FstExonStart_0base;
$secondExonContextWidth = $contextWidth * ($secondExonGenomeSpan/$asGenomeSpan);
$secondExonHeight = $contextHeightUnit * 2;
$secondExonStartX = $firstIntronEndX;
$secondExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $secondExonStartX . "\' y=\'" . $secondExonStartY . "\' width=\'" . $secondExonContextWidth . "\' height=\'" . $secondExonHeight . "\' style=\"fill:#" . $secondExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#########################
# draw second intron
#########################
my ($secondIntronGenomeSpan, $secondIntronContextWidth, $secondIntronStartX, $secondIntronStartY, $secondIntronEndX, $secondIntronEndY);
$secondIntronGenomeSpan = $ScndExonStart_0base - ($FstExonEnd + 1) + 1;
$secondIntronContextWidth = $contextWidth * ($secondIntronGenomeSpan/$asGenomeSpan);
$secondIntronStartX = $paddingLeft + $firstExonContextWidth + $firstIntronContextWidth + $secondExonContextWidth;
$secondIntronStartY = $intronContextY;
$secondIntronEndX = $secondIntronStartX + $secondIntronContextWidth;
$secondIntronEndY = $intronContextY;;
print WW "<line x1=\"" . $secondIntronStartX . "\" y1=\"" . $secondIntronStartY . "\" x2=\"" . $secondIntronEndX . "\" y2=\"" . $secondIntronEndY . "\" style=\"stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . "\"/>\n";

#############################
# draw third exon
#############################
my ($thirdExonGenomeSpan, $thirdExonContextWidth, $thirdExonHeight, $thirdExonStartX, $thirdExonStartY);
$thirdExonGenomeSpan = $ScndExonEnd - $ScndExonStart_0base;
$thirdExonContextWidth = $contextWidth * ($thirdExonGenomeSpan/$asGenomeSpan);
$thirdExonHeight = $contextHeightUnit * 2;
$thirdExonStartX = $secondIntronEndX;
$thirdExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $thirdExonStartX . "\' y=\'" . $thirdExonStartY . "\' width=\'" . $thirdExonContextWidth . "\' height=\'" . $thirdExonHeight . "\' style=\"fill:#" . $thirdExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";

#########################
# draw third intron
#########################
my ($thirdIntronGenomeSpan, $thirdIntronContextWidth, $thirdIntronStartX, $thirdIntronStartY, $thirdIntronEndX, $thirdIntronEndY);
$thirdIntronGenomeSpan = $downstreamES - ($ScndExonEnd + 1) +1;
$thirdIntronContextWidth = $contextWidth * ($thirdIntronGenomeSpan/$asGenomeSpan);
$thirdIntronStartX = $paddingLeft + $firstExonContextWidth + $firstIntronContextWidth + $secondExonContextWidth + $secondIntronContextWidth + $thirdExonContextWidth;
$thirdIntronStartY = $intronContextY;
$thirdIntronEndX = $thirdIntronStartX + $thirdIntronContextWidth;
$thirdIntronEndY = $intronContextY;;
print WW "<line x1=\"" . $thirdIntronStartX . "\" y1=\"" . $thirdIntronStartY . "\" x2=\"" . $thirdIntronEndX . "\" y2=\"" . $thirdIntronEndY . "\" style=\"stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . "\"/>\n";

#############################
# draw fourth exon
#############################
my ($fourthExonGenomeSpan, $fourthExonContextWidth, $fourthExonHeight, $fourthExonStartX, $fourthExonStartY);
$fourthExonGenomeSpan = $downstreamEE - $downstreamES;
$fourthExonContextWidth = $contextWidth * ($fourthExonGenomeSpan/$asGenomeSpan);
$fourthExonHeight = $contextHeightUnit * 2;
$fourthExonStartX = $thirdIntronEndX;
$fourthExonStartY = $contextHeightUnit * 4 + $paddingTop;
print WW "<rect x=\'". $fourthExonStartX . "\' y=\'" . $fourthExonStartY . "\' width=\'" . $fourthExonContextWidth . "\' height=\'" . $fourthExonHeight . "\' style=\"fill:#" . $fourthExonFillColor . "; stroke:#" . $exonBorderColor . ";stroke-width:" . $exonBorderWidth . "\"/>\n";


#############################
# draw right external intron
############################
my $externalIntronStartX = $fourthExonStartX + $fourthExonContextWidth;
my $externalIntronEndX = $containerWidth;
print WW "<line x1=\"" . $externalIntronStartX . "\" y1=\"" . $intronContextY . "\" x2=\"" . $externalIntronEndX . "\" y2=\"" . $intronContextY . "\" style=\"stroke-dasharray:5;stroke:#" . $middleLineColor . ";stroke-width:" . $middleLineWidth . ";\"/>\n";


###########################
# draw upperLine1
###########################
my ($upperLine1StartX, $upperLine1StartY, $upperLine1EndX, $upperLine1EndY);
$upperLine1StartX = $firstIntronStartX;
$upperLine1StartY = $contextHeightUnit * 4 + $paddingTop;
$upperLine1EndX = $upperLine1StartX + ($firstIntronContextWidth)/2;
$upperLine1EndY = $paddingTop;
print WW "<line x1=\"" . $upperLine1StartX . "\" y1=\"" . $upperLine1StartY . "\" x2=\"" . $upperLine1EndX . "\" y2=\"" . $upperLine1EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";

##########################
# draw upperLine2
##########################
my ($upperLine2StartX, $upperLine2StartY, $upperLine2EndX, $upperLine2EndY);
$upperLine2StartX = $upperLine1EndX;
$upperLine2StartY = $upperLine1EndY;
$upperLine2EndX = $firstIntronEndX;
$upperLine2EndY = $contextHeightUnit * 4 + $paddingTop;
print WW "<line x1=\"" . $upperLine2StartX . "\" y1=\"" . $upperLine2StartY . "\" x2=\"" . $upperLine2EndX . "\" y2=\"" . $upperLine2EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";

###########################
# draw upperLine3
###########################
my ($upperLine3StartX, $upperLine3StartY, $upperLine3EndX, $upperLine3EndY);
$upperLine3StartX = $secondIntronStartX;
$upperLine3StartY = $contextHeightUnit * 4 + $paddingTop;
$upperLine3EndX = $upperLine3StartX + ($secondIntronContextWidth + $thirdExonContextWidth + $thirdIntronContextWidth)/2;
$upperLine3EndY = $paddingTop;
print WW "<line x1=\"" . $upperLine3StartX . "\" y1=\"" . $upperLine3StartY . "\" x2=\"" . $upperLine3EndX . "\" y2=\"" . $upperLine3EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";

##########################
# draw upperLine4
##########################
my ($upperLine4StartX, $upperLine4StartY, $upperLine4EndX, $upperLine4EndY);
$upperLine4StartX = $upperLine3EndX;
$upperLine4StartY = $upperLine3EndY;
$upperLine4EndX = $thirdIntronEndX;
$upperLine4EndY = $contextHeightUnit * 4 + $paddingTop;
print WW "<line x1=\"" . $upperLine4StartX . "\" y1=\"" . $upperLine4StartY . "\" x2=\"" . $upperLine4EndX . "\" y2=\"" . $upperLine4EndY . "\" style=\"stroke:#" . $upperLineColor . ";stroke-width:" . $upperLineWidth . "\"/>\n";

##########################
# draw belowLine1
##########################
my ($belowLine1StartX, $belowLine1StartY, $belowLine1EndX, $belowLine1EndY);
$belowLine1StartX = $firstIntronStartX;
$belowLine1StartY = $contextHeightUnit * 6 + $paddingTop;
$belowLine1EndX = $firstIntronStartX + ($firstIntronContextWidth + $secondExonContextWidth + $secondIntronContextWidth)/2;
$belowLine1EndY = $contextHeightUnit * 10 + $paddingTop;
print WW "<line x1=\"" . $belowLine1StartX . "\" y1=\"" . $belowLine1StartY . "\" x2=\"" . $belowLine1EndX . "\" y2=\"" . $belowLine1EndY . "\" style=\"stroke:#" . $belowLineColor . ";stroke-width:" . $belowLineWidth . "\"/>\n";

###########################
# draw belowLine2
###########################
my ($belowLine2StartX, $belowLine2StartY, $belowLine2EndX, $belowLine2EndY);
$belowLine2StartX = $belowLine1EndX;
$belowLine2StartY = $belowLine1EndY;
$belowLine2EndX = $thirdExonStartX;
$belowLine2EndY = $contextHeightUnit * 6 + $paddingTop;
print WW "<line x1=\"" . $belowLine2StartX . "\" y1=\"" . $belowLine2StartY . "\" x2=\"" . $belowLine2EndX . "\" y2=\"" . $belowLine2EndY . "\" style=\"stroke:#" . $belowLineColor . ";stroke-width:" . $belowLineWidth . "\"/>\n";

##########################
# draw belowLine3
##########################
my ($belowLine3StartX, $belowLine3StartY, $belowLine3EndX, $belowLine3EndY);
$belowLine3StartX = $thirdIntronStartX;
$belowLine3StartY = $contextHeightUnit * 6 + $paddingTop;
$belowLine3EndX = $thirdIntronStartX + ($thirdIntronContextWidth)/2;
$belowLine3EndY = $contextHeightUnit * 10 + $paddingTop;
print WW "<line x1=\"" . $belowLine3StartX . "\" y1=\"" . $belowLine3StartY . "\" x2=\"" . $belowLine3EndX . "\" y2=\"" . $belowLine3EndY . "\" style=\"stroke:#" . $belowLineColor . ";stroke-width:" . $belowLineWidth . "\"/>\n";

##########################
# draw belowLine4
##########################
my ($belowLine4StartX, $belowLine4StartY, $belowLine4EndX, $belowLine4EndY);
$belowLine4StartX = $belowLine3EndX;
$belowLine4StartY = $belowLine3EndY;
$belowLine4EndX = $fourthExonStartX;
$belowLine4EndY = $contextHeightUnit * 6 + $paddingTop;
print WW "<line x1=\"" . $belowLine4StartX . "\" y1=\"" . $belowLine4StartY . "\" x2=\"" . $belowLine4EndX . "\" y2=\"" . $belowLine4EndY . "\" style=\"stroke:#" . $belowLineColor . ";stroke-width:" . $belowLineWidth . "\"/>\n";

#######################
#
# draw variant position
#
#######################

my ($line, @fields, $varPositionInAS, $varPositionInSvgX, $varPositionInSvgY, $varFillColor, $genomeBeg);
open FF, "<$inputVarFile";

if($strand eq "+"){
	$genomeBeg = $upstreamES;
}else{
	$genomeBeg = $upstreamES;
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
