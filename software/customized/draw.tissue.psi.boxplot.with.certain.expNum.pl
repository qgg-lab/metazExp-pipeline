#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asTissuePsiFile ASxa \\\n" .
		"--minExpNum 3 \\\n" .
                "--Rscript /bin/Rscript \\\n" .
                "--psiBoxPlotRscript draw.tissue.psi.boxplot.R \\\n" .
                "--outputBoxplotDir psiBoxPlot\n";
	exit;
}

my ($asTissuePsiFile, $Rscript, $psiBoxPlotRscript, $outputBoxplotDir, $xvfb_run, $minExpNum);
GetOptions(
        'asTissuePsiFile=s'=>\$asTissuePsiFile,
        'minExpNum=s'=>\$minExpNum,
        'Rscript=s'=>\$Rscript,
        'psiBoxPlotRscript=s'=>\$psiBoxPlotRscript,
        'outputBoxplotDir=s'=>\$outputBoxplotDir,
);



system("mkdir -p " . $outputBoxplotDir);

my ($asId, @tissuePsi, $tissuePsi, $tissue, $psi, @psi,  $line, $cmd, $tissueNum, $expNum);
open FF, "<$asTissuePsiFile";
while($line=<FF>){
	chomp($line);
	@tissuePsi = ();
	@tissuePsi = split(/\t/, $line);
	$tissueNum = $#tissuePsi;
	$asId = shift(@tissuePsi);

#	next if($asId ne "BTAUMXE0000010058");

	# generate psi data file
	open WW, ">$outputBoxplotDir" . "/" . $asId . ".psi.tsv";
	print WW join("\t", "Tissue", "Psi") . "\n";
	foreach $tissuePsi(@tissuePsi){
		@psi = ();
		@psi = split(/,/, $tissuePsi);
		$expNum = $#psi;
		$tissue = shift(@psi);
		$tissue = "unknown" if($tissue eq ".");
		if($tissue=~/'/){
			$tissue=~s/'//g;
		}
		if($tissue=~/"/){
			$tissue=~s/"//g;
		}
		foreach $psi(@psi){
			print WW join("\t", $tissue . "(" . $expNum . ")", $psi) . "\n";
		}
	}
	close WW;

	# draw boxplot
	$cmd = $Rscript . " " . $psiBoxPlotRscript . " -p " . $outputBoxplotDir . "/" . $asId . ".psi.tsv" . " -i " . $outputBoxplotDir . "/" . $asId . " -f png" . " -n " . $tissueNum . " -e " . $minExpNum;
	system($cmd);

	system("rm -rf " . $outputBoxplotDir . "/" . $asId . ".psi.tsv");
}
close FF;
