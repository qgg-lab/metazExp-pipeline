#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asTissuePsiFile ASxa \\\n" .
                "--Rscript /bin/Rscript \\\n" .
                "--psiBoxPlotRscript draw.tissue.psi.boxplot.R \\\n" .
                "--outputBoxplotDir psiBoxPlot\n";
	exit;
}

my ($asTissuePsiFile, $Rscript, $psiBoxPlotRscript, $outputBoxplotDir, $xvfb_run);
GetOptions(
        'asTissuePsiFile=s'=>\$asTissuePsiFile,
        'Rscript=s'=>\$Rscript,
        'psiBoxPlotRscript=s'=>\$psiBoxPlotRscript,
        'outputBoxplotDir=s'=>\$outputBoxplotDir,
);



system("mkdir -p " . $outputBoxplotDir);

my ($asId, @tissue, $tissue, $speciesList, @species, $species, $psi, @psi,  $line, $cmd, $tissueNum, $expNum, $tissueName, $speciesName, $speciesNum, $maxSpeciesNum);
open FF, "<$asTissuePsiFile";

#ORTHSESE0000030350      liver:Bos taurus,0.942557918121231,0.946185688941455,0.939941690962099,0.912221144519884,0.956480153036824,0.950356555128908,0.89179548156956,0.971356227251147,0.880724876441516,0.882479552303056,0.926362896663954,0.902896995708155#Ovis aries,0.97507863537382,0.98147148767764,0.975446960667461,0.986344955587962,0.97571613815867,0.980377214707563,0.989910862964051,0.936850738108256,0.964287483282976,0.925277230801119,0.978880459298749#  muscle:Equus caballus,0.977105368857946,0.891273247496423,0.935786706721742,0.959051724137931,0.964121674321865,0.949801956187859,0.934695436318503,0.910306845003934,0.929996618194116,0.981733145075891,0.95816001957426,0.979719800137161,0.912352639671963,0.959906213364595#Bos taurus,0.95379075815163,0.902896995708155,0.970795892169448,0.96396603167317,0.958624324680459,0.956480153036824,0.884020618556701,0.936744560838034,0.867786705624544,0.976308253059099,0.985871759043627,0.93546086646461,0.945954016124216,0.940349544072948,0.979293062516486,0.945954016124216,0.96160409556314,0.977480821578817,0.907794192562405,0.949059052563271,0.944307692307692#Ovis aries,0.973569969356486,0.954095866091808,0.955093736375527,0.932835820895522,0.964898330114734,0.93331175137585,0.968384279475982,0.889521437466643,0.972240847525815,0.929996618194116,0.929996618194116,0.976612695965048,0.925257266654631,0.980377214707563,0.921955403087478,0.95311695002576,0.950646861523718,0.938967136150235,0.961177794448612# 

while($line=<FF>){
	chomp($line);
	@tissue = ();
	@tissue = split(/\t/, $line);
	$tissueNum = $#tissue;
	$asId = shift(@tissue);

	$maxSpeciesNum = 0;	
	open WW, ">$outputBoxplotDir" . "/" . $asId . ".psi.tsv";
	print WW join("\t", "Tissue", "Species", "Psi") . "\n";
	foreach $tissue(@tissue){

		($tissueName, $speciesList) = split(/:/, $tissue);;
		$tissueName = "unknown" if($tissueName eq ".");
		
		@species = ();
		@species = split(/#/, $speciesList);
		$speciesNum = $#species + 1;
		if($speciesNum > $maxSpeciesNum){
			$maxSpeciesNum = $speciesNum;
		}
		foreach $species(@species){
			@psi = ();
			@psi = split(/,/, $species);
			$speciesName = shift(@psi);
			$expNum = $#psi + 1;
			foreach $psi(@psi){
				print WW join("\t", $tissueName, $speciesName, $psi) . "\n";
			}
		}
	}
	close WW;
	
	my $fontSize = 14;
	my $ncol = 0;
	# draw boxplot
	if($tissueNum > 36){
		$fontSize = 12;
		$ncol = 7;
	}elsif($tissueNum > 25){
		$fontSize = 12;
		$ncol = 6;
	}elsif($tissueNum > 16){
		$fontSize = 14;
		$ncol = 5;
	}elsif($tissueNum > 9){
		$fontSize = 14;
		$ncol = 4;
	}elsif($tissueNum > 4){
		$fontSize = 16;
		$ncol = 3;
	}elsif($tissueNum > 1){
		$fontSize = 16;
		$ncol = 2;
	}else{
		$fontSize = 16;
		$ncol = 1;
	}

	$cmd = $Rscript . " " . $psiBoxPlotRscript . " -p " . $outputBoxplotDir . "/" . $asId . ".psi.tsv" . " -i " . $outputBoxplotDir . "/" . $asId . " -f png" . " -n " . $tissueNum  . " -e " . $maxSpeciesNum . " -c " . $ncol . " -s " . $fontSize;
	system($cmd);

	system("rm -rf " . $outputBoxplotDir . "/" . $asId . ".psi.tsv");

}
close FF;
