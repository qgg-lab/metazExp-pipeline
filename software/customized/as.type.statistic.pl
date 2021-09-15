#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdList \\\n" .
                "--baseDir \\\n" .
                "--outputAsTypePerTsv\n";
	exit;
}

my ($taxonIdList, $baseDir, $outputAsTypePerTsv, $asNumInExperimentTsv);

GetOptions(
        'taxonIdList=s'=>\$taxonIdList,
        'baseDir=s'=>\$baseDir,
        'outputAsTypePerTsv=s'=>\$outputAsTypePerTsv,
	'outputAsNumInExperimentTsv=s'=>\$asNumInExperimentTsv,
);



my (@taxonId, $taxonId, @name, @value, %tmpHash, $line, $nameList, $valueList);
my ($experimentInfoGz);

open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;

open ASNUMINEXPT, ">$asNumInExperimentTsv";
print ASNUMINEXPT join("\t", "TaxonId", "ExperimentId", "Tissue", "mapped Bases(GB)", "AS Num", "A3SS(%)", "A5SS(%)", "MXE(%)", "RI(%)", "SE(%)") . "\n";
open WW, ">$outputAsTypePerTsv";
print WW join("\t", "percentage", "type") . "\n";
foreach $taxonId(@taxonId){
	chomp($taxonId);

	$experimentInfoGz="$baseDir/999999-make-mysqlTsv/$taxonId/015-experimentTable/$taxonId.curated.experiment.mysql.tsv.gz";
	`cat $experimentInfoGz | gunzip > experiment.tsv`;

	open FF, "<experiment.tsv";
	# ecotype, cultivar, genotype, tissue, subTissue, development, tissueGroup, treatment, treatmentGroup, experimentId, studyId, dataSource, totalBasedVolume, layout, spotNum, readNum, splotLen, readLen, assemblyTag, runIdList, phenotype, alignPercentage, mappedBaseVolume, mappedReadNum, detectedReadLen, libraryType, phredScore, fileterTrsptMinTrsptReadCov, filterTrsptMinExonReadCov, asNum, a3ssNum, a5ssNum, mxeNum, riNum, seNum___NA, NA, NA, seedling, NA, NA, NA, NA, NA, DRX014494, DRP001761, NA, 2.50, PAIRED, 12.10, 24.30, 202, 101, N, DRR016125, NA, 98.51, 2.46, 23.94, 100, RF, 33, -1, -1, 14673, 3126, 1958, 23, 8299, 1267 
	while($line=<FF>){
		chomp($line);
		($nameList, $valueList) = split(/___/, $line);
		@name = split(/, /, $nameList);
		@value = split(/, /, $valueList);
		for(my $i=0; $i<=$#name; $i++){
			$tmpHash{$name[$i]} = $value[$i];
		}
		print WW join("\t", sprintf("%.2f", $tmpHash{"a3ssNum"}/$tmpHash{"asNum"}*100), "A3SS") . "\n";
		print WW join("\t", sprintf("%.2f", $tmpHash{"a5ssNum"}/$tmpHash{"asNum"}*100), "A5SS") . "\n";
		print WW join("\t", sprintf("%.2f", $tmpHash{"mxeNum"}/$tmpHash{"asNum"}*100), "MXE") . "\n";
		print WW join("\t", sprintf("%.2f", $tmpHash{"riNum"}/$tmpHash{"asNum"}*100), "RI") . "\n";
		print WW join("\t", sprintf("%.2f", $tmpHash{"seNum"}/$tmpHash{"asNum"}*100), "SE") . "\n";

		print ASNUMINEXPT join("\t", $taxonId, $tmpHash{"experimentId"}, $tmpHash{"tissue"}, $tmpHash{"mappedBaseVolume"}, $tmpHash{"asNum"}, sprintf("%.2f", $tmpHash{"a3ssNum"}/$tmpHash{"asNum"}*100), sprintf("%.2f", $tmpHash{"a5ssNum"}/$tmpHash{"asNum"}*100), sprintf("%.2f", $tmpHash{"mxeNum"}/$tmpHash{"asNum"}*100), sprintf("%.2f", $tmpHash{"riNum"}/$tmpHash{"asNum"}*100), sprintf("%.2f", $tmpHash{"seNum"}/$tmpHash{"asNum"}*100)) . "\n";
	}
}

close WW;
close ASNUMINEXPT;
