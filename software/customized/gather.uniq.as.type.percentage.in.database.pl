#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdList \\\n" .
                "--baseDir \\\n" .
		"--outputUniqAsTypePctgInEachPlant \\n" .
                "--outputUniqAsTypePctgInDatabase\n";
	exit;
}

my ($taxonIdList, $baseDir, 
$outputUniqAsTypePctgInEachPlant, $outputUniqAsTypePctgInDatabase, 
$outputOrthUniqAsTypePctgInDatabase, $outputOrthUniqAsTypePctgInEachPlant,
$outputTrsptUniqAsTypePctgInDatabase, $outputTrsptUniqAsTypePctgInEachPlant,
$outputAsInEachSpecies,
);

GetOptions(
        'taxonIdList=s'=>\$taxonIdList,
        'baseDir=s'=>\$baseDir,
        'outputUniqAsTypePctgInDatabase=s'=>\$outputUniqAsTypePctgInDatabase,
	'outputUniqAsTypePctgInEachPlant=s'=>\$outputUniqAsTypePctgInEachPlant,

        'outputTrsptUniqAsTypePctgInDatabase=s'=>\$outputTrsptUniqAsTypePctgInDatabase,
	'outputTrsptUniqAsTypePctgInEachPlant=s'=>\$outputTrsptUniqAsTypePctgInEachPlant,

        'outputOrthUniqAsTypePctgInDatabase=s'=>\$outputOrthUniqAsTypePctgInDatabase,
	'outputOrthUniqAsTypePctgInEachPlant=s'=>\$outputOrthUniqAsTypePctgInEachPlant,

	'outputAsInEachSpecies=s'=>\$outputAsInEachSpecies,
);



my (@taxonId, $taxonId, @name, @value, %tmpHash, $line, $nameList, $valueList);
my ($speciesMysqlTsv, %asTypeNum, %orthAsTypeNum, %basedOnTrsptAsTypeNum);

open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;

open ASNUMINSPECIES, ">$outputAsInEachSpecies";
open WW, ">$outputUniqAsTypePctgInEachPlant";
open TRSPT, ">$outputTrsptUniqAsTypePctgInEachPlant";
open ORTH, ">$outputOrthUniqAsTypePctgInEachPlant";

print WW join("\t", "taxonId", "percentage", "type") . "\n";
print TRSPT join("\t", "taxonId", "percentage", "type") . "\n";
print ORTH join("\t", "taxonId", "percentage", "type") . "\n";
foreach $taxonId(@taxonId){
	chomp($taxonId);

	$speciesMysqlTsv="$baseDir/999999-make-mysqlTsv/$taxonId/014-make-speciesTsv-to-insertMysql/$taxonId.species.tsv";

	open FF, "<$speciesMysqlTsv";
	#taxonId, speciesName, speciesAbbr, genomeAssemblyVersion, genomeAssemblyAddress, genomeAnnoVersion, genomeAnnoGtfAddress, improveGtfMinMappedVolum, improveGtfMinAligPer, improveGtfExptUnitVolume, improveGtfExptUnitTrsptMinCov, improveGtfExptUnitExonMinCov, improveGtfExptNumForImproveGtf, identifyAsMinMappedVolum, identifyAsMinAlignPercentage, identifyAsExptNum, origGtfGeneNum, origGtfTrsptNum, origAllAsNum, origA3ssNum, origA5ssNum, origMxeNum, origRiNum, origSeNum, improvedGtfTrsptNum, improvedAllAsNum, improvedA3ssNum, improvedA5ssNum, improvedMxeNum, improvedRiNum, improvedSeNum, finalAllAsNum, finalA3ssNum, finalA5ssNum, finalMxeNum, finalRiNum, finalSeNum, origInnerConserAllAsNum, origInnerConserA3ssNum, origInnerConserA5ssNum, origInnerConserMxeNum, origInnerConserRiNum, origInnerConserSeNum, origOrthConserAllAsNum, origOrthConserA3ssNum, origOrthConserA5ssNum, origOrthConserMxeNum, origOrthConserRiNum, origOrthConserSeNum, improvedInnerConserAllAsNum, improvedInnerConserA3ssNum, improvedInnerConserA5ssNum, improvedInnerConserMxeNum, improvedInnerConserRiNum, improvedInnerConserSeNum, improvedOrthConserAllAsNum, improvedOrthConserA3ssNum, improvedOrthConserA5ssNum, improvedOrthConserMxeNum, improvedOrthConserRiNum, improvedOrthConserSeNum, finalInnerConserAllAsNum, finalInnerConservedA3ssNum, finalInnerConservedA5ssNum, finalInnerConservedMxeNum, finalInnerConservedRiNum, finalInnerConservedSeNum, finalOrthConserAllAsNum, finalOrthConserA3ssNum, finalOrthConserA5ssNum, finalOrthConserMxeNum, finalOrthConserRiNum, finalOrthConserSeNum, improvedGtfExptList, identifyAsExptList, origAnnoName, speciesDesptText, annoDesptText, identifyAsDesptText, origExonNum, origSpliceNum, improvedExonNum, improvedSpliceNum___4081, Solanum lycopersicum, ...
	while($line=<FF>){
		chomp($line);
		($nameList, $valueList) = split(/___/, $line);
		@name = split(/, /, $nameList);
		@value = split(/, /, $valueList);
		for(my $i=0; $i<=$#name; $i++){
			$tmpHash{$name[$i]} = $value[$i];
		}

		print ASNUMINSPECIES join("\t", $taxonId, $tmpHash{"speciesName"}, $tmpHash{"finalAllAsNum"}, sprintf("%.2f", $tmpHash{"finalA3ssNum"}/$tmpHash{"finalAllAsNum"}*100), sprintf("%.2f", $tmpHash{"finalA5ssNum"}/$tmpHash{"finalAllAsNum"}*100), sprintf("%.2f", $tmpHash{"finalMxeNum"}/$tmpHash{"finalAllAsNum"}*100), sprintf("%.2f", $tmpHash{"finalRiNum"}/$tmpHash{"finalAllAsNum"}*100), sprintf("%.2f", $tmpHash{"finalSeNum"}/$tmpHash{"finalAllAsNum"}*100)) . "\n";

		print WW join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalA3ssNum"}/$tmpHash{"finalAllAsNum"}*100), "A3SS") . "\n";
		print WW join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalA5ssNum"}/$tmpHash{"finalAllAsNum"}*100), "A5SS") . "\n";
		print WW join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalMxeNum"}/$tmpHash{"finalAllAsNum"}*100), "MXE") . "\n";
		print WW join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalRiNum"}/$tmpHash{"finalAllAsNum"}*100), "RI") . "\n";
		print WW join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalSeNum"}/$tmpHash{"finalAllAsNum"}*100), "SE") . "\n";
		$asTypeNum{"A3SS"}+=$tmpHash{"finalA3ssNum"};
		$asTypeNum{"A5SS"}+=$tmpHash{"finalA5ssNum"};
		$asTypeNum{"MXE"}+=$tmpHash{"finalMxeNum"};
		$asTypeNum{"RI"}+=$tmpHash{"finalRiNum"};
		$asTypeNum{"SE"}+=$tmpHash{"finalSeNum"};



		print ORTH join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalOrthConserA3ssNum"}/$tmpHash{"finalOrthConserAllAsNum"}*100), "A3SS") . "\n";
		print ORTH join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalOrthConserA5ssNum"}/$tmpHash{"finalOrthConserAllAsNum"}*100), "A5SS") . "\n";
		print ORTH join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalOrthConserMxeNum"}/$tmpHash{"finalOrthConserAllAsNum"}*100), "MXE") . "\n";
		print ORTH join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalOrthConserRiNum"}/$tmpHash{"finalOrthConserAllAsNum"}*100), "RI") . "\n";
		print ORTH join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $tmpHash{"finalOrthConserSeNum"}/$tmpHash{"finalOrthConserAllAsNum"}*100), "SE") . "\n";
		$orthAsTypeNum{"A3SS"}+=$tmpHash{"finalOrthConserA3ssNum"};
		$orthAsTypeNum{"A5SS"}+=$tmpHash{"finalOrthConserA5ssNum"};
		$orthAsTypeNum{"MXE"}+=$tmpHash{"finalOrthConserMxeNum"};
		$orthAsTypeNum{"RI"}+=$tmpHash{"finalOrthConserRiNum"};
		$orthAsTypeNum{"SE"}+=$tmpHash{"finalOrthConserSeNum"};



		my $basedonTrsptAllNum = $tmpHash{"improvedAllAsNum"};
		my $basedonTrsptA3ssNum = $tmpHash{"improvedA3ssNum"};
		$basedOnTrsptAsTypeNum{"A3SS"} += $basedonTrsptA3ssNum;
		my $basedonTrsptA5ssNum = $tmpHash{"improvedA5ssNum"};
		$basedOnTrsptAsTypeNum{"A5SS"} += $basedonTrsptA5ssNum;
		my $basedonTrsptMxeNum = $tmpHash{"improvedMxeNum"};
		$basedOnTrsptAsTypeNum{"MXE"} += $basedonTrsptMxeNum;
		my $basedonTrsptRiNum = $tmpHash{"improvedRiNum"};
		$basedOnTrsptAsTypeNum{"RI"} += $basedonTrsptRiNum;
		my $basedonTrsptSeNum = $tmpHash{"improvedSeNum"};
		$basedOnTrsptAsTypeNum{"SE"} += $basedonTrsptSeNum;
		print TRSPT join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $basedonTrsptA3ssNum/$basedonTrsptAllNum*100), "A3SS") . "\n";
		print TRSPT join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $basedonTrsptA5ssNum/$basedonTrsptAllNum*100), "A5SS") . "\n";
		print TRSPT join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $basedonTrsptMxeNum/$basedonTrsptAllNum*100), "MXE") . "\n";
		print TRSPT join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $basedonTrsptRiNum/$basedonTrsptAllNum*100), "RI") . "\n";
		print TRSPT join("\t", $tmpHash{"taxonId"}, sprintf("%.2f", $basedonTrsptSeNum/$basedonTrsptAllNum*100), "SE") . "\n";
	}
}

close WW;
close ORTH;
close TRSPT;
close ASNUMINSPECIES;

my $totalAsNum = $asTypeNum{"A3SS"}+$asTypeNum{"A5SS"}+$asTypeNum{"MXE"}+$asTypeNum{"RI"}+$asTypeNum{"SE"};
open WW, ">$outputUniqAsTypePctgInDatabase";
print WW join("\t", "A3SS", "A5SS", "MXE", "RI", "SE") . "\n";
print WW join("\t", sprintf("%.2f", 100*$asTypeNum{"A3SS"}/$totalAsNum), sprintf("%.2f", 100*$asTypeNum{"A5SS"}/$totalAsNum), sprintf("%.2f", 100*$asTypeNum{"MXE"}/$totalAsNum), sprintf("%.2f", 100*$asTypeNum{"RI"}/$totalAsNum), sprintf("%.2f", 100*$asTypeNum{"SE"}/$totalAsNum)) . "\n";
close WW;

my $orthAsNum = $orthAsTypeNum{"A3SS"} + $orthAsTypeNum{"A5SS"} + $orthAsTypeNum{"MXE"} + $orthAsTypeNum{"RI"} + $orthAsTypeNum{"SE"};
open WW, ">$outputOrthUniqAsTypePctgInDatabase";
print WW join("\t", "A3SS", "A5SS", "MXE", "RI", "SE") . "\n";
print WW join("\t", sprintf("%.2f", 100*$orthAsTypeNum{"A3SS"}/$orthAsNum), sprintf("%.2f", 100*$orthAsTypeNum{"A5SS"}/$orthAsNum), sprintf("%.2f", 100*$orthAsTypeNum{"MXE"}/$orthAsNum), sprintf("%.2f", 100*$orthAsTypeNum{"RI"}/$orthAsNum), sprintf("%.2f", 100*$orthAsTypeNum{"SE"}/$orthAsNum)) . "\n";
close WW;

my $trsptAsNum = $basedOnTrsptAsTypeNum{"A3SS"} + $basedOnTrsptAsTypeNum{"A5SS"} + $basedOnTrsptAsTypeNum{"MXE"} + $basedOnTrsptAsTypeNum{"RI"} + $basedOnTrsptAsTypeNum{"SE"};
open WW, ">$outputTrsptUniqAsTypePctgInDatabase";
print WW join("\t", "A3SS", "A5SS", "MXE", "RI", "SE") . "\n";
print WW join("\t", sprintf("%.2f", $basedOnTrsptAsTypeNum{"A3SS"}/$trsptAsNum*100), sprintf("%.2f", $basedOnTrsptAsTypeNum{"A5SS"}/$trsptAsNum*100), sprintf("%.2f", $basedOnTrsptAsTypeNum{"MXE"}/$trsptAsNum*100), sprintf("%.2f", $basedOnTrsptAsTypeNum{"RI"}/$trsptAsNum*100), sprintf("%.2f", $basedOnTrsptAsTypeNum{"SE"}/$trsptAsNum*100));
close WW;

print "database中搜集的总AS数量为：$totalAsNum, ";
print "A3SS:" . $asTypeNum{"A3SS"} . ", A5SS:" . $asTypeNum{"A5SS"} . ", MXE:" . $asTypeNum{"MXE"} . ", RI:" . $asTypeNum{"RI"} . ", SE:" . $asTypeNum{"SE"} . "\n";

print "基于转录本结构推测的AS数量为：$trsptAsNum, ";
print "A3SS:" . $basedOnTrsptAsTypeNum{"A3SS"} . ", A5SS:" . $basedOnTrsptAsTypeNum{"A5SS"} . ", MXE:" . $basedOnTrsptAsTypeNum{"MXE"} . ", RI:" . $basedOnTrsptAsTypeNum{"RI"} . ", SE:" . $basedOnTrsptAsTypeNum{"SE"} . "\n";

print "直系同源AS数量为：$orthAsNum, ";
print "A3SS:" . $orthAsTypeNum{"A3SS"} . ", A5SS:" . $orthAsTypeNum{"A5SS"} . ", MXE:" . $orthAsTypeNum{"MXE"} . ", RI:" . $orthAsTypeNum{"RI"} . ", SE:" . $orthAsTypeNum{"SE"} . "\n";
