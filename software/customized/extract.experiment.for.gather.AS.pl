#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--assemblyDir psiOutput \\\n" .
                "--assembledExptTsv alignment.info.of.assembled.experiment.tsv \\\n" .
		"--minAligPercentage 50 \\\n" .
		"--minMappedBase 4 \\\n" .
		"--unitVolume 6 \\\n" .
		"--unitTrsptMinCov 2 \\\n" .
		"--unitExonMinCov 1 \\\n" .
		"--filteredAssembledExptTsv filtered.alignment.info.of.assembled.experiment.tsv \\\n" .
                "--outputExptCutoffTsv cutoff.info.of.assembled.experiment.tsv \n";
	exit;
}

my ($assemblyDir, $assembledExptTsv, $outputExptCutoffTsv, $unitVolume, $unitTrsptMinCov, $unitExonMinCov, $filteredAssembledExptTsv, $minAligPercentage, $minMappedBase);

GetOptions(
        'assemblyDir=s'=>\$assemblyDir,
        'assembledExptTsv=s'=>\$assembledExptTsv,
        'minAligPercentage=s'=>\$minAligPercentage,
        'minMappedBase=s'=>\$minMappedBase,
	'filteredAssembledExptTsv=s'=>\$filteredAssembledExptTsv,
        'outputExptCutoffTsv=s'=>\$outputExptCutoffTsv,
);

my (@fieldName, @field, $i, $exptId, $exptIdPos);
my ($experimentLine, %experiment);

open FF, "<$assembledExptTsv";
open FILTER, ">$filteredAssembledExptTsv";
open WW, ">$outputExptCutoffTsv";
$experimentLine=<FF>;
print FILTER $experimentLine;
chomp($experimentLine);
@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}

while($experimentLine=<FF>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);
        $exptId = $field[$exptIdPos];

        for($i=0; $i<=$#field; $i++){
                ${$experiment{$exptId}}{$fieldName[$i]} = $field[$i];
        }

	next if(${$experiment{$exptId}}{"mappedBases"}<$minMappedBase or ${$experiment{$exptId}}{"alignPercent"}<$minAligPercentage);

	print FILTER $experimentLine . "\n";
	print WW join("\t", $assemblyDir . "/" . $exptId . "/transcriptomeByStringtie.gtf", ${$experiment{$exptId}}{"mappedBases"}) . "\n";

}
close FF;
close WW;
close FILTER;
