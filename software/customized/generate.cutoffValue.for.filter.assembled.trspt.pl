#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--assemblyDir psiOutput \\\n" .
                "--assembledExptTsv alignment.info.of.assembled.experiment.tsv \\\n" .
		"--unitVolume 6 \\\n" .
		"--unitTrsptMinCov 2 \\\n" .
		"--unitExonMinCov 1 \\\n" .
                "--outputExptCutoffTsv cutoff.info.of.assembled.experiment.tsv \n";
	exit;
}

my ($assemblyDir, $assembledExptTsv, $outputExptCutoffTsv, $unitVolume, $unitTrsptMinCov, $unitExonMinCov);

GetOptions(
        'assemblyDir=s'=>\$assemblyDir,
        'assembledExptTsv=s'=>\$assembledExptTsv,
        'outputExptCutoffTsv=s'=>\$outputExptCutoffTsv,
	'unitVolume=s'=>\$unitVolume,
	'unitTrsptMinCov=s'=>\$unitTrsptMinCov,
	'unitExonMinCov=s'=>\$unitExonMinCov,
);

my (@fieldName, @field, $i, $exptId, $exptIdPos);
my ($experimentLine, %experiment);

open FF, "<$assembledExptTsv";
open WW, ">$outputExptCutoffTsv";
$experimentLine=<FF>;
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

	print WW join("\t", $assemblyDir . "/" . $exptId . "/transcriptomeByStringtie.gtf", ${$experiment{$exptId}}{"mappedBases"}, sprintf("%.2f", ${$experiment{$exptId}}{"mappedBases"} / $unitVolume * $unitTrsptMinCov), sprintf("%.2f", ${$experiment{$exptId}}{"mappedBases"} / $unitVolume * $unitExonMinCov)) . "\n";
}
close FF;
close WW;

