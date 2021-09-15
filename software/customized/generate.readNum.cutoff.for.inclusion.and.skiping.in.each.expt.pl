#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--unitVolume 6 \\\n" .
		"--inclusionMinReadNumPerUnit 10 \\\n" .
		"--skipingMinReadNumPerUnit 10 \\\n" .
		"--exptInfoList expt.list.tsv \\\n" .
		"--exptReadNumCutoff exptReadNumCutoff.txt \n";
	exit;
}

my ($unitVolume, $inclusionMinReadNumPerUnit, $skipingMinReadNumPerUnit, $exptInfoList, $exptReadNumCutoff);

GetOptions(
        'unitVolume=s'=>\$unitVolume,
        'inclusionMinReadNumPerUnit=s'=>\$inclusionMinReadNumPerUnit,
        'skipingMinReadNumPerUnit=s'=>\$skipingMinReadNumPerUnit,
        'exptInfoList=s'=>\$exptInfoList,
        'exptReadNumCutoff=s'=>\$exptReadNumCutoff,
        );

my (@fieldName, @field, $i, $exptId, $exptIdPos);
my ($experimentLine, %experiment);

open FF, "<$exptInfoList";
open WW, ">$exptReadNumCutoff";
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

	print WW join("\t", $exptId, ${$experiment{$exptId}}{"mappedBases"}, sprintf("%.2f", ${$experiment{$exptId}}{"mappedBases"} / $unitVolume * $inclusionMinReadNumPerUnit), sprintf("%.2f", ${$experiment{$exptId}}{"mappedBases"} / $unitVolume * $skipingMinReadNumPerUnit)) . "\n";

}
close FF;
close WW;
