#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/sum/;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--filteredExptTsv \\\n" .
                "--outputStudyVolumeLevel \n";
	exit;
}

my ($filteredExptTsv, $outputStudyVolumeLevel,  $levelList);

GetOptions(
        'filteredExptTsv=s'=>\$filteredExptTsv,
	'levelList=s'=>\$levelList,
        'outputStudyVolumeLevel=s'=>\$outputStudyVolumeLevel,
);

my (%studyVolume, $studyValueHref, $studyId, @studyId);
$studyValueHref=\%studyVolume;
my (@valueField, @nameField, %tmpHash, $line, @volume);
open FF, "<$filteredExptTsv";
# Ecotype Cultivar Genotype Tissue SubTissue Development Treatment treatmentGroup Experiment Study dataSource Base Layout SpotsNum ReadNum SpotLen ReadLen Gather AS Assemble RunList Phenotype alignPercent mappedBases mappedReadNum detectedReadLen libraryType phredScore
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for(my $i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	$studyId = $tmpHash{"Study"};
	$studyValueHref->{$studyId}->{"mappedVolumeList"} .= $tmpHash{"mappedBases"} . ",";
}
close FF;

my ($level, %level, $covertLevel);
my @level=split(/,/, $levelList);
@studyId = keys(%studyVolume);
foreach $studyId(@studyId){
	@volume = ();
	@volume = split(/,/, $studyValueHref->{$studyId}->{"mappedVolumeList"});
	$studyValueHref->{$studyId}->{"meanMappedVolum"} = sum(@volume)/($#volume+1);
	$studyValueHref->{$studyId}->{"exptNum"} = $#volume + 1;
	foreach $level(@level){
		if($studyValueHref->{$studyId}->{"meanMappedVolum"} <= $level){
			$covertLevel = sprintf("%.2f", $level);
			if(index($covertLevel, ".") == 1){
				$covertLevel = "0" . $covertLevel;
			}
			$studyValueHref->{$studyId}->{"level"} = $covertLevel . "GB";
			$level{$level}++;
			last;
		}
	}
}

=beg
foreach $level(@level){
	print join("\t", $level, $level{$level}) . "\n" if(exists($level{$level}));
}
=cut;
open WW, ">$outputStudyVolumeLevel";
print WW "studyId\texptNum\tmeanMappedvolum\tlevel\tscope\n";
foreach $studyId(@studyId){
        print WW join("\t", $studyId, $studyValueHref->{$studyId}->{"exptNum"}, sprintf("%.4f", $studyValueHref->{$studyId}->{"meanMappedVolum"}), $studyValueHref->{$studyId}->{"level"}, "NA") . "\n";
}
close WW;

