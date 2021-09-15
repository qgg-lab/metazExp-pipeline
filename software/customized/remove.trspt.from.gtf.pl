#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV<0){
	print "\n\t perl $0 \\\n" .
		"\t\t --inputGtf 		final.combined.gtf \\\n" .
		"\t\t --removedTrsptIdList 	removed.trsptId.list \\\n" .
		"\t\t --outputRemainedGtf	noSorted.remained.gtf \\\n" .
		"\t\t 1>			log.o.remove.trspt.from.gtf \\\n" .
		"\t\t 2>			log.e.remove.trspt.from.gtf \n\n\n";
	exit();
}

my ($inputGtf, $removedTrsptIdList, $outputRemainedGtf);
GetOptions(
	'inputGtf=s'=>\$inputGtf,
	'removedTrsptIdList=s'=>\$removedTrsptIdList,
	'outputRemainedGtf=s'=>\$outputRemainedGtf,

);

# read gtf into trspt hash table
my (%referenceTrsptFullFeature);
my ($featureLine, @fields, $transcriptId, $geneId);
open FF, "<$inputGtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register feature into hash %referenceTrsptFullFeature
		$referenceTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}elsif($fields[2] eq "transcript"){
		#register feature into hash %referenceTrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$geneId = &getGeneId($fields[8]);
		$referenceTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}
}
close FF;

# read removed trspt list
my (%rmTrsptHash);
open FF, "<$removedTrsptIdList";
while(my $trsptId=<FF>){
	chomp($trsptId);
	$rmTrsptHash{$trsptId}=1;
}
close FF;

# output trspt from input gtf but without removed trspt
my @trsptId = keys(%referenceTrsptFullFeature);
open WW, ">$outputRemainedGtf";
foreach my $trsptId(@trsptId){
	next if(exists($rmTrsptHash{$trsptId}));
	print WW $referenceTrsptFullFeature{$trsptId};
}
close WW;



sub getTranscriptId{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_id "(.*)"/){
                        return $1;
                }
        }
        return "";
}
sub getGeneId{
        my ($singleTrsptAnnoText) = $_[0];
        my (@featureLines, $featureLine, @attrs);
        @featureLines = split(/\n/, $singleTrsptAnnoText);
        foreach $featureLine(@featureLines){
                @attrs = ();
                @attrs = split(/;/, $featureLine);
                foreach my $attr(@attrs){
                        if($attr=~/gene_id "(.*)"/){
                                return $1;
                        }
                }
        }
        return "";
}

