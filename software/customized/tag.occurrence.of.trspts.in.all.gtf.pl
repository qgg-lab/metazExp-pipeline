#!/usr/bin/perl
use strict;
use Getopt::Long;
 
if($#ARGV<0){
	print "\n\tperl $0 \\\n" . 
		"\t\t --InputCmprGtfFile final.complete.trspt.anno.gtf \\\n" . 
		"\t\t --InputOrigGtfFileList ensembl.gtf,refseq.gtf\\\n" . 
		"\t\t --InputOrigGtfLabelList ensembl,refseq \\\n" . 
		"\t\t --outputTagFile	final.complete.trspt.with.origTag.tsv\n\n";
	print "\t\t This script is used to tag occurence of each trsipt in original gtf\n\n";
	exit;
}

my ($InputCmprGtfFile, $InputOrigGtfFileList, $InputOrigGtfLabelList, $outputTagFile);
GetOptions(
	'InputCmprGtfFile=s'=>\$InputCmprGtfFile,
	'InputOrigGtfFileList=s'=>\$InputOrigGtfFileList,
	'InputOrigGtfLabelList=s'=>\$InputOrigGtfLabelList,
	'outputTagFile=s'=>\$outputTagFile,
);

# read concat exon coord of comprGtfFile into hash. The key of hash is the concat exon coord.
my (%cmprConcatExonCoorToTrsptId, $concatExonCoorText, $concatExonCoorLine, @concatExonCoorLine);
my (@tt);

$concatExonCoorText=&getConcatExonCoorMultiText($InputCmprGtfFile);

@concatExonCoorLine = ();
@concatExonCoorLine = split(/\n/, $concatExonCoorText);
foreach $concatExonCoorLine(@concatExonCoorLine){
	@tt = ();
	@tt = split(/\t/, $concatExonCoorLine);
	${$cmprConcatExonCoorToTrsptId{$tt[1]}}{"trsptId"} = $tt[0];
	${$cmprConcatExonCoorToTrsptId{$tt[1]}}{"scl"} = $tt[2];
}


# read concate exon coord of each inputOrigGtfFile
my (@inputOrigGtfFile, @inputOrigGtfLabel, %concatExonCoor, @concatExonCoor, $concatExonCoor);

@inputOrigGtfFile = split(/,/, $InputOrigGtfFileList);
@inputOrigGtfLabel = split(/,/, $InputOrigGtfLabelList);

for(my $i=0; $i<=$#inputOrigGtfFile; $i++){
	$concatExonCoorText=&getConcatExonCoorMultiText($inputOrigGtfFile[$i]);
	# read concat exon coor into %concatExonCoor
	@concatExonCoorLine = ();
	@concatExonCoorLine = split(/\n/, $concatExonCoorText);
	%concatExonCoor = ();
	# each concatExonCoorLine correspond to trspt
	foreach $concatExonCoorLine(@concatExonCoorLine){
		@tt = ();
		@tt = split(/\t/, $concatExonCoorLine);	
		$concatExonCoor{$tt[1]}=$tt[0];
	}

	# read each trspt in cmprGtfFile and judge whether its concat exon coor exists in concatExonCoor hash
	@concatExonCoor = ();
	@concatExonCoor = keys(%cmprConcatExonCoorToTrsptId);
	foreach $concatExonCoor(@concatExonCoor){
		if(exists($concatExonCoor{$concatExonCoor})){
			${$cmprConcatExonCoorToTrsptId{$concatExonCoor}}{"orign"}.= $concatExonCoor{$concatExonCoor} . "\t";
		}else{
			${$cmprConcatExonCoorToTrsptId{$concatExonCoor}}{"orign"}.= "-\t";
		}
	}
}

# output each trspt id in cmprGtf with origi label
@concatExonCoor = ();
@concatExonCoor = keys(%cmprConcatExonCoorToTrsptId);
open WW, ">$outputTagFile";
print WW "TrsptId\t";
for(my $i=0; $i<=$#inputOrigGtfLabel; $i++){
	print WW $inputOrigGtfLabel[$i] . "\t";
}
print WW  "\n";

foreach $concatExonCoor(@concatExonCoor){
	print WW join("\t", ${$cmprConcatExonCoorToTrsptId{$concatExonCoor}}{"trsptId"}, 
			${$cmprConcatExonCoorToTrsptId{$concatExonCoor}}{"scl"}, 
			${$cmprConcatExonCoorToTrsptId{$concatExonCoor}}{"orign"}) . "\n";
}
close WW;



sub getConcatExonCoorMultiText{
	my ($gtfFile)=$_[0];
	my (%trspt, $rltText);
	my ($featureLine, @fields, $transcriptId, @transcriptId);

	$rltText = "";
	open FF, "<$gtfFile";
	while($featureLine =<FF>){
		chomp($featureLine);
		@fields = ();
		@fields = split(/\t/, $featureLine);
		$transcriptId = &getTranscriptId($fields[8]);

		if($fields[2] eq "exon"){
			#concat transcript exon coord
			$trspt{$transcriptId} .= join("___", $fields[0], $fields[3], $fields[4], $fields[6]) . "!";
		}
	}
	close FF;

	#sort exon according to coordinate and load i
	my (@exonCoordinateFullLine, $exonCoordinateFullLine);
	my (@srtExon);
	@transcriptId = ();
	@transcriptId = keys(%trspt);
	foreach $transcriptId(@transcriptId){

		@exonCoordinateFullLine = ();
		@srtExon = ();

		# split exon coordinate full into 1-dim exon array
		@exonCoordinateFullLine = split(/!/, $trspt{$transcriptId});

		# split exon coordinate into 2-dim exon array
		foreach my $exonCoordinate(@exonCoordinateFullLine){
			push @srtExon, [split '___', $exonCoordinate];	
		}

		# sort by coordinate
		if($srtExon[0][3] eq "+"){
			@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
		}else{
			@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
		}

		# obtain concat exon coor text of a trspt
		my $concatExonCoorText = "";
		for(my $i=0; $i<=$#srtExon; $i++){
			$concatExonCoorText.=join("___", $srtExon[$i][0], $srtExon[$i][1],
						 $srtExon[$i][2], $srtExon[$i][3]) . "!";
		}

		# append concat exon coordinate line of a trspt into rltText
		my $trsptScl = "";
		if($srtExon[0][3] eq "+"){
			$trsptScl = $srtExon[0][0] . ":" . $srtExon[0][1] . "-" . $srtExon[$#srtExon][2] . " (" .  $srtExon[0][3] . ")";
		}else{
			$trsptScl = $srtExon[0][0] . ":" . $srtExon[$#srtExon][2] . "-" . $srtExon[0][1] . " (" .  $srtExon[0][3] . ")";
		}

		$rltText.= $transcriptId . "\t" . $concatExonCoorText . "\t" . $trsptScl . "\n";
	}
	return $rltText;
}


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

