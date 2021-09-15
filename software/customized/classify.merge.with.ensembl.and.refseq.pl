#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV < 0){
	print "perl $0 \\\n" . 
		"\t\t--mergedAnno  /mnt/home/liujind1/workAS/01-cattle/database/merged.annotation.gtf \\\n" .
		"\t\t--anno1   /mnt/home/liujind1/workAS/01-cattle/database/ensembl.annotation.gtf \\\n" . 
		"\t\t--anno2   /mnt/home/liujind1/workAS/01-cattle/database/refSeq.annotation.with_ensemblSeqId.gtf \\\n" . 
		"\t\t--outputAnno1 /mnt/home/liujind1/workAS/01-cattle/database/ensembl.in.merged.with.oriId.gtf \\\n" .
		"\t\t--outputAnno2 /mnt/home/liujind1/workAS/01-cattle/database/refSeq.in.merged.with.oriId.gtf \\\n" .
		"\t\t--outputOverlap /mnt/home/liujind1/workAS/01-cattle/database/ensembl.overlap.refSeq.in.merged.with.combined.OriId.gtf \\\n" .
		"\t\t--outputNewAnno  /mnt/home/liujind1/workAS/01-cattle/database/new.gtf \\\n" .
		"\t\t--outputInAnno2NotInAnno1  /mnt/home/liujind1/workAS/01-cattle/database/in.anno2.not.in.anno2.gtf \\\n" .
		"\t\t--outputMergedAnno /mnt/home/liujind1/workAS/01-cattle/database/mergedAnno.with.origAnno1Id.andThen.origAnno2Id.gtf \n";
	exit(0);  
}

my ($mergedAnno, $anno1, $anno2, $outputAnno1, 
	$outputAnno2, $outputOverlap, 
	$outputNewAnno,$outputMergedAnno, $outputInAnno2NotInAnno1);

GetOptions(
	'mergedAnno=s'=>\$mergedAnno, 
	'anno1=s'=>\$anno1, 
	'anno2=s'=>\$anno2, 
	'outputAnno1=s'=>\$outputAnno1, 
	'outputAnno2=s'=>\$outputAnno2, 
	'outputOverlap=s'=>\$outputOverlap, 
	'outputNewAnno=s'=>\$outputNewAnno,
	'outputInAnno2NotInAnno1=s'=>\$outputInAnno2NotInAnno1,
	'outputMergedAnno=s'=>\$outputMergedAnno,
);

my (%anno1TrsptIdToExonCoord, %anno2TrsptIdToExonCoord, %mergeTrsptIdToExonCoord);
my (%anno1ExonCoordToTrsptId, %anno2ExonCoordToTrsptId, %mergeExonCoordToTrsptId);
my (%anno1TrsptFullFeature, %anno2TrsptFullFeature, %mergeTrsptFullFeature);
my (@fields, @attrs, $featureLine);
my ($transcriptId);
my (@transcriptId, @exonCoordinateFullLine, @srtExon);

#read anno1 exon coordinates into hash
#read anno1 tanscript and exon into hash %anno1Transcript
open FF, "<$anno1";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register coordinate into hash %anno1TrsptIdToExonCoord
		$anno1TrsptIdToExonCoord{$transcriptId}.= join("\t", $fields[0], $fields[3], $fields[4], $fields[6]) . "\n";
		#register feature into hash %anno1TrsptFullFeature
		$anno1TrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}elsif($fields[2] eq "transcript"){
		#register feature into hash %anno1TrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$anno1TrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}
}
close FF;
#sort exon according to coordinate and load into hash ExonCoordToTrsptId
@transcriptId = ();
@transcriptId = keys(%anno1TrsptIdToExonCoord);
foreach $transcriptId(@transcriptId){
	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/\n/, $anno1TrsptIdToExonCoord{$transcriptId});
	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '\t', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#re-generate exon coordinate full line string
	my $coordinateFullLineString = "";
	for(my $i=0; $i<=$#srtExon; $i++){
		$coordinateFullLineString.=join("#", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][2], $srtExon[$i][3]) . "#";
	}
	$coordinateFullLineString = substr($coordinateFullLineString, 0, length($coordinateFullLineString) -1);

	$anno1ExonCoordToTrsptId{$coordinateFullLineString} = $transcriptId;
}

#read anno2 exon coordinates into hash
#read anno2 tanscript and exon into hash %anno2Transcript
open FF, "<$anno2";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register coordinate into hash %anno2TrsptIdToExonCoord
		$anno2TrsptIdToExonCoord{$transcriptId}.= join("\t", $fields[0], $fields[3], $fields[4], $fields[6]) . "\n";
		#register feature into hash %anno2TrsptFullFeature
		$anno2TrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}elsif($fields[2] eq "transcript"){
		#register feature into hash %anno2TrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$anno2TrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}
}
close FF;
#sort exon according to coordinate and load into hash ExonCoordToTrsptId
@transcriptId = ();
@transcriptId = keys(%anno2TrsptIdToExonCoord);
foreach $transcriptId(@transcriptId){
	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/\n/, $anno2TrsptIdToExonCoord{$transcriptId});
	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '\t', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#re-generate exon coordinate full line string
	my $coordinateFullLineString = "";
	for(my $i=0; $i<=$#srtExon; $i++){
		$coordinateFullLineString.=join("#", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][2], $srtExon[$i][3]) . "#";
	}
	$coordinateFullLineString = substr($coordinateFullLineString, 0, length($coordinateFullLineString) -1);

	$anno2ExonCoordToTrsptId{$coordinateFullLineString} = $transcriptId;
}

#read merge exon coordinates into hash
#read merge tanscript and exon into hash %mergeTranscript
open FF, "<$mergedAnno";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register coordinate into hash %merge2TrsptIdToExonCoord
		$mergeTrsptIdToExonCoord{$transcriptId}.= join("\t", $fields[0], $fields[3], $fields[4], $fields[6]) . "\n";
		#register feature into hash %mergeTrsptFullFeature
		$mergeTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}elsif($fields[2] eq "transcript"){
		#register feature into hash %mergeTrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$mergeTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
	}
}
close FF;
#sort exon according to coordinate and load into hash ExonCoordToTrsptId
@transcriptId = ();
@transcriptId = keys(%mergeTrsptIdToExonCoord);
foreach $transcriptId(@transcriptId){
	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/\n/, $mergeTrsptIdToExonCoord{$transcriptId});
	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '\t', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#re-generate exon coordinate full line string
	my $coordinateFullLineString = "";
	for(my $i=0; $i<=$#srtExon; $i++){
		$coordinateFullLineString.=join("#", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][2], $srtExon[$i][3]) . "#";
	}
	$coordinateFullLineString = substr($coordinateFullLineString, 0, length($coordinateFullLineString) -1);

	$mergeExonCoordToTrsptId{$coordinateFullLineString} = $transcriptId;
}

#output anno1 in merge
#check transcript both in anno1 that also exist in merge anno
#Output transcript full feature in anno1 
#the cooresponding original gene id should be change merge id.
my @coordinateFullLineString = keys(%anno1ExonCoordToTrsptId);
my ($outputAnno1AllTrsptText, $anno1SingleTrsptAnnoText);
my ($anno1TrsptId, $mergeTrsptId, $mergeGeneId, $anno1GeneId);
my ($anno1SingleTrsptAnnoText);
$outputAnno1AllTrsptText = "";
foreach my $coordinateFullLineString(@coordinateFullLineString){
	#check coordinate of a transcript in anno1 both in merge anno
	if(exists($mergeExonCoordToTrsptId{$coordinateFullLineString})){
		
		$anno1TrsptId = $anno1ExonCoordToTrsptId{$coordinateFullLineString};
		$mergeTrsptId = $mergeExonCoordToTrsptId{$coordinateFullLineString};
		$anno1SingleTrsptAnnoText = $anno1TrsptFullFeature{$anno1TrsptId};
		
		$anno1GeneId = &getGeneId($anno1TrsptFullFeature{$anno1TrsptId});
		$mergeGeneId = &getGeneId($mergeTrsptFullFeature{$mergeTrsptId});
		
		$anno1GeneId=~s/\./\\\./g;
		
		$anno1SingleTrsptAnnoText=~s/gene_id "$anno1GeneId";/gene_id "$mergeGeneId";/g;

		$outputAnno1AllTrsptText.=$anno1SingleTrsptAnnoText;
	}
}
open WW, ">" . $outputAnno1;
print WW $outputAnno1AllTrsptText;
close WW;

#output anno2 in merge
#check transcript both in anno2 that also exist in merge anno
#Output transcript full feature in anno2 
#the cooresponding original gene id should be change merge id.
my @coordinateFullLineString = keys(%anno2ExonCoordToTrsptId);
my ($outputAnno2AllTrsptText, $anno2SingleTrsptAnnoText);
my ($anno2TrsptId, $mergeTrsptId, $mergeGeneId, $anno2GeneId);
my ($anno2SingleTrsptAnnoText);
$outputAnno2AllTrsptText = "";
foreach my $coordinateFullLineString(@coordinateFullLineString){
	#check coordinate of a transcript in anno2 both in merge anno
	if(exists($mergeExonCoordToTrsptId{$coordinateFullLineString})){
		
		$anno2TrsptId = $anno2ExonCoordToTrsptId{$coordinateFullLineString};
		$mergeTrsptId = $mergeExonCoordToTrsptId{$coordinateFullLineString};
		$anno2SingleTrsptAnnoText = $anno2TrsptFullFeature{$anno2TrsptId};
		
		$anno2GeneId = &getGeneId($anno2TrsptFullFeature{$anno2TrsptId});
		$mergeGeneId = &getGeneId($mergeTrsptFullFeature{$mergeTrsptId});
		
		$anno2GeneId=~s/\./\\\./g;
		
		$anno2SingleTrsptAnnoText=~s/gene_id "$anno2GeneId";/gene_id "$mergeGeneId";/g;

		$outputAnno2AllTrsptText.=$anno2SingleTrsptAnnoText;
	}
}
open WW, ">" . $outputAnno2;
print WW $outputAnno2AllTrsptText;
close WW;

#output overlap of anno1 and anno2 in merge
#check transcript both in anno1 and in anno2
#Output transcript full feature in merge 
#the cooresponding original gene id should be change merge id.
my @coordinateFullLineString = keys(%mergeExonCoordToTrsptId);
my ($outputOverlapAllTrsptText, $overlapSingleTrsptAnnoText);
my ($anno1TrsptId, $anno2TrsptId, $mergeTrsptId, $mergeGeneId, $anno1GeneId, $anno2GeneId, $combinedTrsptId);

$outputOverlapAllTrsptText = "";
foreach my $coordinateFullLineString(@coordinateFullLineString){
	#check coordinate of a transcript in anno1 both in merge anno
	if(exists($mergeExonCoordToTrsptId{$coordinateFullLineString}) and exists($anno1ExonCoordToTrsptId{$coordinateFullLineString}) and exists($anno2ExonCoordToTrsptId{$coordinateFullLineString})){
		
		$anno1TrsptId = $anno1ExonCoordToTrsptId{$coordinateFullLineString};
		$anno2TrsptId = $anno2ExonCoordToTrsptId{$coordinateFullLineString};
		$mergeTrsptId = $mergeExonCoordToTrsptId{$coordinateFullLineString};
		$combinedTrsptId = $anno1TrsptId . "__" . $anno2TrsptId;
		
		$mergeGeneId = &getGeneId($mergeTrsptFullFeature{$mergeTrsptId});
		$anno1GeneId = &getGeneId($anno1TrsptFullFeature{$anno1TrsptId});
	
		$overlapSingleTrsptAnnoText = $anno1TrsptFullFeature{$anno1TrsptId};		
		
		$anno1GeneId =~s/\./\\\./g;
		$overlapSingleTrsptAnnoText=~s/gene_id "$anno1GeneId";/gene_id "$mergeGeneId";/g;

		$anno1TrsptId=~s/\./\\\./g;
		$overlapSingleTrsptAnnoText=~s/transcript_id "$anno1TrsptId";/transcript_id "$combinedTrsptId";/g;

		$outputOverlapAllTrsptText.=$overlapSingleTrsptAnnoText;
	}
}
open WW, ">" . $outputOverlap;
print WW $outputOverlapAllTrsptText;
close WW;

#output transcript in anno2 but not anno1
#Output transcript full feature in anno2 
#the cooresponding original gene id should be change merge id.
my @coordinateFullLineString = keys(%mergeExonCoordToTrsptId);
my ($outputInanno2NotInAnno1AllTrsptText, $inAnno2NotInAnno1SingleTrsptAnnoText);
my ($anno1TrsptId, $anno2TrsptId, $mergeTrsptId, $mergeGeneId, $anno1GeneId, $anno2GeneId, $combinedTrsptId);

$outputInanno2NotInAnno1AllTrsptText = "";
foreach my $coordinateFullLineString(@coordinateFullLineString){
	#check coordinate of a transcript in anno1 both in merge anno
	if(exists($mergeExonCoordToTrsptId{$coordinateFullLineString}) and not exists($anno1ExonCoordToTrsptId{$coordinateFullLineString}) and exists($anno2ExonCoordToTrsptId{$coordinateFullLineString})){
		
		$anno2TrsptId = $anno2ExonCoordToTrsptId{$coordinateFullLineString};
		$mergeTrsptId = $mergeExonCoordToTrsptId{$coordinateFullLineString};
		
		$mergeGeneId = &getGeneId($mergeTrsptFullFeature{$mergeTrsptId});
		$anno2GeneId = &getGeneId($anno2TrsptFullFeature{$anno2TrsptId});
	
		$inAnno2NotInAnno1SingleTrsptAnnoText = $anno2TrsptFullFeature{$anno2TrsptId};		
		
		$anno2GeneId =~s/\./\\\./g;
		$inAnno2NotInAnno1SingleTrsptAnnoText =~s/gene_id "$anno2GeneId";/gene_id "$mergeGeneId";/g;

		$outputInanno2NotInAnno1AllTrsptText .=$inAnno2NotInAnno1SingleTrsptAnnoText;
	}
}
open WW, ">" . $outputInAnno2NotInAnno1;
print WW $outputInanno2NotInAnno1AllTrsptText;
close WW;

#output transcript both not in anno1 and anno2
#Output transcript full feature in merge 
my @coordinateFullLineString = keys(%mergeExonCoordToTrsptId);
my ($outputNewAllTrsptText, $newSingleTrsptAnnoText);
my ($anno1TrsptId, $anno2TrsptId, $mergeTrsptId, $mergeGeneId, $anno1GeneId, $anno2GeneId);

$outputNewAllTrsptText = "";
foreach my $coordinateFullLineString(@coordinateFullLineString){
	#check coordinate of a transcript in anno1 both in merge anno
	if(exists($mergeExonCoordToTrsptId{$coordinateFullLineString}) and not exists($anno1ExonCoordToTrsptId{$coordinateFullLineString}) and not exists($anno2ExonCoordToTrsptId{$coordinateFullLineString})){
		
		$mergeTrsptId = $mergeExonCoordToTrsptId{$coordinateFullLineString};
		$newSingleTrsptAnnoText = $mergeTrsptFullFeature{$mergeTrsptId};		
		$outputNewAllTrsptText .=$newSingleTrsptAnnoText;
	}
}
open WW, ">" . $outputNewAnno;
print WW $outputNewAllTrsptText;
close WW;



#output merge with original transcript Id. Prior to use Ensembl Id
open WW, ">$outputMergedAnno";
print WW $outputAnno1AllTrsptText;
print WW $outputNewAllTrsptText;
print WW  $outputInanno2NotInAnno1AllTrsptText;
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
