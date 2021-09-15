#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV < 0){
	print "perl $0 \\\n" . 
		"\t\t--mergedAnnoGtf combinedAnnotation.gtf \\\n" .
		"\t\t--mergeLabel     Merge \\\n" .
		"\t\t--anno1Gtf   sorted.nonRedundant.ensembl.gtf \\\n" . 
		"\t\t--anno1Label     Ensembl \\\n" .
		"\t\t--anno2Gtf   sorted.nonRedundant.refseq.gtf \\\n" . 
		"\t\t--anno2Label     Refseq \\\n" . 
		"\t\t--anno3Gtf sorted.nonRedundant.jgi.gtf \\\n" .
		"\t\t--anno3Label jgi \\\n" .
		"\t\t--xvfb-run  /usr/bin/xvfb-run  \\\n" .
		"\t\t--Rscript   /opt/software/R/3.5.1-foss-2018a-X11-20180131/bin/Rscript \\\n" . 
		"\t\t--drawVennRscript draw.vennDiagram.for.four.sets.R \\\n" . 
		"\t\t--drawBarRscript draw.barDiagram.for.four.sets.R \\\n" .
		"\t\t--drawDensityRscript draw.densityDiagram.for.four.sets.R  \\\n" .
		"\t\t--densityExonMinLen 1 \\\n" .
		"\t\t--densityExonMaxLen 200 \\\n" .
		"\t\t--densityTranscriptMinLen 100 \\\n" .
		"\t\t--densityTranscriptMaxLen 5000 \\\n" .

		"\t\t--outputStaTbl combined.annotation.gtf.full.sta \\\n" .

		"\t\t--outputTrsptEnsRefMergVenn venn.of.trspt.png \\\n" .
		"\t\t--outputExonEnsRefMergVenn  venn.of.exon.png \\\n" .
		"\t\t--outputDonorEnsRefMergVenn venn.of.donor.png \\\n" .
		"\t\t--outputAcceptorEnsRefMergVenn venn.of.acceptor.png \\\n" .
		"\t\t--outputTrsptLenEnsRefMergDensity density.of.trsptLen.png \\\n" .
		"\t\t--outputExonLenEnsRefMergDensity density.of.exonLen.png \\\n" .
		"\t\t--outputTrsptNumEnsRefMergBar bar.of.gene.with.trsptNum.png \\\n" .
		"\t\t--outputExonNumEnsRefMergBar bar.of.trspt.with.exonNum.png \\\n" .
		"\t\t--outputImagType png \\\n" . 
		"\t\t--tmpOutputDir tmpStatistic \n";
	exit(0);  
}

my($mergedAnnoGtf, $anno1Gtf, $anno2Gtf, $anno3Gtf,
	$mergeLabel, $anno1Label, $anno2Label, $anno3Label,
	$xvfb_run,
	$Rscript, $drawVennRscript, $drawDensityRscript, $drawBarRscript,
	$densityExonMinLen, $densityExonMaxLen,
	$densityTranscriptMinLen, $densityTranscriptMaxLen,
	$outputStaTbl,
	$outputTrsptEnsRefMergVenn, $outputExonEnsRefMergVenn,
	$outputDonorEnsRefMergVenn, $outputAcceptorEnsRefMergVenn,
	$outputTrsptLenEnsRefMergDensity, $outputTrsptNumEnsRefMergBar,
	$outputExonLenEnsRefMergDensity, $outputExonNumEnsRefMergBar,
	$outputImagType, $tmpOutputDir);

GetOptions(
	'mergedAnnoGtf=s'=>\$mergedAnnoGtf,
	'mergeLabel=s'=>\$mergeLabel,
	'anno1Gtf=s'=>\$anno1Gtf, 
	'anno1Label=s'=>\$anno1Label,
	'anno2Gtf=s'=>\$anno2Gtf,
	'anno2Label=s'=>\$anno2Label,
	'anno3Gtf=s'=>\$anno3Gtf,
	'anno3Label=s'=>\$anno3Label,

	'xvfb-run=s'=>\$xvfb_run,
	'Rscript=s'=>\$Rscript,
	'drawVennRscript=s'=>\$drawVennRscript,
	'drawDensityRscript=s'=>\$drawDensityRscript,
	'densityExonMinLen=i'=>\$densityExonMinLen,
	'densityExonMaxLen=i'=>\$densityExonMaxLen,

	'densityTranscriptMinLen=i'=>\$densityTranscriptMinLen,
	'densityTranscriptMaxLen=i'=>\$densityTranscriptMaxLen,

	'drawBarRscript=s'=>\$drawBarRscript,

	'outputStaTbl=s'=>\$outputStaTbl,

	'outputTrsptEnsRefMergVenn=s'=>\$outputTrsptEnsRefMergVenn,
	'outputExonEnsRefMergVenn=s'=>\$outputExonEnsRefMergVenn,

	'outputDonorEnsRefMergVenn=s'=>\$outputDonorEnsRefMergVenn,
	'outputAcceptorEnsRefMergVenn=s'=>\$outputAcceptorEnsRefMergVenn,

	'outputTrsptLenEnsRefMergDensity=s'=>\$outputTrsptLenEnsRefMergDensity,
	'outputExonLenEnsRefMergDensity=s'=>\$outputExonLenEnsRefMergDensity,

	'outputTrsptNumEnsRefMergBar=s'=>\$outputTrsptNumEnsRefMergBar,
	'outputExonNumEnsRefMergBar=s'=>\$outputExonNumEnsRefMergBar,

	'outputImagType=s'=>\$outputImagType,
	'tmpOutputDir=s'=>\$tmpOutputDir
);

my (%anno1GeneIdToTrsptNum, %anno2GeneIdToTrsptNum, %anno3GeneIdToTrsptNum, %mergeGeneIdToTrsptNum);

my (%anno1TrsptIdToConcatExonCoord, %anno2TrsptIdToConcatExonCoord, %anno3TrsptIdToConcatExonCoord, %mergeTrsptIdToConcatExonCoord);
my (%anno1TrsptIdToTrsptLen, %anno2TrsptIdToTrsptLen, %anno3TrsptIdToTrsptLen, %mergeTrsptIdToTrsptLen);

my (@anno1ExonCoord, @anno2ExonCoord, @anno3ExonCoord, @mergeExonCoord);
my (@anno1ExonLen, @anno2ExonLen, @anno3ExonLen, @mergeExonLen);
my (@anno1DonorCoord, @anno2DonorCoord, @anno3DonorCoord, @mergeDonorCoord);
my (@anno1AcceptorCoord, @anno2AcceptorCoord, @anno3AcceptorCoord, @mergeAcceptorCoord);
my (%anno1TranscriptNum, %anno2TranscriptNum, %anno3TranscriptNum, %mergeTranscriptNum);
my (%anno1TrsptIdToExonNum, %anno2TrsptIdToExonNum, %anno3TrsptIdToExonNum, %mergeTrsptIdToExonNum);

my (@fields, @attrs, $featureLine);
my ($geneId, $transcriptId);
my (@geneId, @transcriptId, @exonCoordinateFullLine, @srtExon, $exonCoord, $exonLen);
my ($donorCoord, $acceptorCoord);

system("mkdir -p $tmpOutputDir");

open TOTALTRANSCRIPTLEN, ">$tmpOutputDir/total.transcript.length.txt";

################################
#
# ======== anno1 ============ 
#
###############################
##read anno1Gtf into hashes, arrays
open FF, "<$anno1Gtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	$geneId = &getGeneId($fields[8]);
	$transcriptId = &getTranscriptId($fields[8]);
	if($fields[2] eq "transcript"){

		#The num of gene with this transcript was increased 1
		$anno1GeneIdToTrsptNum{$geneId}++;

	}elsif($fields[2] eq "exon"){

		#concat transcript exon coord
		$anno1TrsptIdToConcatExonCoord{$transcriptId} .= join("___", $fields[0], $fields[3], $fields[4], $fields[6]) . "!";
		#increase transcript length with this exon length
		$anno1TrsptIdToTrsptLen{$transcriptId} += $fields[4] - $fields[3] + 1;
		#regist exon coord into array
		$anno1ExonCoord[$#anno1ExonCoord + 1] = join("___", $fields[0], $fields[3], $fields[4], $fields[6]);
		#regist exon length into array
		$anno1ExonLen[$#anno1ExonLen + 1] = $fields[4] - $fields[3] + 1;
		#register exon num of transcript
		$anno1TrsptIdToExonNum{$transcriptId}++;
	}
}
close FF;

#sort exon according to coordinate and load i
my (@exonCoordinateFullLine, $exonCoordinateFullLine);
@transcriptId = ();
@transcriptId = keys(%anno1TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){

	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/!/, $anno1TrsptIdToConcatExonCoord{$transcriptId});

	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '___', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#register donor
	for(my $i=0; $i<=$#srtExon - 1; $i++){
		$anno1DonorCoord[$#anno1DonorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][2], $srtExon[$i][3]);
	}
	#register acceptor
	for(my $i=1; $i<=$#srtExon; $i++){
		$anno1AcceptorCoord[$#anno1AcceptorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][3]);
	}
}

#output transcriptNum, transcriptConcatExonCoord, transcriptLength, exonCoord, exonLen, donorCoord, acceptorCoord
open TRANSCRIPTCONCATEXONCOORD, ">$tmpOutputDir/$anno1Label.transcript.concat.exon.coord.txt";
open EXONCOORD, ">$tmpOutputDir/$anno1Label.exon.coord.txt";
open DONORCOORD, ">$tmpOutputDir/$anno1Label.donor.coord.txt";
open ACCEPTORCOORD, ">$tmpOutputDir/$anno1Label.acceptor.coord.txt";

#transcript num of gene
@geneId = keys(%anno1GeneIdToTrsptNum);
foreach $geneId(@geneId){
	if($anno1GeneIdToTrsptNum{$geneId}<=5){
		$anno1TranscriptNum{$anno1GeneIdToTrsptNum{$geneId}} += 1;
	}else{
		$anno1TranscriptNum{">5"} += 1;
	}
}

#transcript concat exon coord
@transcriptId =();
@transcriptId = keys(%anno1TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){
	print TRANSCRIPTCONCATEXONCOORD $anno1TrsptIdToConcatExonCoord{$transcriptId} . "\n";
}

#transcript len
@transcriptId =();
@transcriptId = keys(%anno1TrsptIdToTrsptLen);
foreach $transcriptId(@transcriptId){
	print TOTALTRANSCRIPTLEN join("\t", $anno1TrsptIdToTrsptLen{$transcriptId}, $anno1Label) . "\n";
}

#exon coordinate
foreach $exonCoord(@anno1ExonCoord){
	print EXONCOORD $exonCoord . "\n";
}

#exon length
#foreach $exonLen(@anno1ExonLen){
#	print TOTALEXONLEN join("\t", $exonLen, $anno1Label ). "\n";
#}

#Donor site
foreach $donorCoord(@anno1DonorCoord){
	print DONORCOORD $donorCoord . "\n";
}

#Acceptor site
foreach $acceptorCoord(@anno1AcceptorCoord){
	print ACCEPTORCOORD $acceptorCoord . "\n";
}

#close TRANSCRIPTNU;
close TRANSCRIPTCONCATEXONCOORD;
#close TRANSCRIPTLEN;
close EXONCOORD;
#close EXONLEN;
close DONORCOORD;
close ACCEPTORCOORD;

################################
#
# ======== anno2 ============ 
#
###############################
##read anno2Gtf into hashes, arrays
open FF, "<$anno2Gtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	$geneId = &getGeneId($fields[8]);
	$transcriptId = &getTranscriptId($fields[8]);
	if($fields[2] eq "transcript"){

		#The num of gene with this transcript was increased 1
		$anno2GeneIdToTrsptNum{$geneId}++;

	}elsif($fields[2] eq "exon"){

		#concat transcript exon coord
		$anno2TrsptIdToConcatExonCoord{$transcriptId} .= join("___", $fields[0], $fields[3], $fields[4], $fields[6]) . "!";
		#increase transcript length with this exon length
		$anno2TrsptIdToTrsptLen{$transcriptId} += $fields[4] - $fields[3] + 1;
		#regist exon coord into array
		$anno2ExonCoord[$#anno2ExonCoord + 1] = join("___", $fields[0], $fields[3], $fields[4], $fields[6]);
		#regist exon length into array
		$anno2ExonLen[$#anno2ExonLen + 1] = $fields[4] - $fields[3] + 1;
		#register exon num of transcript
		$anno2TrsptIdToExonNum{$transcriptId}++;

	}
}
close FF;

#sort exon according to coordinate and load i
my (@exonCoordinateFullLine, $exonCoordinateFullLine);
@transcriptId = ();
@transcriptId = keys(%anno2TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){

	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/!/, $anno2TrsptIdToConcatExonCoord{$transcriptId});

	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '___', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#register donor
	for(my $i=0; $i<=$#srtExon - 1; $i++){
		$anno2DonorCoord[$#anno2DonorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][2], $srtExon[$i][3]);
	}
	#register acceptor
	for(my $i=1; $i<=$#srtExon; $i++){
		$anno2AcceptorCoord[$#anno2AcceptorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][3]);
	}
}

#output transcriptNum, transcriptConcatExonCoord, transcriptLength, exonCoord, exonLen, donorCoord, acceptorCoord
open TRANSCRIPTCONCATEXONCOORD, ">$tmpOutputDir/$anno2Label.transcript.concat.exon.coord.txt";
open EXONCOORD, ">$tmpOutputDir/$anno2Label.exon.coord.txt";
open DONORCOORD, ">$tmpOutputDir/$anno2Label.donor.coord.txt";
open ACCEPTORCOORD, ">$tmpOutputDir/$anno2Label.acceptor.coord.txt";

#transcript num of gene
@geneId = keys(%anno2GeneIdToTrsptNum);
foreach $geneId(@geneId){
	#print TOTALTRANSCRIPTNUM join("\t", $anno2GeneIdToTrsptNum{$geneId}, $anno2Label) . "\n";
	if($anno2GeneIdToTrsptNum{$geneId}<=5){
		$anno2TranscriptNum{$anno2GeneIdToTrsptNum{$geneId}} += 1;
	}else{
		$anno2TranscriptNum{">5"} += 1;
	}
}

#transcript concat exon coord
@transcriptId =();
@transcriptId = keys(%anno2TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){
	print TRANSCRIPTCONCATEXONCOORD $anno2TrsptIdToConcatExonCoord{$transcriptId} . "\n";
}

#transcript len
@transcriptId =();
@transcriptId = keys(%anno2TrsptIdToTrsptLen);
foreach $transcriptId(@transcriptId){
	print TOTALTRANSCRIPTLEN join("\t", $anno2TrsptIdToTrsptLen{$transcriptId}, $anno2Label) . "\n";
}

#exon coordinate
foreach $exonCoord(@anno2ExonCoord){
	print EXONCOORD $exonCoord . "\n";
}

#exon length
#foreach $exonLen(@anno2ExonLen){
#	print TOTALEXONLEN join("\t", $exonLen, $anno2Label) . "\n";
#}

#Donor site
foreach $donorCoord(@anno2DonorCoord){
	print DONORCOORD $donorCoord . "\n";
}

#Acceptor site
foreach $acceptorCoord(@anno2AcceptorCoord){
	print ACCEPTORCOORD $acceptorCoord . "\n";
}

close TRANSCRIPTCONCATEXONCOORD;
close EXONCOORD;
close DONORCOORD;
close ACCEPTORCOORD;


################################
#
# ======== anno3 ============ 
#
###############################
##read anno3Gtf into hashes, arrays
open FF, "<$anno3Gtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	$geneId = &getGeneId($fields[8]);
	$transcriptId = &getTranscriptId($fields[8]);
	if($fields[2] eq "transcript"){

		#The num of gene with this transcript was increased 1
		$anno3GeneIdToTrsptNum{$geneId}++;

	}elsif($fields[2] eq "exon"){

		#concat transcript exon coord
		$anno3TrsptIdToConcatExonCoord{$transcriptId} .= join("___", $fields[0], $fields[3], $fields[4], $fields[6]) . "!";
		#increase transcript length with this exon length
		$anno3TrsptIdToTrsptLen{$transcriptId} += $fields[4] - $fields[3] + 1;
		#regist exon coord into array
		$anno3ExonCoord[$#anno3ExonCoord + 1] = join("___", $fields[0], $fields[3], $fields[4], $fields[6]);
		#regist exon length into array
		$anno3ExonLen[$#anno3ExonLen + 1] = $fields[4] - $fields[3] + 1;
		#register exon num of transcript
		$anno3TrsptIdToExonNum{$transcriptId}++;

	}
}
close FF;

#sort exon according to coordinate and load i
my (@exonCoordinateFullLine, $exonCoordinateFullLine);
@transcriptId = ();
@transcriptId = keys(%anno3TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){

	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/!/, $anno3TrsptIdToConcatExonCoord{$transcriptId});

	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '___', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#register donor
	for(my $i=0; $i<=$#srtExon - 1; $i++){
		$anno3DonorCoord[$#anno3DonorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][2], $srtExon[$i][3]);
	}
	#register acceptor
	for(my $i=1; $i<=$#srtExon; $i++){
		$anno3AcceptorCoord[$#anno3AcceptorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][3]);
	}
}

#output transcriptNum, transcriptConcatExonCoord, transcriptLength, exonCoord, exonLen, donorCoord, acceptorCoord
open TRANSCRIPTCONCATEXONCOORD, ">$tmpOutputDir/$anno3Label.transcript.concat.exon.coord.txt";
open EXONCOORD, ">$tmpOutputDir/$anno3Label.exon.coord.txt";
open DONORCOORD, ">$tmpOutputDir/$anno3Label.donor.coord.txt";
open ACCEPTORCOORD, ">$tmpOutputDir/$anno3Label.acceptor.coord.txt";

#transcript num of gene
@geneId = keys(%anno3GeneIdToTrsptNum);
foreach $geneId(@geneId){
	#print TOTALTRANSCRIPTNUM join("\t", $anno3GeneIdToTrsptNum{$geneId}, $anno3Label) . "\n";
	if($anno3GeneIdToTrsptNum{$geneId}<=5){
		$anno3TranscriptNum{$anno3GeneIdToTrsptNum{$geneId}} += 1;
	}else{
		$anno3TranscriptNum{">5"} += 1;
	}
}

#transcript concat exon coord
@transcriptId =();
@transcriptId = keys(%anno3TrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){
	print TRANSCRIPTCONCATEXONCOORD $anno3TrsptIdToConcatExonCoord{$transcriptId} . "\n";
}

#transcript len
@transcriptId =();
@transcriptId = keys(%anno3TrsptIdToTrsptLen);
foreach $transcriptId(@transcriptId){
	print TOTALTRANSCRIPTLEN join("\t", $anno3TrsptIdToTrsptLen{$transcriptId}, $anno3Label) . "\n";
}

#exon coordinate
foreach $exonCoord(@anno3ExonCoord){
	print EXONCOORD $exonCoord . "\n";
}

#exon length
#foreach $exonLen(@anno3ExonLen){
#	print TOTALEXONLEN join("\t", $exonLen, $anno3Label) . "\n";
#}

#Donor site
foreach $donorCoord(@anno3DonorCoord){
	print DONORCOORD $donorCoord . "\n";
}

#Acceptor site
foreach $acceptorCoord(@anno3AcceptorCoord){
	print ACCEPTORCOORD $acceptorCoord . "\n";
}

close TRANSCRIPTCONCATEXONCOORD;
close EXONCOORD;
close DONORCOORD;
close ACCEPTORCOORD;



################################
#
# ======== merge ============ 
#
###############################
##read mergeGtf into hashes, arrays
open FF, "<$mergedAnnoGtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	$geneId = &getGeneId($fields[8]);
	$transcriptId = &getTranscriptId($fields[8]);
	if($fields[2] eq "transcript"){

		#The num of gene with this transcript was increased 1
		$mergeGeneIdToTrsptNum{$geneId}++;

	}elsif($fields[2] eq "exon"){

		#concat transcript exon coord
		$mergeTrsptIdToConcatExonCoord{$transcriptId} .= join("___", $fields[0], $fields[3], $fields[4], $fields[6]) . "!";
		#increase transcript length with this exon length
		$mergeTrsptIdToTrsptLen{$transcriptId} += $fields[4] - $fields[3] + 1;
		#regist exon coord into array
		$mergeExonCoord[$#mergeExonCoord + 1] = join("___", $fields[0], $fields[3], $fields[4], $fields[6]);
		#regist exon length into array
		$mergeExonLen[$#mergeExonLen + 1] = $fields[4] - $fields[3] + 1;
		#register exon num of transcript
		$mergeTrsptIdToExonNum{$transcriptId}++;

	}
}
close FF;

#sort exon according to coordinate and load i
my (@exonCoordinateFullLine, $exonCoordinateFullLine);
@transcriptId = ();
@transcriptId = keys(%mergeTrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){

	@exonCoordinateFullLine = ();
	@srtExon = ();

	#split exon coordinate full into 1-dim exon array
	@exonCoordinateFullLine = split(/!/, $mergeTrsptIdToConcatExonCoord{$transcriptId});

	#split exon coordinate into 2-dim exon array
	foreach my $exonCoordinate(@exonCoordinateFullLine){
		push @srtExon, [split '___', $exonCoordinate];	
	}

	#sort by coordinate
	if($srtExon[0][3] eq "+"){
		@srtExon = sort{$a->[1]<=>$b->[1]}@srtExon;	
	}else{
		@srtExon = sort{$b->[1]<=>$a->[1]}@srtExon;
	}

	#register donor
	for(my $i=0; $i<=$#srtExon - 1; $i++){
		$mergeDonorCoord[$#mergeDonorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][2], $srtExon[$i][3]);
	}
	#register acceptor
	for(my $i=1; $i<=$#srtExon; $i++){
		$mergeAcceptorCoord[$#mergeAcceptorCoord + 1] = join("___", $srtExon[$i][0], $srtExon[$i][1], $srtExon[$i][3]);
	}
}


########################################################
#
# draw venn for transcript, exon, donor and acceptor
#
########################################################

#output transcriptNum, transcriptConcatExonCoord, transcriptLength, exonCoord, exonLen, donorCoord, acceptorCoord
open TRANSCRIPTCONCATEXONCOORD, ">$tmpOutputDir/$mergeLabel.transcript.concat.exon.coord.txt";
open EXONCOORD, ">$tmpOutputDir/$mergeLabel.exon.coord.txt";
open DONORCOORD, ">$tmpOutputDir/$mergeLabel.donor.coord.txt";
open ACCEPTORCOORD, ">$tmpOutputDir/$mergeLabel.acceptor.coord.txt";

#transcript num of gene
@geneId = keys(%mergeGeneIdToTrsptNum);
foreach $geneId(@geneId){
	if($mergeGeneIdToTrsptNum{$geneId}<=5){
		$mergeTranscriptNum{$mergeGeneIdToTrsptNum{$geneId}} += 1;
	}else{
		$mergeTranscriptNum{">5"} += 1;
	}
}

#transcript concat exon coord
@transcriptId =();
@transcriptId = keys(%mergeTrsptIdToConcatExonCoord);
foreach $transcriptId(@transcriptId){
	print TRANSCRIPTCONCATEXONCOORD $mergeTrsptIdToConcatExonCoord{$transcriptId} . "\n";
}

#transcript len
@transcriptId =();
@transcriptId = keys(%mergeTrsptIdToTrsptLen);
foreach $transcriptId(@transcriptId){
	print TOTALTRANSCRIPTLEN join("\t", $mergeTrsptIdToTrsptLen{$transcriptId}, $mergeLabel ) . "\n";
}

#exon coordinate
foreach $exonCoord(@mergeExonCoord){
	print EXONCOORD $exonCoord . "\n";
}

#exon length
#foreach $exonLen(@mergeExonLen){
#	print TOTALEXONLEN join("\t", $exonLen, $mergeLabel) . "\n";
#}

#Donor site
foreach $donorCoord(@mergeDonorCoord){
	print DONORCOORD $donorCoord . "\n";
}

#Acceptor site
foreach $acceptorCoord(@mergeAcceptorCoord){
	print ACCEPTORCOORD $acceptorCoord . "\n";
}

#close TRANSCRIPTNU;
close TRANSCRIPTCONCATEXONCOORD;
#close TRANSCRIPTLEN;
close EXONCOORD;
#close EXONLEN;
close DONORCOORD;
close ACCEPTORCOORD;






#draw venn for transcript
my $cmd = $Rscript . " " . $drawVennRscript . " " .
	" -A $tmpOutputDir/$anno1Label.transcript.concat.exon.coord.txt" . 
	" -a $anno1Label" .
	" -B $tmpOutputDir/$anno2Label.transcript.concat.exon.coord.txt" .
	" -b $anno2Label" .
        " -C $tmpOutputDir/$anno3Label.transcript.concat.exon.coord.txt" .
        " -c $anno3Label" .
	" -D $tmpOutputDir/$mergeLabel.transcript.concat.exon.coord.txt" .
	" -d $mergeLabel" .
	" -i $outputTrsptEnsRefMergVenn" .
	" -t $outputImagType";
system($cmd);

#draw venn for exon
#remove redundant exon, donor and acceptors
system("sort -u " . $tmpOutputDir . "/" . $anno1Label . ".exon.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno1Label . ".uniq.exon.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno2Label . ".exon.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno2Label . ".uniq.exon.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno3Label . ".exon.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno3Label . ".uniq.exon.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $mergeLabel . ".exon.coord.txt" . " > " .  $tmpOutputDir . "/" . $mergeLabel . ".uniq.exon.coord.txt");
my $cmd = $Rscript . " " . $drawVennRscript . " " .
	" -A $tmpOutputDir/$anno1Label.uniq.exon.coord.txt" . 
	" -a $anno1Label" .
	" -B $tmpOutputDir/$anno2Label.uniq.exon.coord.txt" .
	" -b $anno2Label" .
        " -C $tmpOutputDir/$anno3Label.uniq.exon.coord.txt" .
        " -c $anno3Label" .
	" -D $tmpOutputDir/$mergeLabel.uniq.exon.coord.txt" .
	" -d $mergeLabel" .
	" -i $outputExonEnsRefMergVenn" .
	" -t $outputImagType";
system($cmd);

#draw venn for donor
system("sort -u " . $tmpOutputDir . "/" . $anno1Label . ".donor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno1Label . ".uniq.donor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno2Label . ".donor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno2Label . ".uniq.donor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno3Label . ".donor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno3Label . ".uniq.donor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $mergeLabel . ".donor.coord.txt" . " > " .  $tmpOutputDir . "/" . $mergeLabel . ".uniq.donor.coord.txt");
my $cmd = $Rscript . " " . $drawVennRscript . " " .
	" -A $tmpOutputDir/$anno1Label.uniq.donor.coord.txt" . 
	" -a $anno1Label" .
	" -B $tmpOutputDir/$anno2Label.uniq.donor.coord.txt" .
	" -b $anno2Label" .
	" -C $tmpOutputDir/$anno3Label.uniq.donor.coord.txt" .
	" -c $anno3Label" .
	" -D $tmpOutputDir/$mergeLabel.uniq.donor.coord.txt" .
	" -d $mergeLabel" .
	" -i $outputDonorEnsRefMergVenn" .
	" -t $outputImagType";
system($cmd);

#draw venn for acceptor
system("sort -u " . $tmpOutputDir . "/" . $anno1Label . ".acceptor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno1Label . ".uniq.acceptor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno2Label . ".acceptor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno2Label . ".uniq.acceptor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $anno3Label . ".acceptor.coord.txt" . " > " .  $tmpOutputDir . "/" . $anno3Label . ".uniq.acceptor.coord.txt");
system("sort -u " . $tmpOutputDir . "/" . $mergeLabel . ".acceptor.coord.txt" . " > " .  $tmpOutputDir . "/" . $mergeLabel . ".uniq.acceptor.coord.txt");
my $cmd = $Rscript . " " . $drawVennRscript . " " .
	" -A $tmpOutputDir/$anno1Label.uniq.acceptor.coord.txt" . 
	" -a $anno1Label" .
	" -B $tmpOutputDir/$anno2Label.uniq.acceptor.coord.txt" .
	" -b $anno2Label" .
	" -C $tmpOutputDir/$anno3Label.uniq.acceptor.coord.txt" .
	" -c $anno3Label" .
	" -D $tmpOutputDir/$mergeLabel.uniq.acceptor.coord.txt" .
	" -d $mergeLabel" .
	" -i $outputAcceptorEnsRefMergVenn" .
	" -t $outputImagType";
system($cmd);

##############################################
#
# draw density for exon and transcript length
#
#############################################

#draw density for exon length: $tmpOutputDir/total.exon.length.txt
#awk -F '___' '{print $3-$2+1"\tCombine"}'
system("rm -rf " . $tmpOutputDir . "/total.exon.length.txt");
my $cmd = "awk -F \'___\' \'{print \$3-\$2+1\"\\t$anno1Label\"}\' " . $tmpOutputDir . "/" . $anno1Label . ".uniq.exon.coord.txt" . " > " . $tmpOutputDir . "/total.exon.length.txt";
system($cmd);

$cmd = "awk -F \'___\' \'{print \$3-\$2+1\"\\t$anno2Label\"}\' " . $tmpOutputDir . "/" . $anno2Label . ".uniq.exon.coord.txt" . " >> " . $tmpOutputDir . "/total.exon.length.txt";

$cmd = "awk -F \'___\' \'{print \$3-\$2+1\"\\t$anno3Label\"}\' " . $tmpOutputDir . "/" . $anno3Label . ".uniq.exon.coord.txt" . " >> " . $tmpOutputDir . "/total.exon.length.txt";
system($cmd);

$cmd = "awk -F \'___\' \'{print \$3-\$2+1\"\\t$mergeLabel\"}\' " . $tmpOutputDir . "/" . $mergeLabel . ".uniq.exon.coord.txt" . " >> " . $tmpOutputDir . "/total.exon.length.txt";
system($cmd);

#xvfb-run Rscript --no-save ~/software/customized/draw.densityDiagram.for.three.sets.R -d total.exon.length.txt -X "Exon len" -m 1 -x 200 -i exon.dentisity.png -t png

my $cmd = $xvfb_run . " " . $Rscript . " --no-save " . $drawDensityRscript . " -d " . "$tmpOutputDir/total.exon.length.txt " . " -X \"Exon length\" " . " -m $densityExonMinLen -x $densityExonMaxLen " . " -i " . $outputExonLenEnsRefMergDensity . " -t " . $outputImagType;
#print $cmd . "\n";
system($cmd);

#draw density for transcript length
my $cmd = $xvfb_run . " " . $Rscript . " --no-save " . $drawDensityRscript . " -d " . "$tmpOutputDir/total.transcript.length.txt" . " -X \"Exon length\" " . " -m $densityTranscriptMinLen -x $densityTranscriptMaxLen " . " -i " . $outputTrsptLenEnsRefMergDensity . " -t " . $outputImagType;
#print $cmd . "\n";
system($cmd);

###################################################
#
#  draw bar for gene with different transcript num
#
####################################################

#draw bar for percentage of transcript num of gene
open TOTALTRANSCRIPTNUM, ">$tmpOutputDir/total.transcript.num.txt";
#print TOTALTRANSCRIPTNUM join("\t", "TrsptNum", "Annotation") . "\n";
my @transcriptNum = ();
@transcriptNum = keys(%anno1TranscriptNum);
my $totalGeneNum = 0;
my $nnn = 0;
#calculate total gene num
foreach my $transcriptNum(@transcriptNum){
	$totalGeneNum += $anno1TranscriptNum{$transcriptNum};
}
#adjust and output each class gene in percentage format
foreach my $transcriptNum(@transcriptNum){
	my $outputGeneNum = int($anno1TranscriptNum{$transcriptNum}/$totalGeneNum*1000+0.5);
	for($nnn=0; $nnn<$outputGeneNum;  $nnn+=1){
		print TOTALTRANSCRIPTNUM join("\t", $transcriptNum, $anno1Label) . "\n";
	}
}

@transcriptNum = ();
@transcriptNum = keys(%anno2TranscriptNum);
$totalGeneNum = 0;
#calculate total gene num
foreach my $transcriptNum(@transcriptNum){
	$totalGeneNum += $anno2TranscriptNum{$transcriptNum};
}
#adjust and output each class gene in percentage format
foreach my $transcriptNum(@transcriptNum){
	my $outputGeneNum = int($anno2TranscriptNum{$transcriptNum}/$totalGeneNum*1000+0.5);
	for($nnn=0; $nnn<$outputGeneNum;  $nnn+=1){
		print TOTALTRANSCRIPTNUM join("\t", $transcriptNum, $anno2Label) . "\n";
	}
}

@transcriptNum = ();
@transcriptNum = keys(%anno3TranscriptNum);
$totalGeneNum = 0;
#calculate total gene num
foreach my $transcriptNum(@transcriptNum){
	$totalGeneNum += $anno3TranscriptNum{$transcriptNum};
}
#adjust and output each class gene in percentage format
foreach my $transcriptNum(@transcriptNum){
	my $outputGeneNum = int($anno3TranscriptNum{$transcriptNum}/$totalGeneNum*1000+0.5);
	for($nnn=0; $nnn<$outputGeneNum;  $nnn+=1){
		print TOTALTRANSCRIPTNUM join("\t", $transcriptNum, $anno3Label) . "\n";
	}
}


@transcriptNum = ();
@transcriptNum = keys(%mergeTranscriptNum);
$totalGeneNum = 0;
#calculate total gene num
foreach my $transcriptNum(@transcriptNum){
	$totalGeneNum += $mergeTranscriptNum{$transcriptNum};
}
#adjust and output each class gene in percentage format
foreach my $transcriptNum(@transcriptNum){
	my $outputGeneNum = int($mergeTranscriptNum{$transcriptNum}/$totalGeneNum*1000+0.5);
	for(my $nnn=0; $nnn<$outputGeneNum;  $nnn+=1){
		print TOTALTRANSCRIPTNUM join("\t", $transcriptNum, $mergeLabel) . "\n";
	}
}
close TOTALTRANSCRIPTNUM;

#draw bar for genes with different truanscript num
# xvfb-run Rscript --no-save  ~/software/customized/draw.barDiagram.for.three.sets.R -d total.transcript.num.txt -X Annotation -Y "Permill of gene with various transcript num" -L "Gene type with different transcript num" -i transcript.num.png -t png
my $cmd = $xvfb_run . " " . $Rscript . " --no-save  " . $drawBarRscript .  " -d " . $tmpOutputDir . "/total.transcript.num.txt " . " -X Annotation -Y \"Permill of gene with various transcript num\" -L \"Gene type with different transcript num\" -i $outputTrsptNumEnsRefMergBar -t png";
system($cmd);

####################################################
#
# draw bar for transcript with different exon num
#
###################################################

open TOTALEXONNUM, ">$tmpOutputDir/total.exon.num.txt";

#anno1
#register transcript num into the exon num hash table
my %exonNumTypeToTrspIdNum = ();
@transcriptId = ();
@transcriptId = keys(%anno1TrsptIdToExonNum);
foreach $transcriptId(@transcriptId){
	if($anno1TrsptIdToExonNum{$transcriptId}<=20){
		$exonNumTypeToTrspIdNum{$anno1TrsptIdToExonNum{$transcriptId}}+=1;
	}else{
		$exonNumTypeToTrspIdNum{">20"}+=1;
	}
}
#calculate total transcript num
my $totalTrsptNum = 0;
$nnn = 0;
my @exonNumType = ();
@exonNumType = keys(%exonNumTypeToTrspIdNum);
#calculate total gene num
foreach my $exonNumType(@exonNumType){
	$totalTrsptNum += $exonNumTypeToTrspIdNum{$exonNumType};
}
#adjust and output each class gene in percentage format
foreach my $exonNumType(@exonNumType){
	my $outputTrsptNum = int($exonNumTypeToTrspIdNum{$exonNumType}/$totalTrsptNum*1000+0.5);
	for($nnn=0; $nnn<$outputTrsptNum;  $nnn+=1){
		print TOTALEXONNUM join("\t", $exonNumType, $anno1Label) . "\n";
	}
}

#anno2
%exonNumTypeToTrspIdNum = ();
@transcriptId = ();
@transcriptId = keys(%anno2TrsptIdToExonNum);
foreach $transcriptId(@transcriptId){
	if($anno2TrsptIdToExonNum{$transcriptId}<=20){
		$exonNumTypeToTrspIdNum{$anno2TrsptIdToExonNum{$transcriptId}}+=1;
	}else{
		$exonNumTypeToTrspIdNum{">20"}+=1;
	}
}
#calculate total transcript num
$totalTrsptNum = 0;
$nnn = 0;
@exonNumType = ();
@exonNumType = keys(%exonNumTypeToTrspIdNum);
#calculate total gene num
foreach my $exonNumType(@exonNumType){
	$totalTrsptNum += $exonNumTypeToTrspIdNum{$exonNumType};
}
#adjust and output each class gene in percentage format
foreach my $exonNumType(@exonNumType){
	my $outputTrsptNum = int($exonNumTypeToTrspIdNum{$exonNumType}/$totalTrsptNum*1000+0.5);
	for($nnn=0; $nnn<$outputTrsptNum;  $nnn+=1){
		print TOTALEXONNUM join("\t", $exonNumType, $anno2Label) . "\n";
	}
}

#anno3
%exonNumTypeToTrspIdNum = ();
@transcriptId = ();
@transcriptId = keys(%anno3TrsptIdToExonNum);
foreach $transcriptId(@transcriptId){
	if($anno3TrsptIdToExonNum{$transcriptId}<=20){
		$exonNumTypeToTrspIdNum{$anno3TrsptIdToExonNum{$transcriptId}}+=1;
	}else{
		$exonNumTypeToTrspIdNum{">20"}+=1;
	}
}
#calculate total transcript num
$totalTrsptNum = 0;
$nnn = 0;
@exonNumType = ();
@exonNumType = keys(%exonNumTypeToTrspIdNum);
#calculate total gene num
foreach my $exonNumType(@exonNumType){
	$totalTrsptNum += $exonNumTypeToTrspIdNum{$exonNumType};
}
#adjust and output each class gene in percentage format
foreach my $exonNumType(@exonNumType){
	my $outputTrsptNum = int($exonNumTypeToTrspIdNum{$exonNumType}/$totalTrsptNum*1000+0.5);
	for($nnn=0; $nnn<$outputTrsptNum;  $nnn+=1){
		print TOTALEXONNUM join("\t", $exonNumType, $anno3Label) . "\n";
	}
}


#merge
%exonNumTypeToTrspIdNum = ();
@transcriptId = ();
@transcriptId = keys(%mergeTrsptIdToExonNum);
foreach $transcriptId(@transcriptId){
	if($mergeTrsptIdToExonNum{$transcriptId}<=20){
		$exonNumTypeToTrspIdNum{$mergeTrsptIdToExonNum{$transcriptId}}+=1;
	}else{
		$exonNumTypeToTrspIdNum{">20"}+=1;
	}
}
#calculate total transcript num
my $totalTrsptNum = 0;
$nnn = 0;
@exonNumType = ();
@exonNumType = keys(%exonNumTypeToTrspIdNum);
#calculate total gene num
foreach my $exonNumType(@exonNumType){
	$totalTrsptNum += $exonNumTypeToTrspIdNum{$exonNumType};
}
#adjust and output each class gene in percentage format
foreach my $exonNumType(@exonNumType){
	my $outputTrsptNum = int($exonNumTypeToTrspIdNum{$exonNumType}/$totalTrsptNum*1000+0.5);
	for($nnn=0; $nnn<$outputTrsptNum;  $nnn+=1){
		print TOTALEXONNUM join("\t", $exonNumType, $mergeLabel) . "\n";
	}
}
close TOTALEXONNUM;

my $cmd = $xvfb_run . " " . $Rscript . " --no-save  " . $drawBarRscript .  " -d " . $tmpOutputDir . "/total.exon.num.txt " . " -X Annotation -Y \"Permill of transcript with different exon num\" -L \"Transcript with different exon num\" -i $outputExonNumEnsRefMergBar -t png";
system($cmd);

#generate statistic table
my ($geneNum, $transcriptNum, $exonNum);
my (@exon, $donorNum, $acceptorNum);
my ($totalLen, $avgTrsptLen, $avgExonLen, $donorNum);
open STA, ">$outputStaTbl";

#anno1
#statistic gene num
print STA join("\t", "Annotation", "Gene_num", "Trspt_num", "Exon_num", "Trspt_len(avg)", "Exon_len(avg)", "Donor_num", "Acceptor_num") . "\n";
@geneId = keys(%anno1GeneIdToTrsptNum);
$geneNum = $#geneId + 1;
@transcriptId = ();
@transcriptId = keys(%anno1TrsptIdToConcatExonCoord);
$transcriptNum = $#transcriptId+1;

#exonNum
$exonNum = "cat " . $tmpOutputDir . "/" . $anno1Label . ".uniq.exon.coord.txt | wc -l";
$exonNum = `$exonNum`;
chomp($exonNum);

#calculate avg transcript length
@transcriptId = ();
@transcriptId = keys(%anno1TrsptIdToTrsptLen);
$totalLen = 0;
foreach $transcriptId(@transcriptId){
	$totalLen+=$anno1TrsptIdToTrsptLen{$transcriptId};
}
$avgTrsptLen = $totalLen/($#transcriptId+1);

#calculate ave exon length
$totalLen = 0;
my $uniqExonNum = 0;
open FEXON, "<" . $tmpOutputDir . "/" . $anno1Label . ".uniq.exon.coord.txt";
while(my $line=<FEXON>){
#1___350267___350389___-
	my @tttt = ();
	@tttt = split(/___/, $line);
	$totalLen+=$tttt[2] - $tttt[1] + 1;
	$uniqExonNum+=1;
}
close FEXON;
$avgExonLen = $totalLen/$uniqExonNum;

#donor num
$donorNum = "cat " . $tmpOutputDir . "/" . $anno1Label . ".uniq.donor.coord.txt | wc -l";
$donorNum = `$donorNum`;
chomp($donorNum);

#acceptor num
$acceptorNum = "cat " . $tmpOutputDir . "/" . $anno1Label . ".uniq.acceptor.coord.txt | wc -l";
$acceptorNum = `$acceptorNum`;
chomp($acceptorNum);

print STA join("\t", $anno1Label, $geneNum, $transcriptNum, $exonNum, sprintf("%.1f",$avgTrsptLen), sprintf("%.1f", $avgExonLen),  $donorNum, $acceptorNum) . "\n";

#anno2
#statistic gene num
@geneId = keys(%anno2GeneIdToTrsptNum);
$geneNum = $#geneId + 1;
@transcriptId = ();
@transcriptId = keys(%anno2TrsptIdToConcatExonCoord);
$transcriptNum = $#transcriptId+1;

#exonNum
$exonNum = "cat " . $tmpOutputDir . "/" . $anno2Label . ".uniq.exon.coord.txt | wc -l";
$exonNum = `$exonNum`;
chomp($exonNum);

#calculate avg transcript length
@transcriptId = ();
@transcriptId = keys(%anno2TrsptIdToTrsptLen);
$totalLen = 0;
foreach $transcriptId(@transcriptId){
	$totalLen+=$anno2TrsptIdToTrsptLen{$transcriptId};
}
$avgTrsptLen = $totalLen/($#transcriptId+1);

#calculate ave exon length
$totalLen = 0;
my $uniqExonNum = 0;
open FEXON, "<" . $tmpOutputDir . "/" . $anno2Label . ".uniq.exon.coord.txt";
while(my $line=<FEXON>){
#1___350267___350389___-
	my @tttt = ();
	@tttt = split(/___/, $line);
	$totalLen+=$tttt[2] - $tttt[1] + 1;
	$uniqExonNum+=1;
}
close FEXON;
$avgExonLen = $totalLen/$uniqExonNum;


#donor num
$donorNum = "cat " . $tmpOutputDir . "/" . $anno2Label . ".uniq.donor.coord.txt | wc -l";
$donorNum = `$donorNum`;
chomp($donorNum);

#acceptor num
$acceptorNum = "cat " . $tmpOutputDir . "/" . $anno2Label . ".uniq.acceptor.coord.txt | wc -l";
$acceptorNum = `$acceptorNum`;
chomp($acceptorNum);

print STA join("\t", $anno2Label, $geneNum, $transcriptNum, $exonNum, sprintf("%.1f",$avgTrsptLen), sprintf("%.1f", $avgExonLen),  $donorNum, $acceptorNum) . "\n";

#anno3
#statistic gene num
@geneId = keys(%anno3GeneIdToTrsptNum);
$geneNum = $#geneId + 1;
@transcriptId = ();
@transcriptId = keys(%anno3TrsptIdToConcatExonCoord);
$transcriptNum = $#transcriptId+1;

#exonNum
$exonNum = "cat " . $tmpOutputDir . "/" . $anno3Label . ".uniq.exon.coord.txt | wc -l";
$exonNum = `$exonNum`;
chomp($exonNum);

#calculate avg transcript length
@transcriptId = ();
@transcriptId = keys(%anno3TrsptIdToTrsptLen);
$totalLen = 0;
foreach $transcriptId(@transcriptId){
	$totalLen+=$anno3TrsptIdToTrsptLen{$transcriptId};
}
$avgTrsptLen = $totalLen/($#transcriptId+1);

#calculate ave exon length
$totalLen = 0;
my $uniqExonNum = 0;
open FEXON, "<" . $tmpOutputDir . "/" . $anno3Label . ".uniq.exon.coord.txt";
while(my $line=<FEXON>){
#1___350267___350389___-
	my @tttt = ();
	@tttt = split(/___/, $line);
	$totalLen+=$tttt[2] - $tttt[1] + 1;
	$uniqExonNum+=1;
}
close FEXON;
$avgExonLen = $totalLen/$uniqExonNum;


#donor num
$donorNum = "cat " . $tmpOutputDir . "/" . $anno3Label . ".uniq.donor.coord.txt | wc -l";
$donorNum = `$donorNum`;
chomp($donorNum);

#acceptor num
$acceptorNum = "cat " . $tmpOutputDir . "/" . $anno3Label . ".uniq.acceptor.coord.txt | wc -l";
$acceptorNum = `$acceptorNum`;
chomp($acceptorNum);

print STA join("\t", $anno3Label, $geneNum, $transcriptNum, $exonNum, sprintf("%.1f",$avgTrsptLen), sprintf("%.1f", $avgExonLen),  $donorNum, $acceptorNum) . "\n";


#merge
#statistic gene num
@geneId = keys(%mergeGeneIdToTrsptNum);
$geneNum = $#geneId + 1;
@transcriptId = ();
@transcriptId = keys(%mergeTrsptIdToConcatExonCoord);
$transcriptNum = $#transcriptId+1;

#exonNum
$exonNum = "cat " . $tmpOutputDir . "/" . $mergeLabel . ".uniq.exon.coord.txt | wc -l";
$exonNum = `$exonNum`;
chomp($exonNum);


#calculate avg transcript length
@transcriptId = ();
@transcriptId = keys(%mergeTrsptIdToTrsptLen);
$totalLen = 0;
foreach $transcriptId(@transcriptId){
	$totalLen+=$mergeTrsptIdToTrsptLen{$transcriptId};
}
$avgTrsptLen = $totalLen/($#transcriptId+1);

#calculate ave exon length
$totalLen = 0;
$uniqExonNum = 0;
open FEXON, "<" . $tmpOutputDir . "/" . $mergeLabel . ".uniq.exon.coord.txt";
while(my $line=<FEXON>){
#1___350267___350389___-
	my @tttt = ();
	@tttt = split(/___/, $line);
	$totalLen+=$tttt[2] - $tttt[1] + 1;
	$uniqExonNum+=1;
}
close FEXON;
$avgExonLen = $totalLen/$uniqExonNum;

#donor num
$donorNum = "cat " . $tmpOutputDir . "/" . $mergeLabel . ".uniq.donor.coord.txt | wc -l";
$donorNum = `$donorNum`;
chomp($donorNum);

#acceptor num
$acceptorNum = "cat " . $tmpOutputDir . "/" . $mergeLabel . ".uniq.acceptor.coord.txt | wc -l";
$acceptorNum = `$acceptorNum`;
chomp($acceptorNum);
print STA join("\t", $mergeLabel, $geneNum, $transcriptNum, $exonNum, sprintf("%.1f",$avgTrsptLen), sprintf("%.1f", $avgExonLen), $donorNum, $acceptorNum) . "\n";
close STA;



#system("rm -rf $tmpOutputDir");


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
