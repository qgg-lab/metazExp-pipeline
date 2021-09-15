#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--entrezMapFile \\\n" .
                "--refseqMapFile \\\n" .
                "--uniprotMapFile \\\n" .
		"--outputGeneIdMap \\\n" .
		"--outputTrsptIdMap \n";
	exit;
}

my ($entrezMapFile, $refseqMapFile, $uniprotMapFile, $outputGeneIdMap, $outputTrsptIdMap);

GetOptions(
        'entrezMapFile=s'=>\$entrezMapFile,
        'refseqMapFile=s'=>\$refseqMapFile,
        'uniprotMapFile=s'=>\$uniprotMapFile,
        'outputGeneIdMap=s'=>\$outputGeneIdMap,
	'outputTrsptIdMap=s'=>\$outputTrsptIdMap,
);

my (%geneToEntrezIdList, %geneToRefseqIdList, %geneToUniprotIdList);
my (%trsptToEntrezIdList, %trsptToRefseqIdList, %trsptToUniprotIdList);
my (@geneId, $geneId, @trsptId, $trsptId, $proteinId, $line, @field, $field);
my ($entrezId, $refseqId, $uniprotId);
my (%geneAnno, $geneAnnoHref, %trsptAnno, $trsptAnnoHref);
$geneAnnoHref = \%geneAnno;
$trsptAnnoHref =\%trsptAnno;

open FF, "<$entrezMapFile";
<FF>;
# gene_stable_id  transcript_stable_id    protein_stable_id       xref    db_name    info_type       source_identity xref_identity   linkage_type
# Os01g0100100    Os01t0100100-01         Os01t0100100-01         4326813 EntrezGene DEPENDENT       -       -       -
# Os01g0100600    Os01t0100600-01         Os01t0100600-01         4326456 EntrezGene DEPENDENT       -       -       -
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$geneId = $field[0];
	$trsptId = $field[1];
	$entrezId = $field[3];
	if(not exists($geneToEntrezIdList{$geneId})){
		$geneToEntrezIdList{$geneId} = "#$entrezId#";
	}elsif(index($geneToEntrezIdList{$geneId}, "#$entrezId#")<0){
		$geneToEntrezIdList{$geneId} .= "$entrezId#";
	}
	if(not exists($trsptToEntrezIdList{$trsptId})){
		$trsptToEntrezIdList{$trsptId} = "#$entrezId#";
	}elsif(index($trsptToEntrezIdList{$trsptId}, "#$entrezId#")<0){
		$trsptToEntrezIdList{$trsptId} .= "$entrezId#";
	}

}
close FF;
@geneId = ();
@geneId = keys(%geneToEntrezIdList);
foreach $geneId(@geneId){
	$geneAnnoHref->{$geneId}->{"EntrezIdList"} = $geneToEntrezIdList{$geneId};
}
@trsptId = ();
@trsptId = keys(%trsptToEntrezIdList);
foreach $trsptId(@trsptId){
	$trsptAnnoHref->{$trsptId}->{"EntrezIdList"} = $trsptToEntrezIdList{$trsptId};
}




open FF, "<$refseqMapFile";
#gene_stable_id  transcript_stable_id    protein_stable_id       xref    db_name info_type       source_identity xref_identity   linkage_type
#Os01g0100100    Os01t0100100-01 Os01t0100100-01 XP_015622096.1  RefSeq_peptide  DEPENDENT       -       -       -
#Os01g0100100    Os01t0100100-01 Os01t0100100-01 XM_015766610.1  RefSeq_dna      DEPENDENT       -       -       -
<FF>;
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$geneId = $field[0];
	$trsptId = $field[1];
	$refseqId = $field[3];
	if(not exists($geneToRefseqIdList{$geneId})){
		$geneToRefseqIdList{$geneId} = "#$refseqId#";
	}elsif(index($geneToRefseqIdList{$geneId}, "#$refseqId#")<0){
		$geneToRefseqIdList{$geneId} .= "$refseqId#";
	}
	if(not exists($trsptToRefseqIdList{$trsptId})){
		$trsptToRefseqIdList{$trsptId} = "#$refseqId#";
	}elsif(index($trsptToRefseqIdList{$trsptId}, "#$refseqId#")<0){
		$trsptToRefseqIdList{$trsptId} .= "$refseqId#";
	}

}
close FF;

@geneId = ();
@geneId = keys(%geneToRefseqIdList);
foreach $geneId(@geneId){
	$geneAnnoHref->{$geneId}->{"RefseqIdList"} = $geneToRefseqIdList{$geneId};
}
@trsptId = ();
@trsptId = keys(%trsptToRefseqIdList);
foreach $trsptId(@trsptId){
	$trsptAnnoHref->{$trsptId}->{"RefseqIdList"} = $trsptToRefseqIdList{$trsptId};
}



open FF, "<$uniprotMapFile";
# gene_stable_id  transcript_stable_id    protein_stable_id       xref    db_name info_type       source_identity xref_identity   linkage_type
# Os01g0100100    Os01t0100100-01 Os01t0100100-01 A0A0P0UX28      Uniprot/SPTREMBL        DEPENDENT       -       -       -
# Os01g0100200    Os01t0100200-01 Os01t0100200-01 Q655L9  Uniprot/SPTREMBL        DEPENDENT       -      
<FF>;
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$geneId = $field[0];
	$trsptId = $field[1];
	$uniprotId = $field[3];
	if(not exists($geneToUniprotIdList{$geneId})){
		$geneToUniprotIdList{$geneId} = "#$uniprotId#";
	}elsif(index($geneToUniprotIdList{$geneId}, "#$uniprotId#")<0){
		$geneToUniprotIdList{$geneId} .= "$uniprotId#";
	}
	if(not exists($trsptToUniprotIdList{$trsptId})){
		$trsptToUniprotIdList{$trsptId} = "#$uniprotId#";
	}elsif(index($trsptToUniprotIdList{$trsptId}, "#$uniprotId#")<0){
		$trsptToUniprotIdList{$trsptId} .= "$uniprotId#";
	}

}
close FF;

@geneId = ();
@geneId = keys(%geneToUniprotIdList);
foreach $geneId(@geneId){
	$geneAnnoHref->{$geneId}->{"UniprotIdList"} = $geneToUniprotIdList{$geneId};
}
@trsptId = ();
@trsptId = keys(%trsptToUniprotIdList);
foreach $trsptId(@trsptId){
	$trsptAnnoHref->{$trsptId}->{"UniprotIdList"} = $trsptToUniprotIdList{$trsptId};
}


my ($EntrezIdList, $RefseqIdList, $UniprotIdList);
open WW, ">$outputGeneIdMap";
print WW join("\t", "geneId", "EntrezIdList", "RefseqIdList", "UniprotIdList") . "\n";
@geneId = ();
@geneId = keys(%geneAnno);
foreach $geneId(@geneId){
	if(exists($geneAnnoHref->{$geneId}->{"EntrezIdList"})){
		$EntrezIdList=$geneAnnoHref->{$geneId}->{"EntrezIdList"};
		$EntrezIdList=substr($EntrezIdList, 1, length($EntrezIdList) -2);
		$EntrezIdList=~s/#/,/g;
	}else{
		$EntrezIdList="NA";
	}

	if(exists($geneAnnoHref->{$geneId}->{"RefseqIdList"})){
		$RefseqIdList=$geneAnnoHref->{$geneId}->{"RefseqIdList"};
		$RefseqIdList=substr($RefseqIdList, 1, length($RefseqIdList) -2);
		$RefseqIdList=~s/#/,/g;
	}else{
		$RefseqIdList="NA";
	}

	if(exists($geneAnnoHref->{$geneId}->{"UniprotIdList"})){
		$UniprotIdList=$geneAnnoHref->{$geneId}->{"UniprotIdList"};
		$UniprotIdList=substr($UniprotIdList, 1, length($UniprotIdList) -2);
		$UniprotIdList=~s/#/,/g;
	}else{
		$UniprotIdList="NA";
	}

	print WW join("\t", $geneId,  $EntrezIdList, $RefseqIdList, $UniprotIdList) . "\n";
}
close WW;

my ($EntrezIdList, $RefseqIdList, $UniprotIdList);
open WW, ">$outputTrsptIdMap";
print WW join("\t", "trsptId", "EntrezIdList", "RefseqIdList", "UniprotIdList") . "\n";
@trsptId = ();
@trsptId = keys(%trsptAnno);
foreach $trsptId(@trsptId){
	if(exists($trsptAnnoHref->{$trsptId}->{"EntrezIdList"})){
		$EntrezIdList=$trsptAnnoHref->{$trsptId}->{"EntrezIdList"};
		$EntrezIdList=substr($EntrezIdList, 1, length($EntrezIdList) -2);
		$EntrezIdList=~s/#/,/g;
	}else{
		$EntrezIdList="NA";
	}

	if(exists($trsptAnnoHref->{$trsptId}->{"RefseqIdList"})){
		$RefseqIdList=$trsptAnnoHref->{$trsptId}->{"RefseqIdList"};
		$RefseqIdList=substr($RefseqIdList, 1, length($RefseqIdList) -2);
		$RefseqIdList=~s/#/,/g;
	}else{
		$RefseqIdList="NA";
	}

	if(exists($trsptAnnoHref->{$trsptId}->{"UniprotIdList"})){
		$UniprotIdList=$trsptAnnoHref->{$trsptId}->{"UniprotIdList"};
		$UniprotIdList=substr($UniprotIdList, 1, length($UniprotIdList) -2);
		$UniprotIdList=~s/#/,/g;
	}else{
		$UniprotIdList="NA";
	}

	print WW join("\t", $trsptId, $EntrezIdList, $RefseqIdList, $UniprotIdList) . "\n";
}
close WW;
