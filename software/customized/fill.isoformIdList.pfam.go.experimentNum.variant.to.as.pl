#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--species \"Bos taurus\"\\\n" . 
		"--speciesAbbr BTAU \\\n" .
		"--jcecExpNumFile jcec.experimentNum.txt \\\n".
		"--jcExpNumFile jc.experimentNum.txt \\\n".
		"--genomeFile genome.fa \\\n" .
		"--isoformTsvFile 9925.isoform.tsv \\\n" .
		"--variantTsvFile 9925.variant.tsv \\\n" .
		"--asTsvFile 9925.as.tsv \\\n" .
		"--tmpDir ./tmpDir \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($jcecExpNumFile, $jcExpNumFile);
my ($species, $speciesAbbr, $genomeFile);
my ($isoformTsvFile, $variantTsvFile, $asTsvFile);
my ($tmpDir);
my ($dbh, $query, $i);

GetOptions(
        'species=s'=>\$species,
	'speciesAbbr=s'=>\$speciesAbbr,
	'jcecExpNumFile=s'=>\$jcecExpNumFile,
	'jcExpNumFile=s'=>\$jcExpNumFile,
	'genomeFile=s'=>\$genomeFile,
	'isoformTsvFile=s'=>\$isoformTsvFile,
	'variantTsvFile=s'=>\$variantTsvFile,
	'asTsvFile=s'=>\$asTsvFile,
	'tmpDir=s'=>\$tmpDir,
);

my ($line, $fieldNameString, @fieldName, $fieldName, $valueString, @value, $value);
my (@tmp);

###### --- 1 --- ###########
# load isoforms into hash to obtain pfam, go, exonSeries and annoOrigin to gene
my (%isoform, %isoformOrigin, %geneToIsoform);
open FF, "<$isoformTsvFile";
while($line=<FF>){
        chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

        $valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);

	%isoform = ();

        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $isoform{$fieldName[$i]}=$value[$i];
        }

	$geneToIsoform{$isoform{"geneId"}} .= $isoform{"isoformId"} . "#" . $isoform{"exonSeries"} . "#" . $isoform{"isoformPfam"} . "#" . $isoform{"isoformGo"} . "\n";
	$isoformOrigin{$isoform{"isoformId"}} = $isoform{"annoOrigin"};
}
close FF;




my (%asExperimentNum);
print &getLocalTime();
print ": Begin load jcecExpnum of AS into hash.\n";
open FF, "<$jcecExpNumFile";
#313 ECABA3SS0000000003
while($line = <FF>){
	if($line=~/^ *(\d+) +($speciesAbbr.*)\n/){
		${$asExperimentNum{$2}}{"jcecExpNum"} = $1;
	}
}
close FF;

print &getLocalTime();
print ": Finish load jcecExpnum of AS into hash.\n";

print &getLocalTime();
print ": Begin load jcExpnum of AS into hash.\n";
open FF, "<$jcExpNumFile";
# 45 ECABA3SS0000000001
while($line = <FF>){
	if($line=~/^ *(\d+) +($speciesAbbr.*)\n/){
		${$asExperimentNum{$2}}{"jcExpNum"} = $1;
	}
}
close FF;

print &getLocalTime();
print ": Finish load jcExpnum of AS into hash.\n";





my (%genomeSeq, $line, @tt, $id);
print &getLocalTime();
print ": Begin load genome sequence into hash.\n";
open FF, "<$genomeFile";
while($line = <FF>){
        chomp($line);
        if($line=~/>/){
                @tt = ();
                @tt = split(/ /, $line);
                if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
                        $id = $1;
                }
        }else{
                $genomeSeq{$id} .=uc($line);
        }
}
close FF;

print &getLocalTime();
print ": Finish load genome sequence into hash.\n";




print &getLocalTime();
print ": Begin load AS into hash and mask AS regeions in genome sequences with spaces.\n";
my (%cmplAs, $asId, $sblen, $i);
open FF, "<$asTsvFile";
while($line=<FF>){
        chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

        $valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);
	
	$asId = $value[0];
	if($asId=~/"(.*)"/){
		$asId = $1;
	}
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }

		${$cmplAs{$asId}}{$fieldName[$i]} = $value[$i];

        }

        $sblen = ${$cmplAs{$asId}}{"end"} - ${$cmplAs{$asId}}{"start"} + 1;
        substr($genomeSeq{${$cmplAs{$asId}}{"chr"}}, ${$cmplAs{$asId}}{"start"}-1, $sblen) = " "x$sblen;
}
close FF;

print &getLocalTime();
print ": Finish load AS into hash and mask AS regeions in genome sequences with spaces.\n";


print &getLocalTime();
print ": Begin flag variants in AS region in genome sequences.\n";
my (%tmpVar);
open FF, "<$variantTsvFile";
while($line=<FF>){
        chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

        $valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);

	%tmpVar = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $tmpVar{$fieldName[$i]}=$value[$i];
        }
        if($tmpVar{"varType"} eq "SNV"){
                substr($genomeSeq{$tmpVar{"chr"}}, $tmpVar{"pos"}-1, 1) = "S";
        }elsif($tmpVar{"varType"} eq "deletion"){
                substr($genomeSeq{$tmpVar{"chr"}}, $tmpVar{"pos"}-1, 1) = "D";
        }elsif($tmpVar{"varType"} eq "insertion"){
                substr($genomeSeq{$tmpVar{"chr"}}, $tmpVar{"pos"}-1, 1) = "I";
        }
	
}
close FF;

print &getLocalTime();
print ": Finish flag variants in AS region in genome sequences.\n";



print &getLocalTime();
print ": Begin calculate cmbVarType of AS.\n";
my (%asCmbVarType, $asRegion, @asId, $asId);
@asId = ();
@asId = keys(%cmplAs);

foreach $asId(@asId){

        $sblen = ${$cmplAs{$asId}}{"end"} - ${$cmplAs{$asId}}{"start"} + 1;
	$asRegion = substr($genomeSeq{${$cmplAs{$asId}}{"chr"}}, ${$cmplAs{$asId}}{"start"}-1, $sblen);

	$asRegion=~s/ //g;

        if($asRegion=~/S/ and $asRegion=~/D/ and $asRegion=~/I/){
                ${$asCmbVarType{$asId}}{"type"} = "SDI";
        }elsif($asRegion=~/S/ and $asRegion=~/D/){
		${$asCmbVarType{$asId}}{"type"} = "SD";
        }elsif($asRegion=~/S/ and $asRegion=~/I/){
		${$asCmbVarType{$asId}}{"type"} = "SI";
        }elsif($asRegion=~/D/ and $asRegion=~/I/){
		${$asCmbVarType{$asId}}{"type"} = "DI";
        }elsif($asRegion=~/S/){
		${$asCmbVarType{$asId}}{"type"} = "S";
        }elsif($asRegion=~/D/){
		${$asCmbVarType{$asId}}{"type"} = "D";
        }elsif($asRegion=~/I/){
		${$asCmbVarType{$asId}}{"type"} = "I";
        }else{
		${$asCmbVarType{$asId}}{"type"} = "None";
        }
	
        if( ${$asCmbVarType{$asId}}{"type"} eq "None"){
		${$asCmbVarType{$asId}}{"num"} = 0;
        }else{
		${$asCmbVarType{$asId}}{"num"} = length($asRegion);
        }
}

print &getLocalTime();
print ": Finish calculate cmbVarType of AS:\n";




print &getLocalTime();
print ": Assign pfam, go, cmbVarType, varNum, jcecExpNum, jcExpNum, descoverApproach to each AS.\n";

open WW, ">$tmpDir" . "/tmp.as.with.complete.infor.tsv";
my ($firstIsoformIdList, $firstIsoformNum, $secondIsoformIdList, $secondIsoformNum, $pfam, $go);
my ($variantPointTypeCmb, $variantPointNum, $discoveryApproach, $jcecExperimentNum, $jcExperimentNum);

foreach $asId(@asId){

	# calculate firstIsoformIdList, firstIsoformNum, $secondIsoformIdList, $secondIsoformNum, pfam, go
	($firstIsoformIdList, $firstIsoformNum, $secondIsoformIdList, $secondIsoformNum, $pfam, $go) = ("", "", "", "", "", "");
	&getIsoformIdListPfamGo( $geneToIsoform{${$cmplAs{$asId}}{"geneID"}}, ${$cmplAs{$asId}}{"firstAltExonSeries"}, ${$cmplAs{$asId}}{"secondAltExonSeries"}, \$firstIsoformIdList, \$firstIsoformNum, \$secondIsoformIdList, \$secondIsoformNum, \$pfam, \$go);

	# calculate discoveryAproach
	$discoveryApproach = "";
	$discoveryApproach = &judgeDiscoveryApproach(\%isoformOrigin, $firstIsoformIdList, $secondIsoformIdList);

	# calculate jcecExpNum and jcExpNum
	if(exists(${$asExperimentNum{$asId}}{"jcecExpNum"})){
		$jcecExperimentNum = ${$asExperimentNum{$asId}}{"jcecExpNum"};
	}else{
		$jcecExperimentNum = 0;
	}

	if(exists(${$asExperimentNum{$asId}}{"jcExpNum"})){

		$jcExperimentNum = ${$asExperimentNum{$asId}}{"jcExpNum"};

	}else{

		$jcExperimentNum = 0;

	}
	
	# calculate cmbVarType of AS
	$variantPointTypeCmb = "None";
	$variantPointNum = 0;
	$variantPointTypeCmb = ${$asCmbVarType{$asId}}{"type"};
	if(exists(${$asCmbVarType{$asId}}{"type"})){
		$variantPointNum = ${$asCmbVarType{$asId}}{"num"};
	}else{		
		$variantPointNum = 0;
	}
	
	print WW "asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach" . "_____" . "\"" . ${$cmplAs{$asId}}{"asId"} . "\", \"" . ${$cmplAs{$asId}}{"asType"} . "\", \"" . ${$cmplAs{$asId}}{"species"} . "\", \"" . ${$cmplAs{$asId}}{"chr"} . "\", \"" . ${$cmplAs{$asId}}{"strand"} . "\", " . ${$cmplAs{$asId}}{"start"} . ", " . ${$cmplAs{$asId}}{"end"} . ", " . ${$cmplAs{$asId}}{"1stExonEnd"} . ", " . ${$cmplAs{$asId}}{"1stExonStart_0base"} . ", " . ${$cmplAs{$asId}}{"2ndExonEnd"} . ", " . ${$cmplAs{$asId}}{"2ndExonStart_0base"} . ", " . ${$cmplAs{$asId}}{"downstreamEE"} . ", " . ${$cmplAs{$asId}}{"downstreamES"} . ", " . ${$cmplAs{$asId}}{"exonEnd"} . ", " . ${$cmplAs{$asId}}{"exonStart_0base"} . ", " . ${$cmplAs{$asId}}{"flankingEE"} . ", " . ${$cmplAs{$asId}}{"flankingES"} . ", " . ${$cmplAs{$asId}}{"longExonEnd"} . ", " . ${$cmplAs{$asId}}{"longExonStart_0base"} . ", " . ${$cmplAs{$asId}}{"riExonEnd"} . ", " . ${$cmplAs{$asId}}{"riExonStart_0base"} . ", " . ${$cmplAs{$asId}}{"shortEE"} . ", " . ${$cmplAs{$asId}}{"shortES"} . ", " . ${$cmplAs{$asId}}{"upstreamEE"} . ", " . ${$cmplAs{$asId}}{"upstreamES"} . ", \"" . ${$cmplAs{$asId}}{"geneID"} . "\", \"" . ${$cmplAs{$asId}}{"geneSymbol"} . "\", \"" . ${$cmplAs{$asId}}{"firstAltExonSeries"} . "\", \"" . $firstIsoformIdList . "\", " . $firstIsoformNum . ", \"" . ${$cmplAs{$asId}}{"secondAltExonSeries"} . "\", \"" . $secondIsoformIdList . "\", " . $secondIsoformNum . ", \"" . $pfam . "\", \"" . $go . "\", " . $variantPointNum . ", \"" . $variantPointTypeCmb . "\", \"" . ${$cmplAs{$asId}}{"asOrthId"} . "\", " . ${$cmplAs{$asId}}{"conservedSpeciesNum"} . ", " . $jcecExperimentNum  . ", " . $jcExperimentNum  . ", \"" . $discoveryApproach . "\"\n";


}

close WW;

system("rm -rf " . $asTsvFile);
system("mv " . $tmpDir . "/tmp.as.with.complete.infor.tsv " . $asTsvFile);

# get isoformIds
sub getIsoformIdListPfamGo{
	my ($multiIsoformIdAndExonSeries, $firstAltExonSeries, $secondAltExonSeries, $firstIsoformIdListX, $firstIsoformNumX, $secondIsoformIdListX, $secondIsoformNumX,  $pfamTermsX, $goTermsX) = @_;
	my ($returnPfam, $returnGo);
	my (@isoformIdExonSeries, $isoformId, $exonSeries, $pfam, $go);
	my (@term, $term,  @tt);

	$$firstIsoformNumX = 0;
	$$secondIsoformNumX = 0;
	$$firstIsoformIdListX = "";
	$$secondIsoformIdListX = "";

	@isoformIdExonSeries = split(/\n/, $multiIsoformIdAndExonSeries);
	foreach my $isoformIdExonSeries(@isoformIdExonSeries){

		$isoformId = "";
		$exonSeries = "";

		($isoformId, $exonSeries, $pfam, $go) = split(/#/, $isoformIdExonSeries);

		if(index($exonSeries, $firstAltExonSeries)>=0){

			$$firstIsoformIdListX .= $isoformId . ",";
			$$firstIsoformNumX = $$firstIsoformNumX + 1;

			@term = ();
			@term = split(/,/, $pfam);
			foreach $term(@term){
				$returnPfam .= $term . "," if(index($returnPfam, $term)<0);
			}

			@term = ();
			@term = split(/,/, $go);
			foreach $term(@term){
				$returnGo .= $term . "," if(index($returnGo, $term)<0);
			}


		}elsif(index($exonSeries, $secondAltExonSeries)>=0){
			
			$$secondIsoformIdListX .= $isoformId . ",";
			$$secondIsoformNumX = $$secondIsoformNumX + 1;
	
			@term = ();
			@term = split(/,/, $pfam);
			foreach $term(@term){
				$returnPfam .= $term . "," if(index($returnPfam, $term)<0);
			}

			@term = ();
			@term = split(/,/, $go);
			foreach $term(@term){
				$returnGo .= $term . "," if(index($returnGo, $term)<0);
			}
		
		}
	}

	if($$firstIsoformIdListX ne ""){
		$$firstIsoformIdListX = substr($$firstIsoformIdListX, 0, length($$firstIsoformIdListX)-1);
	}

	if($$secondIsoformIdListX ne ""){
		$$secondIsoformIdListX = substr($$secondIsoformIdListX, 0, length($$secondIsoformIdListX)-1);
	}

	if($returnPfam ne ""){
		$$pfamTermsX = substr($returnPfam, 0, length($returnPfam) -1 );
	}
	if($returnGo ne ""){
		$$goTermsX = substr($returnGo, 0, length($returnGo) - 1);
	}
}

sub judgeDiscoveryApproach{

	my ($isoformOriginX, $firstIsoformIdListX, $secondIsoformIdListX) = @_;

	my %isoformOriginX = %$isoformOriginX;

	my ($firstIsoformOrigin, $secondIsoformOrigin, @tmp, $isoformIdX);

	@tmp = ();
	@tmp = split(/,/, $firstIsoformIdListX);
	foreach $isoformIdX(@tmp){
		$firstIsoformOrigin .= $isoformOriginX{$isoformIdX} . ",";
	}

	@tmp = ();
	@tmp = split(/,/, $secondIsoformIdListX);
	foreach $isoformIdX(@tmp){
		$secondIsoformOrigin.=$isoformOriginX{$isoformIdX};
	}

	if($firstIsoformIdListX eq "" or $secondIsoformIdListX eq ""){

		return "Novel";

	}elsif(index(uc($firstIsoformOrigin), "ENSEMBL") >=0 and index(uc($secondIsoformOrigin), "ENSEMBL") >=0){

		return "Ensembl";

	}elsif(index(uc($firstIsoformOrigin), "RNASEQ") <0 and index(uc($secondIsoformOrigin), "RNASEQ") <0){

		return "Ensembl_Refseq";

	}else{

		return "Ensembl_Refseq_RNAseq";

	}
}

sub getLocalTime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	return sprintf("%02d:%02d:%02d %02d-%02d-%2d", $hour, $min, $sec, $mon, $mday, $year);
}
