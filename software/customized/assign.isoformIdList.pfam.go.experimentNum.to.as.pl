#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--dbName asdb \\\n" . 
		"--dbUser lsas \\\n" . 
		"--dbPWD njaulsas2019 \\\n" . 
		"--species \"Bos taurus\"\\\n" . 
		"--speciesAbbr BTAU \\\n" .
		"--jcecExpNumFile jcec.experimentNum.txt \\\n".
		"--jcExpNumFile jc.experimentNum.txt \\\n".
		"--genomeFile genome.fa \\\n" .
		"--tmpDir ./tmpDir \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($jcecExpNumFile, $jcExpNumFile);
my ($species, $speciesAbbr, $genomeFile);
my ($tmpDir);
my ($dbh, $query);

GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'species=s'=>\$species,
	'speciesAbbr=s'=>\$speciesAbbr,
	'jcecExpNumFile=s'=>\$jcecExpNumFile,
	'jcExpNumFile=s'=>\$jcExpNumFile,
	'genomeFile=s'=>\$genomeFile,
	'tmpDir=s'=>\$tmpDir,
);


###### --- 1 --- ###########
# search all as and write into asTable.tsv file

my (%as, %geneToIsoform, %asToJcec, %asToJc);
my (@row, $line, $sql);
my ($asId, $geneId, $isoformId, $firstAltExonSeries, $secondAltExonSeries, $exonSeries);

$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;


$sql = "select asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base,	downstreamEE, downstreamES,exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE,shortES, upstreamEE, upstreamES, firstAltExonSeries, secondAltExonSeries, geneID, geneSymbol, firstIsoformIdList, firstIsoformNum, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum,discoveryApproach from asTable where species=\"" . $species . "\"";
$query=$dbh->prepare($sql);
$query->execute();

print "has executed query asTable.\n begin write:\n";
system("mkdir $tmpDir");
open WW, ">$tmpDir" . "/asTable.tsv";
while(@row=$query->fetchrow_array){
	print WW join("\t", @row) . "\n";
}
$query->finish();
$dbh->disconnect();
close WW;

print "Finish write as into file.\n";


########### --- 2 --- ###########################
# search all isoform and assign them to gene

#print "Begin search isoform\n";
my (%isoformOrigin);
$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $row;
$query = $dbh->prepare("select isoformId, geneId, exonSeries, isoformPfam, isoformGo, annoOrigin from isoformTable where species=\"" . $species . "\"");
$query->execute();
while($row=$query->fetchrow_hashref()){
	$geneToIsoform{$row->{"geneId"}} .= $row->{"isoformId"} . "#" . $row->{"exonSeries"} . "#" . $row->{"isoformPfam"} . "#" . $row->{"isoformGo"} . "\n";
	$isoformOrigin{$row->{"isoformId"}} = $row->{"annoOrigin"};
}
$query->finish();
$dbh->disconnect();
print "Finish load isoform\n";


############# --- 3 --- ###############################
# load jcecExperimentNum and jcExperimentNum into hash
open FF, "<$jcecExpNumFile";
#313 ECABA3SS0000000003
while($line = <FF>){
	if($line=~/^ *(\d+) +($speciesAbbr.*)\n/){
		${$as{$2}}{"jcecExpNum"} = $1;
	}
}
close FF;

open FF, "<$jcExpNumFile";
# 45 ECABA3SS0000000001
while($line = <FF>){
	if($line=~/^ *(\d+) +($speciesAbbr.*)\n/){
		${$as{$2}}{"jcExpNum"} = $1;
	}
}
close FF;

######## --- 4 : mask genome sequence with snp type #########
my (%genomeSeq, $line, @tt, $id);
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
print "Finish load genome sequence into hash.\n";

$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $sql = "select chr, start, end from asTable where species=\"$species\"";
$query = $dbh->prepare($sql);
$query->execute();
my $sblen;
while(my $row=$query->fetchrow_hashref()){
        $sblen = $row->{"end"} - $row->{"start"} + 1;
        substr($genomeSeq{$row->{"chr"}}, $row->{"start"}-1, $sblen) = " "x$sblen;
}
print "Finish mask as region with space.\n";
$query->finish();
$dbh->disconnect();

# flag snv type into genome sequence
$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
$sql = "select chr, pos, varType from variantTable where species=\"$species\"";
$query = $dbh->prepare($sql);
$query->execute();
while(my $row=$query->fetchrow_hashref()){
	if($row->{"varType"} eq "SNV"){
		substr($genomeSeq{$row->{"chr"}}, $row->{"pos"}-1, 1) = "S";
	}elsif($row->{"varType"} eq "deletion"){
		substr($genomeSeq{$row->{"chr"}}, $row->{"pos"}-1, 1) = "D";
	}elsif($row->{"varType"} eq "insertion"){
		substr($genomeSeq{$row->{"chr"}}, $row->{"pos"}-1, 1) = "I";
	}
}
$query->finish();
$dbh->disconnect();

print "Finish flag variant type in as regtion.\n";

my (%asVar, $asRegion);
$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
$sql = "select asId, chr, start, end from asTable where species=\"$species\"";
$query = $dbh->prepare($sql);
$query->execute();
my $sblen;
while(my $row=$query->fetchrow_hashref()){

        $sblen = $row->{"end"} - $row->{"start"} + 1;

	$asRegion = substr($genomeSeq{$row->{"chr"}}, $row->{"start"}-1, $sblen);
#	print "AS:" . $row->{"asId"} . " LEN:$sblen\n" . $asRegion . "\n";
#	<STDIN>;
	$asRegion=~s/ //g;
#	print "Removed space:$asRegion\n";
#	<STDIN>;

	if($asRegion=~/S/ and $asRegion=~/D/ and $asRegion=~/I/){
		${$asVar{$row->{"asId"}}}{"type"} = "SDI";
	}elsif($asRegion=~/S/ and $asRegion=~/D/){
		${$asVar{$row->{"asId"}}}{"type"} = "SD";
	}elsif($asRegion=~/S/ and $asRegion=~/I/){
		${$asVar{$row->{"asId"}}}{"type"} = "SI";
	}elsif($asRegion=~/D/ and $asRegion=~/I/){
		${$asVar{$row->{"asId"}}}{"type"} = "DI";
	}elsif($asRegion=~/S/){
		${$asVar{$row->{"asId"}}}{"type"} = "S";
	}elsif($asRegion=~/D/){
		${$asVar{$row->{"asId"}}}{"type"} = "D";
	}elsif($asRegion=~/I/){
		${$asVar{$row->{"asId"}}}{"type"} = "I";
	}else{
		${$asVar{$row->{"asId"}}}{"type"} = "None";
	}

#	print ${$asVar{$row->{"asId"}}}{"type"};
#	<STDIN>;

	if( ${$asVar{$row->{"asId"}}}{"type"} eq "None"){
		${$asVar{$row->{"asId"}}}{"num"} = 0;
	}else{
		${$asVar{$row->{"asId"}}}{"num"} = length($asRegion);
	}

}
$query->finish();
$dbh->disconnect();
print "Finish calculate as combined variant type.\n";

######### --- 5 ----
# delete as from asTable
$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
$sql ="delete from asTable where species=\"$species\"";
#print $sql;
#<STDIN>;
$query = $dbh->prepare($sql);
$query->execute();
$dbh->disconnect();
print "Finish delete alternative splicing events.\n";


############# --- 6 ---- ##########
# read each as from asTable.tsv
# obtain corresponding isoformIdList, pfam and go
# Insert into asTable
#print "Begin insert modified alternative splicing events into asTable.\n";

my ($asId, $asType, $species, $chr, $strand, $start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $firstAltExonSeries, $secondAltExonSeries, $geneID, $geneSymbol, $firstIsoformIdList, $firstIsoformNum, $secondIsoformIdList, $secondIsoformNum, $pfam, $go, $variantPointNum, $variantPointTypeCmb, $asOrthId, $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum, $discoveryApproach);

$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
open FF, "<$tmpDir" . "/asTable.tsv";

while(my $line=<FF>){
	($asId, $asType, $species, $chr, $strand, $start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $firstAltExonSeries, $secondAltExonSeries, $geneID, $geneSymbol, $firstIsoformIdList, $firstIsoformNum, $secondIsoformIdList, $secondIsoformNum, $pfam, $go, $variantPointNum, $variantPointTypeCmb, $asOrthId, $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum, $discoveryApproach) = split(/\t/, $line);

	&getIsoformIdListPfamGo($geneToIsoform{$geneID}, $firstAltExonSeries, $secondAltExonSeries, \$firstIsoformIdList, \$firstIsoformNum, \$secondIsoformIdList, \$secondIsoformNum, \$pfam, \$go);

	if(exists(${$as{$asId}}{"jcecExpNum"})){
		$jcecExperimentNum = ${$as{$asId}}{"jcecExpNum"};
	}else{
		$jcecExperimentNum = 0;
	}
	if(exists(${$as{$asId}}{"jcExpNum"})){
		$jcExperimentNum = ${$as{$asId}}{"jcExpNum"};
	}else{
		$jcExperimentNum = 0;
	}

	$discoveryApproach = &judgeDiscoveryApproach(\%isoformOrigin, $firstIsoformIdList, $secondIsoformIdList);


	if($asType eq "A5SS"){
		if($strand eq "+"){
			$start = $longExonStart_0base+1;
			$end = $flankingEE;
		}else{
			$start = $flankingES+1;
			$end = $longExonEnd
		}
	}elsif($asType eq "A3SS"){
		if($strand eq "+"){
			$start = $flankingES+1;
			$end = $longExonEnd;
		}else{
			$start = $longExonStart_0base + 1;
			$end = $flankingEE;
		}
	}
	if($asType eq "RI"){
		$start = $upstreamES + 1;
		$end = $downstreamEE;
	}
	if($asType eq "SE"){
		$start = $upstreamES + 1;
		$end = $downstreamEE;
	}
	if($asType eq "MXE"){
		$start = $upstreamES + 1;
		$end = $downstreamEE;
	}

	$variantPointTypeCmb = ${$asVar{$asId}}{"type"};

	if(exists(${$asVar{$asId}}{"type"})){
		$variantPointNum = ${$asVar{$asId}}{"num"};
	}else{		
		$variantPointNum = 0;
	}

	$sql = "insert into asTable (asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES,exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE,shortES, upstreamEE, upstreamES, firstAltExonSeries, secondAltExonSeries, geneID, geneSymbol, firstIsoformIdList, firstIsoformNum, secondIsoformIdList, secondIsoformNum, pfam, go,variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach) values (\"$asId\", \"$asType\", \"$species\", \"$chr\", \"$strand\", $start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, \"$firstAltExonSeries\", \"$secondAltExonSeries\", \"$geneID\", \"$geneSymbol\", \"$firstIsoformIdList\", $firstIsoformNum, \"$secondIsoformIdList\", $secondIsoformNum, \"$pfam\", \"$go\", $variantPointNum, \"$variantPointTypeCmb\", \"$asOrthId\", $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum, \"$discoveryApproach\")";
	#<STDIN>;
	$query = $dbh->prepare($sql);
	$query->execute();
}
close FF;
$dbh->disconnect();

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
