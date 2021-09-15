#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD lsasNjau2019 \\\n" .
		"--speciesList \"Gallus gallus,Equus caballus,Sus scrofa,Bos taurus,Ovis aries\" \\\n" .
		"--tissueList \"Gallus gallus:liver,embryos,spleen,adipose,bursa,muscle,pituitary gland,lung,macrophages,heart,mesenchymal stem cells,utricle,retina,thymus,kidney;Equus caballus:PBMCs,blood,muscle,oviduct,gluteus medius,dermal cells,spinal cord,cerebellum,brainstem,inner cell mass,trophectoderm,liver,trophoblasts,skin,laminae;Sus scrofa:blood,muscle,adipose,endometrium,fetal thymus,heart,liver,ganglion,ovary,spleen,lung,hippocampus,knee synovial membrane,fibroblasts,uterus;Bos taurus:muscle,liver,blood,endometrium,macrophages,lymph nodes,mammary gland,pituitary gland,granulosa cells,embryos,adipose,rumen,fibroblasts,hypothalamus,kidney;Ovis aries:ovary,muscle,macrophages,adipose,lymph nodes,milk cells,abomasum,colon,conceptus,myocardium,liver,ileum,rumen,ruminant reticulum,mammary gland\" \\\n" .
		"--Rscript /usr/local/bin/Rscript \\\n" .
		"--boxPlotRscript psi.Boxplot.R \\\n" .
		"--outputDir /var/www/html/lsas/boxplot \n\n"; 
	exit;
}

my ($dbName, $dbUser, $dbPWD, $speciesList, $tissueList, $Rscript, $boxPlotRscript, $outputDir);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'speciesList=s'=>\$speciesList,
        'tissueList=s'=>\$tissueList,
        'Rscript=s'=>\$Rscript,
        'boxPlotRscript=s'=>\$boxPlotRscript,
        'outputDir=s'=>\$outputDir,
);

my ($dbh, $sqlAs, $queryAs, $rowAs);
my ($sqlPsi, $queryPsi, $rowPsi);
my (@tissue, $tissue, $condition, $asId);

@tissue = split(/,/, $tissueList);
foreach $tissue(@tissue){
	$condition .= "tissue =\"" . $tissue . "\" or ";
}
$condition = substr($condition, 0, length($condition) - 3);


#select asId from asTable where FIND_IN_SET(species, "Gallus gallus,Equus caballus,Sus scrofa,Bos taurus,Ovis aries");
$dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);
$dbh->{mysql_auto_reconnect} = 1;
$sqlAs = "select asId from asTable where FIND_IN_SET(species, \"$speciesList\")";
$queryAs = $dbh->prepare($sqlAs);
$queryAs->execute();
while($rowAs=$queryAs->fetchrow_hashref()){
	$asId = $rowAs->{"asId"};
	
	$sqlPsi = "select tissue, JCECpsi from psiTable where asId=\"" . $asId . "\" and (" . $condition . ")";
	$queryPsi = $dbh->prepare($sqlPsi);
	$queryPsi->execute();

	open WW, ">$outputDir/psi.tsv";
	print WW join("\t", "Tissue", "Psi") . "\n";
	while($rowPsi=$queryPsi->fetchrow_hashref()){
		print WW join("\t", $rowPsi->{"tissue"}, $rowPsi->{"JCECpsi"}) . "\n";
	}
	close WW;

}
$queryAs->finish();
$dbh->disconnect();
