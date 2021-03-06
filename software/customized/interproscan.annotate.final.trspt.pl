#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--interproscan interproscan.sh \\\n" .
		"--outputDir . \\\n" .
		"--inputCdnaSeqFile ./final.complete.trspt.cDNA.fa \\\n".
		"--inputPepSeqFile ./final.complete.trspt.pep.fa \\\n" . 
		"--taxonId 9913 \\\n";
	exit;
}
my ($outputDir);
my ($interproscan, $inputPepSeqFile, $inputCdnaSeqFile, $taxonId);

GetOptions(
	'interproscan=s'=>\$interproscan,
	'outputDir=s'=>\$outputDir,
	'inputCdnaSeqFile=s'=>\$inputCdnaSeqFile,
	'inputPepSeqFile=s'=>\$inputPepSeqFile,
	'taxonId=s'=>\$taxonId,
);

system("mkdir -p " . $outputDir);

my (%cDNAseq, %pepSeq, $id);
open FF, "<$inputCdnaSeqFile";
while(my $line =<FF>){
        chomp($line);
        if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*)/){
                $id = $1;
        }else{
                $cDNAseq{$id}.=$line;
        }
}
close FF;

open FF, "<$inputPepSeqFile";
while(my $line =<FF>){
        chomp($line);
        if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*) protein_id:(.*)/){
                $id = $1;
        }else{
                $pepSeq{$id}.=$line;
        }
}
close FF;

my @cDNAid=keys(%cDNAseq);
my $cDNAseqNotInPepFile = "$outputDir" . "/cDNAseqNotInPepFile.fa";
open WW, ">$cDNAseqNotInPepFile";
foreach $id(@cDNAid){
	if(not(exists($pepSeq{$id}))){
		print WW ">$id\n";
		print WW $cDNAseq{$id} . "\n";
	}
}
close WW;

########################################
# split cDNAseqNotInPepFile
#
my ($cmd, $fileNameWithPathList, @fileNameWithPath, @tt, $fileName, $fileNameWithPath);

$cmd = "split -l 4000 " . $outputDir . "/cDNAseqNotInPepFile.fa " . $outputDir . "/cDNA_";
system($cmd);
$cmd = "ls -1 " . $outputDir . "/cDNA_*";
$fileNameWithPathList = `$cmd`;
@fileNameWithPath = ();
@fileNameWithPath = split(/\n/, $fileNameWithPathList);

open SUBMIT, ">" . $outputDir . "/submit.cmd.sh";
foreach $fileNameWithPath(@fileNameWithPath){
	@tt = ();
	@tt = split(/\//, $fileNameWithPath);
	$fileName = $tt[$#tt];

	open WW, ">" . $outputDir . "/" . $fileName . ".sb";
	print WW "#!/bin/bash\n";
	print WW "#SBATCH --job-name=006_" . $taxonId . "_" . $fileName . "\n";
	print WW "#SBATCH --nodes=1\n";
	print WW "#SBATCH --ntasks-per-node=1\n";
	print WW "#SBATCH --cpus-per-task=18\n";
	print WW "#SBATCH --mem=64G\n";
	print WW "#SBATCH --time=03:50:00\n";
	print WW "date > " . $outputDir . "/" . $fileName . ".time.txt\n";
	print WW $interproscan . " -b " . $outputDir . "/" . $fileName . " -cpu 18 -i " . $fileNameWithPath . " -f TSV -goterms -t n -T " . $outputDir . "/tmpDir_" . $fileName . " 1> " . $outputDir . "/log.o." . $fileName . ".sb 2> " . $outputDir . "/log.e." . $fileName . ".sb\n";
	print WW "date >> " . $outputDir . "/" . $fileName . ".time.txt\n";
	close WW;

	print SUBMIT "sbatch " . $outputDir . "/" . $fileName . ".sb\n";
}



###################################
# split cDNAseqInPepFile
#
$cmd = "split -l 4000 " . $inputPepSeqFile . " " . $outputDir . "/pep_";
system($cmd);
$cmd = "ls -1 " . $outputDir . "/pep_*";
$fileNameWithPathList = `$cmd`;
@fileNameWithPath = ();
@fileNameWithPath = split(/\n/, $fileNameWithPathList);

foreach $fileNameWithPath(@fileNameWithPath){
	@tt = ();
	@tt = split(/\//, $fileNameWithPath);
	$fileName = $tt[$#tt];

	open WW, ">" . $outputDir . "/" . $fileName . ".sb";
	print WW "#!/bin/bash\n";
	print WW "#SBATCH --job-name=006_" . $taxonId . "_" . $fileName ."\n";
	print WW "#SBATCH --nodes=1\n";
	print WW "#SBATCH --ntasks-per-node=1\n";
	print WW "#SBATCH --cpus-per-task=18\n";
	print WW "#SBATCH --mem=64G\n";
	print WW "#SBATCH --time=03:50:00\n";
	print WW "date > " . $outputDir . "/" . $fileName . ".time.txt\n";
	print WW $interproscan . " -b " . $outputDir . "/" . $fileName . " -cpu 18 -i " . $fileNameWithPath . " -f TSV -goterms -t p -T " . $outputDir . "/tmpDir_" . $fileName . " 1> " . $outputDir . "/log.o." . $fileName . ".sb 2> " . $outputDir . "/log.e." . $fileName . ".sb\n";
	print WW "date >> " . $outputDir . "/" . $fileName . ".time.txt\n";
	close WW;

	print SUBMIT "sbatch " . $outputDir . "/" . $fileName . ".sb\n";
}

close SUBMIT;

########################
# chmod ##########

system("chmod +x " . $outputDir . "/submit.cmd.sh");

# system($outputDir . "/submit.cmd.sh");
