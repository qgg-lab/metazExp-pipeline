#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--interproscan interproscan.sh \\\n" .
		"--outputDir . \\\n" .
		"--inputPepSeqFile ./final.complete.trspt.pep.fa \\\n" . 
		"--taxonId 9913 \\\n";
	exit;
}
my ($outputDir);
my ($interproscan, $inputPepSeqFile, $inputCdnaSeqFile, $taxonId, $jobName);

GetOptions(
	'interproscan=s'=>\$interproscan,
	'outputDir=s'=>\$outputDir,
	'inputPepSeqFile=s'=>\$inputPepSeqFile,
	'taxonId=s'=>\$taxonId,
	'jobName=s'=>\$jobName,
);

system("mkdir -p " . $outputDir);


########################################
# split cDNAseqNotInPepFile
#
my ($cmd, $fileNameWithPathList, @fileNameWithPath, @tt, $fileName, $fileNameWithPath);

###################################
# split seqInPepFile
#
$cmd = "split -l 4000 " . $inputPepSeqFile . " " . $outputDir . "/pep_";
system($cmd);
$cmd = "ls -1 " . $outputDir . "/pep_*";
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
	print WW "#SBATCH --job-name=$jobName" . "_" . $taxonId . "_" . $fileName ."\n";
	print WW "#SBATCH --nodes=1\n";
	print WW "#SBATCH --ntasks-per-node=1\n";
	print WW "#SBATCH --cpus-per-task=14\n";
	print WW "#SBATCH --mem=64G\n";
	print WW "#SBATCH --time=08:00:00\n";
	print WW "date > " . $outputDir . "/" . $fileName . ".time.txt\n";
	print WW $interproscan . " -b " . $outputDir . "/" . $fileName . " -cpu 14 -i " . $fileNameWithPath . " -f TSV -goterms -t p -T " . $outputDir . "/tmpDir_" . $fileName . " 1> " . $outputDir . "/log.o." . $fileName . ".sb 2> " . $outputDir . "/log.e." . $fileName . ".sb\n";
	print WW "date >> " . $outputDir . "/" . $fileName . ".time.txt\n";
	close WW;

	print SUBMIT "sbatch " . $outputDir . "/" . $fileName . ".sb\n";
}

close SUBMIT;

########################
# chmod ##########

system("chmod +x " . $outputDir . "/submit.cmd.sh");

# system($outputDir . "/submit.cmd.sh");
