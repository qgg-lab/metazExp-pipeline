use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--exptIdToAssemblyList \\\n" .
                "--exptIdToReadLenList \\\n" .
                "--cDNAfa \\\n" .
		"--outputDir \n";
        exit;
}

my ($exptIdToAssemblyList, $exptIdToReadLenList, $cDNAfa,  $outputDir);

GetOptions(
        'exptIdToAssemblyList=s'=>\$exptIdToAssemblyList,
        'exptIdToReadLenList=s'=>\$exptIdToReadLenList,
	'outputDir=s'=>$outputDir,
        'cDNAfa=s'=>\$cDNAfa,
);

my ($line, $cmd, $exptId, $assembly);
open FF, "<$ARGV[1]";
while($line=<FF>){
	chomp($line);
	($exptId, $assembly) = split(/\t/, $line);
}
close FF;
