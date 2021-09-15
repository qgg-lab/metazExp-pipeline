#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origGtf \\\n" .
                "--origPepFasta \\\n" .
                "--simpleTrsptIdGtf  \\\n" .
		"--speciesAbbr \\\n" .
                "--simpleTrsptIdPepFasta \\\n" .
                "--mappingBetweenNewAndOldTrsptIdTsv  \n";
	exit;
}

my ($speciesAbbr, $origGtf, $origPepFasta, $simpleTrsptIdGtf, $simpleTrsptIdPepFasta, $mappingBetweenNewAndOldTrsptIdTsv, $speciesAbbr);

GetOptions(
        'speciesAbbr=s'=>\$speciesAbbr,
        'origGtf=s'=>\$origGtf,
        'origPepFasta=s'=>\$origPepFasta,
        'simpleTrsptIdGtf=s'=>\$simpleTrsptIdGtf,
        'simpleTrsptIdPepFasta=s'=>\$simpleTrsptIdPepFasta,
	'mappingBetweenNewAndOldTrsptIdTsv=s'=>\$mappingBetweenNewAndOldTrsptIdTsv,
);

# 读取原始pep序列编号，构建hash
# >AT2G01420___AT2G01420.2XXXATHAA5SS0000002139#AT2G01420.1
my (%pep, $pepHref, %geneIdToTrsptId, $geneIdToTrsptIdHref);
my (@seqId, $seqId, @trsptId, $trsptIdNum, $geneId, $trsptIdList, $trsptId);
$pepHref=\%pep;
$geneIdToTrsptIdHref = \%geneIdToTrsptId;
open FF, "<$origPepFasta";
while(my $line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$seqId = $1;
		@seqId = ();
		@seqId = split(/ /, $seqId);
		$seqId = $seqId[0];
		($geneId, $trsptIdList) = split(/___/, $seqId);
		$geneIdToTrsptIdHref->{$geneId}->{$trsptIdList} = 1;
		@trsptId = ();
		@trsptId = keys(%{$geneIdToTrsptIdHref->{$geneId}});
		$trsptIdNum = $#trsptId + 1;
#		print $seqId;
#		<STDIN>;
		$pepHref->{$seqId}->{"simpleTrsptId"} = $speciesAbbr . "_" . $geneId . "_pep" . sprintf("%03d", $trsptIdNum);
#		print "simpleTrsptId:" . $pepHref->{$seqId}->{"simpleTrsptId"} . "\n";
#		<STDIN>;
	}else{
		$pepHref->{$seqId}->{"pepSeq"} .= $line;
	}
}
close FF;

# 将新旧trscptId映射关系输出
open MAP, ">$mappingBetweenNewAndOldTrsptIdTsv";
open FASTA, ">$simpleTrsptIdPepFasta";
@trsptId = ();
@trsptId = keys(%pep);
#print $#trsptId;
#<STDIN>;
foreach $trsptId(@trsptId){
	print FASTA ">" . $pepHref->{$trsptId}->{"simpleTrsptId"} . "\n";
	print FASTA $pepHref->{$trsptId}->{"pepSeq"} . "\n";
	print MAP join("\t", $trsptId, $pepHref->{$trsptId}->{"simpleTrsptId"}) . "\n";
}
close MAP;
close FASTA;


# 读取gtf文件，将其中的transcript_id指定的trsptId改变为simpleTrsptId
my (@field, @attr);
open FF, "<$origGtf";
open WW, ">$simpleTrsptIdGtf";
while(my $line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	@attr = ();
	@attr = split(/; /, $field[8]);
	if($attr[1]=~/transcript_id "(.*?)"/){
		$attr[1] = "transcript_id \"" . $pepHref->{$1}->{"simpleTrsptId"} . "\"";
	}
	$field[8] = join("; ", @attr);
	print WW join("\t", @field) . "\n";
}
close FF;
close WW;
