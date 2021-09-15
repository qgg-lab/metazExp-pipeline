#!/usr/bin/perl
use strict;
my ($dataFile);

$dataFile = $ARGV[0];
my (%expt, $line, @field, $exptHref, $exptId);
$exptHref = \%expt;

open FF, "<$dataFile";
# Taxon   ExptId  AsNumName       Num
# 108875  SRX2764785      gtfA3SS 4234
# 108875  SRX2764785      novelA3SS       0
# 108875  SRX2764785      gtfA5SS 2051
# 108875  SRX2764785      novelA5SS       0
# 108875  SRX2764785      gtfMXE  50
# 108875  SRX2764785      novelMXE        154
# 108875  SRX2764785      gtfRI   5169
# 108875  SRX2764785      novelRI 391
# 108875  SRX2764785      gtfSE   1918
# 108875  SRX2764785      novelSE 3778
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	$exptId = $field[1];
#	print $exptId;
#	<STDIN>;
	if($field[2]=~/A3SS/){
		$exptHref->{$exptId}->{"A3SS"}+=$field[3];
	}elsif($field[2]=~/A5SS/){
		$exptHref->{$exptId}->{"A5SS"}+=$field[3];
	}elsif($field[2]=~/SE/){
		$exptHref->{$exptId}->{"SE"}+=$field[3];
	}elsif($field[2]=~/RI/){
		$exptHref->{$exptId}->{"RI"}+=$field[3];
	}elsif($field[2]=~/MXE/){
		$exptHref->{$exptId}->{"MXE"}+=$field[3];
	}
}
close FF;

my $total;
my @exptId = keys(%expt);
open WW, ">" . $ARGV[1];
print WW join("\t", "Type", "Percentage") . "\n";

foreach $exptId(@exptId){

	$total = $exptHref->{$exptId}->{"A3SS"} + $exptHref->{$exptId}->{"A5SS"} + $exptHref->{$exptId}->{"SE"} + $exptHref->{$exptId}->{"RI"} + $exptHref->{$exptId}->{"MXE"};
	next if($total==0);
	print WW join("\t", "A3SS", $exptHref->{$exptId}->{"A3SS"}/$total*100) . "\n";
	print WW join("\t", "A5SS", $exptHref->{$exptId}->{"A5SS"}/$total*100) . "\n";
	print WW join("\t", "SE", $exptHref->{$exptId}->{"SE"}/$total*100) . "\n";
	print WW join("\t", "MXE", $exptHref->{$exptId}->{"MXE"}/$total*100) . "\n";
	print WW join("\t", "RI", $exptHref->{$exptId}->{"RI"}/$total*100) . "\n";
}
close WW;
