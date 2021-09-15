#!/usr/bin/perl
use strict;
print &calPhredScore($ARGV[0], 100000);
sub calPhredScore{
	my ($fastqFile, $readNum) = @_;
	my ($ltZeroBaseCount, $read, $score, $totalBaseCount);
	$ltZeroBaseCount = 0;
	open FF, "<$fastqFile";
	for(my $i=0; $i<$readNum*4; $i++){
		$read = <FF>;
		next if($i-int($i/4)*4 != 3);
		for(my $j=0; $j<=length($read)-2; $j++){
			$score = ord(substr($read, $j, 1))-64;
			$totalBaseCount++;
			$ltZeroBaseCount++ if($score < 0);
		}
	}
	if($ltZeroBaseCount > 0){
		return 33;
	}else{
		return 64;
	}	
}
