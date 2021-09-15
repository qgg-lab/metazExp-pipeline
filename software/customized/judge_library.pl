#!/usr/bin/perl
use strict;
use strict;
use Getopt::Long;
use FindBin;
my $strandSupportPercentage;

if($#ARGV<0){
	print $0 . "map.sam 10000 PE/SE\n";
	exit;
}

#map.sam 10000 PE/SE
my $library_type = &library_type($ARGV[0], $ARGV[1], $ARGV[2], \$strandSupportPercentage);

print $strandSupportPercentage . "\n";

sub library_type{
	my ($map_sam, $num, $pe_se, $supportPercentage) = @_;
	my (%sam, $line, @tmp, @binarr, $binstr, $i, $line);
	my ($total_statistic_num, $RF_num, $FR_num, $R_num, $F_num);
	my ($tempStr);

	goto SE if($pe_se eq "SE");

	open FF, "<$map_sam";
	while($line=<FF>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);

		#It is header
		next if($#tmp<10);

		#parse alignment information
		@binarr = (0);
		$binstr=sprintf("%b", $tmp[1]);
        	for($i=length($binstr); $i>=1; $i--){
			$binarr[$#binarr+1]=substr($binstr, $i-1, 1);
	        }

		#not proper pair end reads
		next if(not($tmp[1]==99 or $tmp[1]==147 or $tmp[1]==83 or $tmp[1]==163));

		#two reads are not in a same refseq
		next if($tmp[6] ne "="); 

		#hit multiple sites
		next if(not($tmp[$#tmp]=~/NH:i:1/));
	
		#register read1 strand and position
		if($binarr[7]==1){
			${${$sam{$tmp[0]}}{"read1"}}{"position"}=$tmp[3];
			if($binarr[5]==1){
				${${$sam{$tmp[0]}}{"read1"}}{"strand"}="-";
			}else{
				${${$sam{$tmp[0]}}{"read1"}}{"strand"}="+";
			}
		}
	
		#register read2 strand and position
		if($binarr[8] == 1){
			${${$sam{$tmp[0]}}{"read2"}}{"position"}=$tmp[3];
			if($binarr[5]==1){
				${${$sam{$tmp[0]}}{"read2"}}{"strand"}="-";
			}else{
				${${$sam{$tmp[0]}}{"read2"}}{"strand"}="+";
			}
		}
		$num--;
		last if($num<=0);
	}
	close FF;

	#calculate numbers of RF and FR pair-ends

	@tmp = ();
	@tmp = keys%sam;

	for($i=0; $i<=$#tmp; $i++){
		if(${${$sam{$tmp[$i]}}{"read1"}}{"strand"} eq "+" and ${${$sam{$tmp[$i]}}{"read2"}}{"strand"} eq "-" and ${${$sam{$tmp[$i]}}{"read1"}}{"position"} <= ${${$sam{$tmp[$i]}}{"read2"}}{"position"}){
		$FR_num++;
		}elsif(${${$sam{$tmp[$i]}}{"read1"}}{"strand"} eq "-" and ${${$sam{$tmp[$i]}}{"read2"}}{"strand"} eq "+" and ${${$sam{$tmp[$i]}}{"read1"}}{"position"} >= ${${$sam{$tmp[$i]}}{"read2"}}{"position"}){
		$RF_num++;
		}
	}

	$total_statistic_num = $#tmp+1;
	#print "Total of $total_statistic_num PEs were statistic\n";
	#print "FR account for " . $FR_num/$total_statistic_num . "\n";
	#print "RF account for " . $RF_num/$total_statistic_num . "\n";
	$tempStr = "";
	$$supportPercentage="";
	if($FR_num/$total_statistic_num>0.85){
		$tempStr = $FR_num/$total_statistic_num;
		$$supportPercentage = "RF:$tempStr";
		return "FR";
	}elsif($RF_num/$total_statistic_num>0.85){
		$tempStr = $RF_num/$total_statistic_num;
		$$supportPercentage = "RF:$tempStr";;
		return "RF";
	}else{
		$tempStr = $FR_num/$total_statistic_num;
		$$supportPercentage="FR:($tempStr),";
		$tempStr = $RF_num/$total_statistic_num;
		$$supportPercentage.="RF:($tempStr)";
		return "UN";
	}

SE:
	open FF, "<$map_sam";
	while($line=<FF>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);

		#It is header
		next if($#tmp<10);

		#parse alignment information
		@binarr = (0);
		$binstr=sprintf("%b", $tmp[1]);
        	for($i=length($binstr); $i>=1; $i--){
			$binarr[$#binarr+1]=substr($binstr, $i-1, 1);
	        }

		#hit multiple sites
		next if(not($tmp[$#tmp]=~/NH:i:1/));
	
		#register alignment strand
		if($binarr[5]==1){
			$R_num++;
		}else{
			$F_num++;
		}

		$num--;
		last if($num<=0);
	}
	close FF;

	#calculate numbers of R and F

	$total_statistic_num = $R_num + $F_num;
	#print "Total of $total_statistic_num SE reads were statistic\n";
	#print "F account for " . $F_num/$total_statistic_num . "\n";
	#print "R account for " . $R_num/$total_statistic_num . "\n";
	$$supportPercentage="";
	$tempStr = "";
	if($F_num/$total_statistic_num>0.85){
		$tempStr = $F_num/$total_statistic_num;
		$$supportPercentage="F($tempStr)";
		return "F";
	}elsif($R_num/$total_statistic_num>0.85){
		$tempStr = $R_num/$total_statistic_num;
		$$supportPercentage = "R($tempStr)";
		return "R";
	}else{
		$tempStr = $F_num/$total_statistic_num;
		$$supportPercentage.="F($tempStr),";
		$tempStr = $R_num/$total_statistic_num;
		$$supportPercentage.="R($tempStr)";
		return "U";
	}
}

