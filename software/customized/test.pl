
#print &calPhredScore("/mnt/scratch/liujind1/fastqDir/SRR2226585_1.fastq", 10000) . "\n";
my $st;
print &obtainLibraryType("/mnt/scratch/liujind1/fastqOutputDir/SRX1302200/SRR2551770.map.cDNA.sam", 1000000, "PE", 85, \$st);

sub obtainLibraryType{
	my ($map_sam, $num, $layout, $supportPercentage) = @_;
	my (%sam, $line, @tmp, @binarr, $binstr, $i, $line);
	my ($total_statistic_num, $RF_num, $FR_num, $R_num, $F_num);
	my ($tempStr);
	print $map_sam . "\n";
	<STDIN>;
	goto SINGLE if($layout eq "SINGLE");

	open FSAM, "<$map_sam";
	while($line=<FSAm>){
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
	close FSAM;

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

	$total_statistic_num = $#tmp + 1;
	if($FR_num/$total_statistic_num> $supportPercentage/100){
		return "FR";
	}elsif($RF_num/$total_statistic_num>$supportPercentage/100){
		return "RF";
	}else{
		return "UN";
	}

SINGLE:
	open FSAM, "<$map_sam";
	while($line=<FSAM>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);

		#It is header
		next if($#tmp<10);
		print $line . "\n";
		<STDIN>;
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
	close FSAM;


	$total_statistic_num = $R_num + $F_num;
	if($F_num/$total_statistic_num>$supportPercentage/100){
		return "F";
	}elsif($R_num/$total_statistic_num>$supportPercentage/100){
		return "R";
	}else{
		return "U";
	}
}


sub calPhredScore{
        my ($fastqFile, $readNum) = @_;
        my ($ltZeroBaseCount, $read, $score, $totalBaseCount);
        $ltZeroBaseCount = 0;
        open FPHRED, "<$fastqFile";
        for(my $i=0; $i<$readNum*4; $i++){
                $read = <FPHRED>;
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
        close FPHRED;
}
