#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--asTsvFile 89462.as.tsv\\\n" .
		"--tmpDir tmpDir \n\n";
        exit;
}

my ($dbName, $dbUser, $dbPWD, $asTsvFile, $tmpDir);
my (%as, $firstExonSeries, $secondExonSeries);
my ($line, @fieldName, @values, $fieldString, $valueString, @tmp, $asId);
my ($fieldName, $value, $i, $start, $end);
GetOptions(
	'asTsvFile=s'=>\$asTsvFile,
	'tmpDir=s'=>\$tmpDir,
);
system("mkdir -p " . $tmpDir);
open WW, ">$tmpDir" . "/tmp.as.tsv";
open FF, "<$asTsvFile";
while($line=<FF>){
	chomp($line);
#	print $line;
#	<STDIN>;
	@tmp = ();
	@tmp = split(/_____/, $line);
	$fieldString = $tmp[0];
#	print $fieldString;
#	<STDIN>;
	@fieldName = ();
	@fieldName = split(/, /, $fieldString);

	$valueString = $tmp[1];
#	print $valueString;
#	<STDIN>;
	@values = ();
	@values = split(/, /, $valueString);

	for($i=0; $i<=$#values; $i++){
		if($values[$i]=~/"(.*)"/){
			$values[$i] = $1;
		}
		$as{$fieldName[$i]}=$values[$i];
	}

#	print $#fieldName;
#	<STDIN>;
#	print $#values;
#	<STDIN>;

	$start = 0;
	$end = 0;
	#print join("\n", $as{"asType"}, $as{"strand"}, $as{"1stExonEnd"}, $as{"1stExonStart_0base"}, $as{"2ndExonEnd"}, $as{"2ndExonStart_0base"}, $as{"downstreamEE"}, $as{"downstreamES"}, $as{"exonEnd"}, $as{"exonStart_0base"}, $as{"flankingEE"}, $as{"flankingES"}, $as{"longExonEnd"}, $as{"longExonStart_0base"}, $as{"riExonEnd"}, $as{"riExonStart_0base"}, $as{"shortEE"}, $as{"shortES"}, $as{"upstreamEE"}, $as{"upstreamES"}) . "\n";
	#<STDIN>;

	&getAsStartAndEnd($as{"asType"}, $as{"strand"}, $as{"1stExonEnd"}, $as{"1stExonStart_0base"}, $as{"2ndExonEnd"}, $as{"2ndExonStart_0base"}, $as{"downstreamEE"}, $as{"downstreamES"}, $as{"exonEnd"}, $as{"exonStart_0base"}, $as{"flankingEE"}, $as{"flankingES"}, $as{"longExonEnd"}, $as{"longExonStart_0base"}, $as{"riExonEnd"}, $as{"riExonStart_0base"}, $as{"shortEE"}, $as{"shortES"}, $as{"upstreamEE"}, $as{"upstreamES"}, \$start, \$end);

	$as{"start"} = $start;
	$as{"end"} = $end;

	#print join("\t", $as{"start"}, $as{"end"});
	#<STDIN>;
	$as{"firstAltExonSeries"} = &generateFirstAltAsExonSeries($as{"asType"}, $as{"strand"}, $as{"1stExonEnd"}, $as{"1stExonStart_0base"}, $as{"2ndExonEnd"}, $as{"2ndExonStart_0base"}, $as{"downstreamEE"}, $as{"downstreamES"}, $as{"exonEnd"}, $as{"exonStart_0base"}, $as{"flankingEE"}, $as{"flankingES"}, $as{"longExonEnd"}, $as{"longExonStart_0base"}, $as{"riExonEnd"}, $as{"riExonStart_0base"}, $as{"shortEE"}, $as{"shortES"}, $as{"upstreamEE"}, $as{"upstreamES"});

	$as{"secondAltExonSeries"} = &generateSecondAltAsExonSeries($as{"asType"}, $as{"strand"}, $as{"1stExonEnd"}, $as{"1stExonStart_0base"}, $as{"2ndExonEnd"}, $as{"2ndExonStart_0base"}, $as{"downstreamEE"}, $as{"downstreamES"}, $as{"exonEnd"}, $as{"exonStart_0base"}, $as{"flankingEE"}, $as{"flankingES"}, $as{"longExonEnd"}, $as{"longExonStart_0base"}, $as{"riExonEnd"}, $as{"riExonStart_0base"}, $as{"shortEE"}, $as{"shortES"}, $as{"upstreamEE"}, $as{"upstreamES"});

	#print join("\n", $as{"firstAltExonSeries"}, $as{"secondAltExonSeries"}) . "\n";
	#<STDIN>;
	print WW $fieldString . "_____";

	print WW join(", ", "\"" .$as{"asId"} . "\"" , "\"" . $as{"asType"} . "\"", "\"" . $as{"species"} . "\"", "\"" . $as{"chr"} . "\"", "\"" . $as{"strand"} . "\"", $as{"start"}, $as{"end"}, $as{"1stExonEnd"}, $as{"1stExonStart_0base"}, $as{"2ndExonEnd"}, $as{"2ndExonStart_0base"}, $as{"downstreamEE"}, $as{"downstreamES"}, $as{"exonEnd"}, $as{"exonStart_0base"}, $as{"flankingEE"}, $as{"flankingES"}, $as{"longExonEnd"}, $as{"longExonStart_0base"}, $as{"riExonEnd"}, $as{"riExonStart_0base"}, $as{"shortEE"}, $as{"shortES"}, $as{"upstreamEE"}, $as{"upstreamES"}, "\"" . $as{"geneID"} . "\"", "\"" . $as{"geneSymbol"} . "\"", "\"" . $as{"firstAltExonSeries"} . "\"", "\"" . $as{"firstIsoformIdList"} . "\"", $as{"firstIsoformNum"}, "\"" . $as{"secondAltExonSeries"} . "\"", "\"" . $as{"secondIsoformIdList"} . "\"", $as{"secondIsoformNum"}, "\"" . $as{"pfam"} . "\"", "\"" . $as{"go"} . "\"", $as{"variantPointNum"}, "\"" . $as{"variantPointTypeCmb"} . "\"", "\"" . $as{"asOrthId"} . "\"", $as{"conservedSpeciesNum"}, $as{"jcecExperimentNum"}, $as{"jcExperimentNum"}, "\"" . $as{"discoveryApproach"} . "\"") . "\n";
}
close FF;
close WW;

system("rm -rf $asTsvFile");
system("mv " . $tmpDir . "/tmp.as.tsv " . $asTsvFile);

sub generateFirstAltAsExonSeries{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $longExonStart_0base . ".." . $longExonEnd . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;
		return $flankingES . ".." . $flankingEE . "," . $longExonStart_0base . ".." . $longExonEnd;	
	}elsif($asType eq "RI"){
		$riExonStart_0base = $riExonStart_0base + 1;
		$upstreamES = $upstreamES + 1;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE;
		}
	}elsif($asType eq "SE"){
		$exonStart_0base = $exonStart_0base + 1;
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $exonStart_0base . ".." . $exonEnd;
		}
	}elsif($asType eq "MXE"){
		$firstExonStart_0base = $firstExonStart_0base + 1;
		$secondExonStart_0base = $secondExonStart_0base + 1;
		$upstreamES = $upstreamES + 1 ;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $firstExonStart_0base . ".." . $firstExonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}
	
}

sub generateSecondAltAsExonSeries{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES)=@_;
	if($asType eq "A5SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;		
		return $shortES . ".." . $shortEE . "," . $flankingES . ".." . $flankingEE;		
	}elsif($asType eq "A3SS"){
		$longExonStart_0base = $longExonStart_0base + 1;
		$shortES = $shortES + 1;
		$flankingES = $flankingES + 1;	
		return $flankingES . ".." . $flankingEE . "," . $shortES . ".." . $shortEE ;	
	}elsif($asType eq "RI"){
		$riExonStart_0base = $riExonStart_0base + 1;
		return $riExonStart_0base . ".." . $riExonEnd;
	}elsif($asType eq "SE"){
		$exonStart_0base = $exonStart_0base + 1;
		$upstreamES = $upstreamES +1;
		$downstreamES = $downstreamES +1;
		if($strand eq "+"){
			return $upstreamES . ".." .  $upstreamEE . "," . $exonStart_0base . ".." . $exonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $upstreamES . ".." .  $upstreamEE . "," . $exonStart_0base . ".." . $exonEnd;
		}
	}elsif($asType eq "MXE"){
		$firstExonStart_0base = $firstExonStart_0base + 1;
		$secondExonStart_0base = $secondExonStart_0base + 1;
		$upstreamES = $upstreamES + 1 ;
		$downstreamES = $downstreamES + 1;
		if($strand eq "+"){
			return $upstreamES . ".." . $upstreamEE . "," . $secondExonStart_0base . ".." . $secondExonEnd . "," . $downstreamES . ".." . $downstreamEE;
		}else{
			return $downstreamES . ".." . $downstreamEE . "," . $secondExonStart_0base . ".." . $secondExonEnd . "," . $upstreamES . ".." . $upstreamEE;
		}
	}		
}


sub getAsStartAndEnd{
my ($asType, $strand, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $start, $end)=@_;
        if($asType eq "A5SS"){
                if($strand eq "+"){
                        $$start = $longExonStart_0base+1;
                        $$end = $flankingEE;
                }else{
                        $$start = $flankingES+1;
                        $$end = $longExonEnd
                }
        }elsif($asType eq "A3SS"){
                if($strand eq "+"){
                        $$start = $flankingES+1;
                        $$end = $longExonEnd;
                }else{
                        $$start = $longExonStart_0base + 1;
                        $$end = $flankingEE;
                }
        }
        if($asType eq "RI"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
        if($asType eq "SE"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
        if($asType eq "MXE"){
                $$start = $upstreamES + 1;
                $$end = $downstreamEE;
        }
}
