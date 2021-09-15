#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \n\n";
        exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($firstExonSeries, $secondExonSeries);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
);

my ($dbh, $query, $update, $hash_ref, $start, $end);

#$dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);
my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$query=$dbh->prepare("select * from asTable");
$query->execute();

while($hash_ref=$query->fetchrow_hashref()){

	$start = 0;
	$end = 0;

	&getAsStartAndEnd($hash_ref->{"asType"}, $hash_ref->{"strand"}, $hash_ref->{"1stExonEnd"}, $hash_ref->{"1stExonStart_0base"}, $hash_ref->{"2ndExonEnd"}, $hash_ref->{"2ndExonStart_0base"}, $hash_ref->{"downstreamEE"}, $hash_ref->{"downstreamES"}, $hash_ref->{"exonEnd"}, $hash_ref->{"exonStart_0base"}, $hash_ref->{"flankingEE"}, $hash_ref->{"flankingES"}, $hash_ref->{"longExonEnd"}, $hash_ref->{"longExonStart_0base"}, $hash_ref->{"riExonEnd"}, $hash_ref->{"riExonStart_0base"}, $hash_ref->{"shortEE"}, $hash_ref->{"shortES"}, $hash_ref->{"upstreamEE"}, $hash_ref->{"upstreamES"}, \$start, \$end);

	$firstExonSeries = &generateFirstAltAsExonSeries($hash_ref->{"asType"}, $hash_ref->{"strand"}, $hash_ref->{"1stExonEnd"}, $hash_ref->{"1stExonStart_0base"}, $hash_ref->{"2ndExonEnd"}, $hash_ref->{"2ndExonStart_0base"}, $hash_ref->{"downstreamEE"}, $hash_ref->{"downstreamES"}, $hash_ref->{"exonEnd"}, $hash_ref->{"exonStart_0base"}, $hash_ref->{"flankingEE"}, $hash_ref->{"flankingES"}, $hash_ref->{"longExonEnd"}, $hash_ref->{"longExonStart_0base"}, $hash_ref->{"riExonEnd"}, $hash_ref->{"riExonStart_0base"}, $hash_ref->{"shortEE"}, $hash_ref->{"shortES"}, $hash_ref->{"upstreamEE"}, $hash_ref->{"upstreamES"});

	$secondExonSeries = &generateSecondAltAsExonSeries($hash_ref->{"asType"}, $hash_ref->{"strand"}, $hash_ref->{"1stExonEnd"}, $hash_ref->{"1stExonStart_0base"}, $hash_ref->{"2ndExonEnd"}, $hash_ref->{"2ndExonStart_0base"}, $hash_ref->{"downstreamEE"}, $hash_ref->{"downstreamES"}, $hash_ref->{"exonEnd"}, $hash_ref->{"exonStart_0base"}, $hash_ref->{"flankingEE"}, $hash_ref->{"flankingES"}, $hash_ref->{"longExonEnd"}, $hash_ref->{"longExonStart_0base"}, $hash_ref->{"riExonEnd"}, $hash_ref->{"riExonStart_0base"}, $hash_ref->{"shortEE"}, $hash_ref->{"shortES"}, $hash_ref->{"upstreamEE"}, $hash_ref->{"upstreamES"});

	$update = $dbh->prepare("update asTable set start=" . $start . ", end=" . $end . ", firstAltExonSeries = \"" . $firstExonSeries . "\", secondAltExonSeries=\"" . $secondExonSeries . "\" where asId=\"" . $hash_ref->{"asId"} . "\"");

	$update->execute();
}

$dbh->disconnect();

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
