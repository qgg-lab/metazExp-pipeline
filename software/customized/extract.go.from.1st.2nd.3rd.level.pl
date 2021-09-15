#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--goDic go-basic.obo/go.obo\\\n" .
                "--outputGoLevel go.with.level.tsv \n";
	exit;
}

my ($goDic, $levelNum, $outputGoLevel);

GetOptions(
        'goDic=s'=>\$goDic,
        'outputGoLevel=s'=>\$outputGoLevel,
);

# 读取GO字典到内存文本中
# id: GO:0000002
# name: mitochondrial genome maintenance
# namespace: biological_process
# def: "The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome." [GOC:ai, GOC:vw]
# is_a: GO:0007005 ! mitochondrion organization
my (%go, $goHref, $goId);
$goHref = \%go;
my $goDicText=`cat $goDic`;
my @term = split(/\[Term\]\n/, $goDicText);
shift(@term);
foreach my $term(@term){
	my @line = ();
	@line = split(/\n/, $term);
	foreach my $line(@line){
		if($line=~/^id: (GO:\d+)/){
			$goId = $1;
			$goHref->{$goId}->{"output"} = 1;
		}elsif($line=~/^name: (.*)/){
			$goHref->{$goId}->{"name"} = $1;
		}elsif($line=~/^namespace: (.*)/){
			$goHref->{$goId}->{"namespace"} = $1;
		}elsif($line=~/^is_a: (GO:\d+) .*/){
			$goHref->{$goId}->{"childList"} .= $1 . ",";
			$goHref->{$1}->{"parentList"} .= $goId . ",";
		}elsif($line=~/^relationship: part_of (GO:\d+)/){
			$goHref->{$goId}->{"childList"} .= $1 . ",";
			$goHref->{$1}->{"parentList"} .= $goId . ",";
		}
	}
}

open WW, ">$outputGoLevel";
# 输出第1层GO
my @goId = keys(%go);
foreach $goId(@goId){
	if(not exists($goHref->{$goId}->{"parentList"})){
		print WW join("\t", $goId, $goHref->{$goId}->{"namespace"}, $goHref->{$goId}->{"name"}, "1") . "\n";
		$goHref->{$goId}->{"output"} = $goHref->{$goId}->{"output"} - 1;
	}
}

# 输出第2层GO
@goId = ();
@goId = keys(%go);
foreach $goId(@goId){
	# 当前的GO没有输出，但是它的双亲有被输出了（如果多个双亲，那么只要有1个双亲输出都属于第2层）
	if($goHref->{$goId}->{"output"}==1 and &checkParentOutput($goHref, $goHref->{$goId}->{"parentList"}) == 1){
		print WW join("\t", $goId, $goHref->{$goId}->{"namespace"}, $goHref->{$goId}->{"name"}, "2") . "\n";
		$goHref->{$goId}->{"output"} = $goHref->{$goId}->{"output"} - 1;
	}	
}

# 输出第3层GO
@goId = ();
@goId = keys(%go);
foreach $goId(@goId){
	# 当前的GO没有输出，但是它的双亲有被输出了（如果多个双亲，那么只要有1个双亲输出都属于第3层）
	if($goHref->{$goId}->{"output"}==1 and &checkParentOutput($goHref, $goHref->{$goId}->{"parentList"}) == 1){
		print WW join("\t", $goId, $goHref->{$goId}->{"namespace"}, $goHref->{$goId}->{"name"}, "3") . "\n";
		$goHref->{$goId}->{"output"} = $goHref->{$goId}->{"output"} - 1;
	}
}

# 输出第4层GO
@goId = ();
@goId = keys(%go);
foreach $goId(@goId){
	# 当前的GO没有输出，但是它的双亲有被输出了（如果多个双亲，那么只要有1个双亲输出都属于第3层）
	if($goHref->{$goId}->{"output"}==1 and &checkParentOutput($goHref, $goHref->{$goId}->{"parentList"}) == 1){
		print WW join("\t", $goId, $goHref->{$goId}->{"namespace"}, $goHref->{$goId}->{"name"}, "4") . "\n";
		$goHref->{$goId}->{"output"} = $goHref->{$goId}->{"output"} - 1;
	}
}
sub checkParentOutput{
	my ($goRef, $parentList) = @_;
#	print "parentList:" . $parentList;
#	<STDIN>;
	my @parent = split(/,/, $parentList);
	foreach my $parent(@parent){
		if($goRef->{$parent}->{"output"} == 0){
			return 1;
		}
	}
	return 0;
}
