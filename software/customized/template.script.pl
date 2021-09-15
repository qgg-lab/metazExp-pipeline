#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "-- \\\n" .
                "-- \\\n" .
                "-- \\\n" .
		"-- \n";
	exit;
}

my ($, $, $, $);

GetOptions(
        '=s'=>\$,
        '=s'=>\$,
        '=s'=>\$,
        '=s'=>\$,
	'=s'=>\$,
);


