#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--originalKeggJson ko00001.json\\\n" .
                "--outputKeggTsv kegg.tsv \n";
	exit;
}
# download kegg json: https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=
# https://www.genome.jp/kegg-bin/get_htext?htext=ko00001&filedir=%2fkegg%2fbrite%2fko&length=
my ($originalKeggJson, $outputKeggTsv);

GetOptions(
        'originalKeggJson=s'=>\$originalKeggJson,
        'outputKeggTsv=s'=>\$outputKeggTsv,
);

open FF, "<$originalKeggJson";
#{
#  "name":"ko00001",
#  "children":[
#  {
#    "name":"09100 Metabolism",
#    "children":[
#    {
#      "name":"09101 Carbohydrate metabolism",
#      "children":[
#      {
#        "name":"00010 Glycolysis \/ Gluconeogenesis [PATH:ko00010]",
#        "children":[
#        {
#          "name":"K00844  HK; hexokinase [EC:2.7.1.1]"
#        },
#        {
#          "name":"K12407  GCK; glucokinase [EC:2.7.1.2]"
#        },

