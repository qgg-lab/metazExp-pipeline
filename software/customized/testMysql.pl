#!/usr/bin/perl
use DBI;
#my $dbh = DBI.connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1");
my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
my $query=$dbh->prepare("create table asTable(asId char(10) primary key, chr char(50), strand char(1) default \"+\", start int default 0, end int default 0)");
$query->execute();
