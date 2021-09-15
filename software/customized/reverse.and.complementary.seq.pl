$st = $ARGV[0];
$st = reverse($st);
$st=~tr/ACGT/TGCA/;
print $st;
