#!/usr/bin/perl
my $gff=$ARGV[0];
my $newGff = $ARGV[1];
# Chr01   phytozomev12    gene    31727   33514   .       -       .       ID=Podel.01G000300.v2.1;Name=Podel.01G000300
# Chr01   phytozomev12    mRNA    31727   33514   .       -       .       ID=Podel.01G000300.1.v2.1;Name=Podel.01G000300.1;pacid=37312175;longest=1;Parent=Podel.01G000300.v2.1
# Chr01   phytozomev12    exon    33052   33514   .       -       .       ID=Podel.01G000300.1.v2.1.exon.1;Parent=Podel.01G000300.1.v2.1;pacid=37312175
# Chr01   phytozomev12    five_prime_UTR  33052   33514   .       -       .       ID=Podel.01G000300.1.v2.1.five_prime_UTR.1;Parent=Podel.01G000300.1.v2.1;pacid=37312175
# Chr01   phytozomev12    exon    31909   32004   .       -       .       ID=Podel.01G000300.1.v2.1.exon.2;Parent=Podel.01G000300.1.v2.1;pacid=37312175
# Chr01   phytozomev12    CDS     31909   31993   .       -       0       ID=Podel.01G000300.1.v2.1.CDS.1;Parent=Podel.01G000300.1.v2.1;pacid=37312175

# >Sevir.4G000200.1 
# pep seq
# >Sevir.4G009601.1.p
# GTF:
# （1）mRNA改成transcript
# （2）Id=* => gene_id "xxxx"; transcript_id "yyyy"
#  (3) geneId 从MD00G1000100.v1.1 => MD00G1000100
#  (4) transcriptId 从MD00G1000100.v1.1.491 => MD00G1000100.491
my ($line, @field, $geneId, $trsptId);
open FF, "<$gff";
open WW, ">$newGff";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	if($field[2] eq "gene"){
		if($field[8]=~/ID=(.*?)\.v2\.1;Name=(.*)/){
			$geneId = $1 . "g";
			pop(@field);
			print WW join("\t", @field, "gene_id \"" . $geneId . "\";") . "\n";
		}
	}elsif($field[2] eq "mRNA"){
		if($field[8]=~/ID=(.*?)\.v2\.1;Name.*/){
			pop(@field);
			$trsptId = $1;
			$field[2] = "transcript";
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}elsif($field[2] eq "exon"){
		if($field[8]=~/ID=(.*?)\.v2\.1\.exon\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\"; exon_number \"" . $exonNum . "\";") . "\n";
		}
	}elsif($field[2] eq "CDS"){
		if($field[8]=~/ID=(.*?)\.v2\.1\.CDS\.(\d+);Parent=.*/){
			pop(@field);
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\"; exon_number \"" . $exonNum . "\";") . "\n";
		}
	}elsif($field[2] eq "five_prime_UTR"){
		# ID=MD00G1000200.v1.1.491.five_prime_UTR.1;Parent=MD00G1000200.v1.1.491;pacid=40123412
		if($field[8]=~/ID=(.*?)\.v2\.1\.five_prime_UTR\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}elsif($field[2] eq "three_prime_UTR"){
		# ID=MD00G1000200.v1.1.491.three_prime_UTR.1;Parent=MD00G1000200.v1.1.491;pacid=40123412
		if($field[8]=~/ID=(.*?)\.v2\.1\.three_prime_UTR\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}

}
close FF;
close WW;
