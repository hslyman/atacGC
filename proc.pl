#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
#use DataBrowser;
use FAlite;


# gid to tid mapping
my %tid;
open(my $fh3, "ensembl_bt6_geneID_to_transcriptID.txt") or die;
my $h3 = <$fh3>;
while (<$fh3>) {
	chomp;
	my ($gid, $tid) = split;
	$tid{$gid} = $tid;
}

# fpkm
my %fpkm;
open(my $fh1, "genes.fpkm_tracking") or die;
my $h1 = <$fh1>;
my %count;
while (<$fh1>) {
	my @f = split;
	my $gid = $f[0];
	my $val = $f[9];
	$fpkm{$gid} = $val;
}

# bed
my %bed;
my %gene;
my $resolution = 100000;
open(my $fh2, "ensGene_promoter.bed") or die;
my $h2 = <$fh2>;
while (<$fh2>) {
	my ($chr, $beg, $end, $text) = split;
	my ($tid) = $text =~ /^(\w+)/;
	my $zip1 = int $beg / $resolution;
	my $zip2 = int $end / $resolution;
	$bed{$chr}{$zip1}{$tid} = {beg => $beg, end => $end};
	$bed{$chr}{$zip2}{$tid} = {beg => $beg, end => $end};
	$gene{$tid} = {chrom => $chr, beg => $beg, end => $end};
}

# genome
my %dna;
open(my $fh4, "genome.fa") or die;
my $fasta = new FAlite($fh4);
while (my $entry = $fasta->nextEntry) {
	my ($n) = $entry->def =~ /^>(\S+)/;
	$dna{"chr$n"} = $entry->seq;
#	last; # testing
}

# ATAC intersection
my %atac;
my %seen;
open(my $fh5, "atac_coverage.bg") or die;
while (<$fh5>) {
	chomp;
	my ($chr, $b1, $e1, $count) = split;
	my $zip1 = int $b1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip1}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip1}{$tid}{beg}, $bed{$chr}{$zip1}{$tid}{end});
		my $overlap = 0;
		if    ($b1 >= $b2 and $b1 <= $e2) {$overlap = 1}
		elsif ($e1 >= $b2 and $e1 <= $e2) {$overlap = 1}
		elsif ($b1 <= $b2 and $e1 >= $e2) {$overlap = 1}
		if ($overlap) {
			for (my $i = $b1; $i <= $e1; $i++) {
				$atac{$tid}{$i} = $count;
			}
		}
	}
	my $zip2 = int $e1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip2}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip2}{$tid}{beg}, $bed{$chr}{$zip2}{$tid}{end});
		my $overlap = 0;
		if    ($b1 >= $b2 and $b1 <= $e2) {$overlap = 1}
		elsif ($e1 >= $b2 and $e1 <= $e2) {$overlap = 1}
		elsif ($b1 <= $b2 and $e1 >= $e2) {$overlap = 1}
		if ($overlap) {
			for (my $i = $b1; $i <= $e1; $i++) {
				$atac{$tid}{$i} = $count;
			}
		}
	}
#	last if $chr ne 'chr1'; # testing
}

# final report
print "gid\tid\tchrom\tbeg\tend\tgc\tfpkm\tatac\n";
foreach my $gid (keys %fpkm) {
	if (not defined $tid{$gid}) {
		print STDERR "skipping $gid\n";
		next; # not sure why some are missing
	}

	my $tid = $tid{$gid};
	my $chr = $gene{$tid}{chrom};
	my $beg = $gene{$tid}{beg};
	my $end = $gene{$tid}{end};
	my $fpkm = $fpkm{$gid};

#	next if $chr ne 'chr1'; # testing

	# GC composition
	next if not defined $dna{$chr}; # ignoring the weird ones
	my $seq = substr($dna{$chr}, $beg -1, $end - $beg+1);
	my $gs = $seq =~ tr/G/G/;
	my $cs = $seq =~ tr/C/C/;
	my $gc = sprintf "%.1f", 100 * ($gs + $cs) / ($end - $beg + 1);
	
	# ATAC coverage
	my $total = 0;
	for (my $i = $beg; $i <= $end; $i++) {
		$total += $atac{$tid}{$i} if defined $atac{$tid}{$i};
	}
	my $atac = sprintf "%.2f", $total / ($end - $beg + 1);
	
	# output
	print "$gid\t$tid\t$chr\t$beg\t$end\t$gc\t$fpkm\t$atac\n";
}

__END__
