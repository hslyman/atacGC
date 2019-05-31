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

# CpG island intersection
my %cpgi;
open(my $fh4, "cpgIslands_bosTau6.bed") or die;
while(<$fh4>) {
	chomp;
	my ($chr, $b1, $e1) = split;
	my $zip1 = int $b1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip1}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip1}{$tid}{beg}, $bed{$chr}{$zip1}{$tid}{end});
		$cpgi{$tid} = bed_overlap($b1, $e1, $b2, $e2);
	}
	my $zip2 = int $e1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip2}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip2}{$tid}{beg}, $bed{$chr}{$zip2}{$tid}{end});
		$cpgi{$tid} = bed_overlap($b1, $e1, $b2, $e2);
	}
}
close $fh4;

# ATAC intersection
my %atac;
open(my $fh5, "lib08_peaks.broadPeak") or die;
while (<$fh5>) {
	chomp;
	my ($chr, $b1, $e1, $name, $score, $strand, $signal, $p, $q) = split;
	my $zip1 = int $b1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip1}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip1}{$tid}{beg}, $bed{$chr}{$zip1}{$tid}{end});
		$atac{$tid} = bed_overlap($b1, $e1, $b2, $e2);
	}
	my $zip2 = int $e1 / $resolution;
	foreach my $tid (keys %{$bed{$chr}{$zip2}}) {
		my ($b2, $e2) = ($bed{$chr}{$zip2}{$tid}{beg}, $bed{$chr}{$zip2}{$tid}{end});
		$atac{$tid} = bed_overlap($b1, $e1, $b2, $e2); 
	}
}

# final report
print "gid\tid\tchrom\tbeg\tend\tcpgi\tfpkm\tatac\n";
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

	# ATAC peak
	my $atac;
	if (not exists $atac{$tid}) {
		$atac = 0;
	}
	else {
		$atac = 1;
	}
	
	# CpG islands
	my $cpgi;
	if (not exists $cpgi{$tid}) {
		$cpgi = 0;
	}
	else {
		$cpgi = $cpgi{$tid};
	}

	# output
	print "$gid\t$tid\t$chr\t$beg\t$end\t$cpgi\t$fpkm\t$atac\n";
}


sub bed_overlap {
	my ($b1, $e1, $b2, $e2) = @_;
	my $overlap = 0;
	if	($b1 >= $b2 and $b1 <= $e2) {$overlap = 1}
	elsif	($e1 >= $b2 and $e1 <= $e2) {$overlap = 1}
	elsif	($b1 <= $b2 and $e1 >= $e2) {$overlap = 1}
	return($overlap);
}
__END__
