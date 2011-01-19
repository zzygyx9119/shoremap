#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $file     = "";	# minor_allele_frequency file (will be replaced by new SNP calling file soon)
my $ref      = "";
my $chrsizes = "";
my $support  = 4;
my $freq     = 0.15;
my $win      = 200000;
my $win_step = 10000;

my %CMD;
GetCom();

my $dist_file = "SHOREmap.markerless.support$support.freq$freq.win$win.distance.txt";
my $slw_file  = "SHOREmap.markerless.support$support.freq$freq.win$win.txt";
my $pdffile   = "SHOREmap.markerless.support$support.freq$freq.win$win.pdf";

### Parse heterozygous predictions from snp file --------------------------------------------
my %chr = ();
my %het_freq = ();
open FILE, $file or die "Cannot open $file\n";
while (<FILE>) {
 	my @a = split "\t", $_;

	# Add new chromosome
	if(! exists $chr{$a[1]}) {
		my @het_array = ();
		my %het_hash  = ();

		$chr{$a[1]} = \@het_array;
		$het_freq{$a[1]} = \%het_hash;
	}

	# Quality filter for marker
	if( 	($a[8] >= $support) &&			# min support
		($a[9] >= $freq) &&			# min frequency of minor allele
		($a[4] ne "-" && $a[7] ne "-")		# No gaps
	) {
		push( @{$chr{$a[1]}}, $a[2] );
		
		if($a[3] eq $a[4]) {			# major allele is ref
			$het_freq{$a[1]}{$a[2]} = $a[9];
		}
		else {					# major allele is not ref
			$het_freq{$a[1]}{$a[2]} = $a[6];
		}
	}
}
close FILE;


### Calculate position wise distance to closest het prediction -----------------------------
my $start = 0;
open OUTFILE, ">$dist_file" or die "Cannot open output file.\n";
foreach my $c ( sort {$a <=> $b} keys %chr ) {

	for(my $i = 0; $i < scalar(@{$chr{$c}}); $i++) {

		for(my $pos = $start + 1; $pos <= $chr{$c}[$i]; $pos++) {
			my $dist_left = $pos - $start;
			my $dist_right = $chr{$c}[$i] - $pos;
			
			if($dist_left < $dist_right) {
				print OUTFILE "$c\t$pos\t$dist_left\n";
			}
			else {
				print OUTFILE "$c\t$pos\t$dist_right\n";
			}
		}
		$start = $chr{$c}[$i];
	}
}
close OUTFILE;

print "FINISHED Calculating distance\n";

### Sliding window ref-call analysis ------------------------------------------------------
my $chr = "NA";
my $genome_pos = 1;
my $front_start = $win_step - $win;
my %ref_density = ();
my %ref_sums = ();

open REF, $ref or die "Cannot open ref file.\n";
while( <REF> ) {
	chomp;
	my @a = split("\t", $_);

	### New chromosome reached
	if( $chr ne $a[1] ) {

		# First chromosome: Initialize
		if($chr eq "NA") {
			$chr = $a[1];
			$genome_pos = 1;
			$front_start = -$win;
			my %chr_hash = ();
			$ref_density{$chr} = \%chr_hash;
		}

		# Later chromosomes
		else {
			# Print last sliding window of chromosome
			if($front_start >= 0) {
				$ref_sums{$front_start} /= $win;
				$ref_density{$chr}{$front_start} = $ref_sums{$front_start};
			}

			# Reset
			$chr = $a[1];
			$genome_pos = 1;
			$front_start = $win_step - $win;
			%ref_sums = ();
		}
	}

	### Walk along the chromosome
	while( $genome_pos < $a[2] ) {
		$genome_pos++;

		# Sliding window end reached: Store and remove window
		if(  ($genome_pos % $win_step) == 0 ) {
			if($front_start >= 0) {
				$ref_sums{$front_start} /= $win;
				$ref_density{$chr}{$front_start} = $ref_sums{$front_start};
			}
			
			delete $ref_sums{$front_start};
			$front_start += $win_step;
		}
	}

	### Update ref_sums
	for( my $i = $front_start; $i < $front_start + $win; $i += $win_step) {
		$ref_sums{$i}++;
	}
}
close(REF);


print "FINISHED READING REF\n";

### Sliding window het distance and density analysis --------------------------------------------------
open OUTFILE, "$dist_file" or die "Cannot open $dist_file file.\n";
open SLW, ">$slw_file" or die "Cannot open file $slw_file.\n";

$chr = "NA";
$genome_pos = 1;
$front_start = $win_step - $win;
my %dist_sum = ();
my %het_freq_sum = ();

while( <OUTFILE> ) {
	chomp;
	my @a = split("\t", $_);
	$genome_pos++;

	### New chromosome reached
	if( $chr ne $a[0] ) {

		print "NOW STARTRTING CHR $a[0]\n";

		# First chromosome
		if($chr eq "NA") {
			$chr = $a[0];
			$genome_pos = 1;
			$front_start = -$win;
		}

		# Later chromosomes
		else {
			if($front_start >= 0) {
				if( (! exists $het_freq_sum{$front_start}) || ($het_freq_sum{$front_start} < 1) ) { $het_freq_sum{$front_start} = 1; }
				my $win_center = $front_start + ($win / 2);
				my $dist_avg = $dist_sum{$front_start} / $het_freq_sum{$front_start} / $win * $ref_density{$chr}{$front_start} * $ref_density{$chr}{$front_start};
				print SLW "$chr\t$win_center\t$dist_sum{$front_start}\t$het_freq_sum{$front_start}\t$dist_avg\n";
			}

			$chr = $a[0];
			$genome_pos = 1;
			$front_start = $win_step - $win;
			%dist_sum = ();
			%het_freq_sum = ();
		}
	}


	### Sliding window end reached
	if( ($genome_pos % $win_step) == 0 )  {
		if($front_start >= 0) {
			if( (! exists $het_freq_sum{$front_start}) || ($het_freq_sum{$front_start} < 1) ) { $het_freq_sum{$front_start} = 1; }
			my $win_center = $front_start + ($win / 2);
			my $dist_avg = $dist_sum{$front_start} / $het_freq_sum{$front_start} / $win * $ref_density{$chr}{$front_start} * $ref_density{$chr}{$front_start};
			print SLW "$chr\t$win_center\t$dist_sum{$front_start}\t$het_freq_sum{$front_start}\t$dist_avg\n";
		}

		delete $dist_sum{$front_start};
		delete $het_freq_sum{$front_start};
		$front_start += $win_step;
	}

	### Update sums
	for( my $i = $front_start; $i < $front_start + $win; $i += $win_step) {
		$dist_sum{$i} += $a[2];
		if( $a[2] == 0 ) {
			$het_freq_sum{$i} += $het_freq{$chr}{$a[1]};
		}
	}
}
close OUTFILE;
close SLW;



### Call R for plotting -----------------------------------------------------------------------
my $cmd = "R --slave --vanilla --args $chrsizes $pdffile $slw_file $win_step $win < $FindBin::Bin/SHOREmap.R";
print STDERR $cmd, "\n";
system($cmd);

my $del_cmd = "rm $dist_file";
system($del_cmd);


exit(0);



### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory:
--snp        STRING     SNP file in Shore format (\"minor_allele_position.txt\")
--refseq     STRING     Reference calls in Shore format (\"reference.txt\")
--chrsizes   STRING     Chromosome sizes file

Optional:
--support    INT        Minimum support for both alleles of a het call       default=4
--freq       INT        Minimum frequency for both alleles of a het call     default=0.15
--winsize    INT        Sliding window size for plot                         default=200000

See documentation for file formats.

SHOREmap is written by Korbinian Schneeberger and Stephan Ossowski.
Max Planck Institute for Developemental Biology, Tübingen, 2009.

\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "snp=s", "refseq=s", "chrsizes=s", "support=s", "freq=s", "winsize=s");

	# Mandatory params
	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify refseq file\n") unless defined($CMD{refseq});
	die("Please specify chr sizes file\n") unless defined($CMD{chrsizes});
	$file = $CMD{snp};
	$ref  = $CMD{refseq};
	$chrsizes = $CMD{chrsizes};

	# Optional params
	if(defined $CMD{support}) { $support = $CMD{support}; }
	if(defined $CMD{freq}) { $freq = $CMD{freq}; }
	if(defined $CMD{winsize}) { $win = $CMD{winsize}; }
}

