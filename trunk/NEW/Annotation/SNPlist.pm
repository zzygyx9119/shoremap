#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use SNP;

package SNPlist;

### Constructor
sub new {
	my $self = {
		chromosome   => "",
		region_start => 0,
		region_end   => 0,
		snp_file     => "",
		snps         => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $chromosome, $region_start, $region_end, $snp_file) = @_;
	$self->{chromosome}   = $chromosome;
	$self->{region_start} = $region_start;
	$self->{region_end}   = $region_end;
	$self->{snp_file}     = $snp_file;


	### Shore SNP format:  sample | chromosome | position | refbase | snp  | support_type | support | concordance | max_qual | avg_hits
	open SNPFILE, $snp_file or die "Cannot open snp file: " . $snp_file . "\n";
	while (<SNPFILE>) {
		my ($sample, $chr, $position, $refbase, $snp, $qual, $support, $concordance, $avg_hits) = split "\t";

		if ( ($chr == $chromosome) && ($position >= $region_start) && ($position <= $region_end) ) {
			$self->{snps}{$position} = new SNP();
			$self->{snps}{$position}->init( $sample, $chromosome, $position, $refbase, $snp, $support, $concordance, $qual );
		}
	}
	close SNPFILE;

	return(1);
}

1;
