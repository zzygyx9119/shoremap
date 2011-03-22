#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package SNP;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => 0,
		position      => 0,
		stype         => 'intergenic',
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		codon_pos     => 0,
		ns_change     => 0,
		new_stop      => 0,
		lost_stop     => 0,
		splicechange  => 0,
		ref_base      => '',
		new_base      => '',
		ref_aa        => '',
		new_aa        => '',
		domain_change => {},
		support       => 0,
		concordance   => 0,
		quality       => 0,
		peak_distance => 999999999,
		marker_ratio  => 0,
		source        => 'IlluminaGA2',
	};
	bless $self;
	return $self;
}

sub init
{
	my ($self, $ecotype, $chromosome, $position, $ref_base, $new_base, $support, $concordance, $quality) = @_;
	$self->{chromosome}  = $chromosome;
	$self->{position}    = $position;
	$self->{ref_base}    = $ref_base;
	$self->{new_base}    = $new_base;
	$self->{support}     = $support;
	$self->{concordance} = $concordance;
	$self->{quality}     = $quality;
}

1;
