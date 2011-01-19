#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package Indel;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => 0,
		begin         => 0,
		end           => 0,
		seq           => '',
		stype         => 'intergenic',
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		codon_pos     => 0,
		new_stop      => 0,
		lost_stop     => 0,
		splicechange  => 0,
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
	my ($self, $ecotype, $chromosome, $begin, $end, $seq, $support, $concordance, $quality) = @_;
	$self->{chromosome}  = $chromosome;
	$self->{begin}       = $begin;
	$self->{end}         = $end;
	$self->{seq}         = $seq;
	$self->{support}     = $support;
	$self->{concordance} = $concordance;
	$self->{quality}     = $quality;

}

1;
