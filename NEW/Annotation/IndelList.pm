#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use Indel;

package IndelList;

### Constructor
sub new {
	my $self = {
		indel_file   => "",
		chromosome   => "",
		region_start => 0,
		region_end   => 0,
		indels       => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $chromosome, $region_start, $region_end, $indel_file) = @_;
	$self->{chromosome}   = $chromosome;
	$self->{region_start} = $region_start;
	$self->{region_end}   = $region_end;
	$self->{indel_file}   = $indel_file;


	### Shore indel format:  sample | chromosome | begin | end | length | seq | support_type | support | concordance | avg_hits
	open INDELFILE, $indel_file or die "Cannot open indel file: " . $indel_file . "\n";
	while (<INDELFILE>) {
		my ($sample, $chr, $begin, $end, $length, $seq , $support_type, $support, $concordance, $avg_hits) = split "\t";

		if ( ($chr == $chromosome) && ($begin >= $region_start) && ($end <= $region_end) ) {
			$self->{indels}{$begin} = new Indel;
			$self->{indels}{$begin}->init( $sample, $chromosome, $begin, $end, $seq, $support, $concordance, 0);
		}
	}
	close INDELFILE;

	return(1);
}

1;
