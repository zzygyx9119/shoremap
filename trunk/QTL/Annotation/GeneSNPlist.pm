#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use Bio::Perl;
use SNPlist;

package GeneSNPlist;

sub new {
	my ($self, $snps) = @_;
	$self = {
		gene_id        => '',
		isoform        => 1, 
		chromosome     => '',
		start          => 0,
		end            => 0,
		cds_length     => 0,
		protein_length => 0,
		orientation    => '',
		SNP_lists      => $snps,
		coding_SNP     => {},
		NS_changes     => {},
		AA_changes     => {},
		CDSexon        => {},
		protein        => {},
		ref_gene       => '',
		ref_coding     => '',
		eco_coding     => '',
	};
	bless $self;
	return $self;
}

#############################################################################################
### Calculates all SNPs for specified gene and isoform
sub get_gene_snps {
	my ($self, $chromosome, $gene_start, $gene_end, $orientation, $gene_id, $isoform, $ref_gene, %cds_ann) = @_;
	$self->{chromosome}  = $chromosome;
	$self->{start}       = $gene_start;
	$self->{end}         = $gene_end;
	$self->{orientation} = $orientation;
	$self->{gene_id}     = $gene_id;
	$self->{isoform}     = $isoform;
	$self->{ref_gene}    = $ref_gene;

	# Get features of gene (range, CDS, RNA)
	foreach my $cds_start ( sort {$a <=> $b} keys %cds_ann ) {
		my $cds_seq = $cds_ann{$cds_start}[3];
		my $cds_end = $cds_ann{$cds_start}[1];
		
		$self->{CDSexon}{$cds_start} = $cds_end;
		$self->{cds_length} = $self->{cds_length} + $cds_end - $cds_start + 1;
		$self->{ref_coding} .= $cds_seq;
	}

	# reverse complement if orientation = "-"
	if($self->{orientation} eq "-") {
		$self->{ref_coding} = reverse($self->{ref_coding});
		$self->{ref_coding} =~ tr/ACGTUacgtu/TGCAATGCAA/;
	}

	# Protein length
	$self->{protein_length} = $self->{cds_length} / 3;

	return($self->{protein_length});
}


#########################################################################################
### Calculates CDS and protein changes for all SNPs in "SNP_lists" and fills SNP object
sub get_protein_changes {
	my ( $self ) = @_;
	$self->{protein}{'ref'} = Bio::Perl::translate_as_string($self->{ref_coding});

	### CDS changes
	my $cds_pos = 1;
	foreach my $x ( sort {$a <=> $b} keys %{$self->{CDSexon}} ) {
		for (my $i = $x; $i <= $self->{CDSexon}{$x}; $i++) {

			my $ref_base = substr($self->{ref_gene}, $i - $self->{start}, 1);
				
			if ( exists $self->{SNP_lists}{snps}{$i} ) {

				$self->{SNP_lists}{snps}{$i}{gene_id} = $self->{gene_id};
				$self->{SNP_lists}{snps}{$i}{stype} = "CDS";
					
				if($self->{orientation} eq "+") {
					$self->{coding_SNP}{$cds_pos} = $i;
					$self->{SNP_lists}{snps}{$i}{cds_pos} = $cds_pos;
					$self->{SNP_lists}{snps}{$i}{gene_pos} = $i - $self->{start} + 1 ;
				}
				elsif($self->{orientation} eq "-") {
					my $real_pos = $self->{cds_length} - $cds_pos + 1;
					$self->{coding_SNP}{$real_pos} = $i;
					$self->{SNP_lists}{snps}{$i}{cds_pos} = $real_pos;
					$self->{SNP_lists}{snps}{$i}{gene_pos} = $self->{end} - $i + 1;
				}
				$self->{eco_coding} .= $self->{SNP_lists}{snps}{$i}{new_base};
			}
			else {
				$self->{eco_coding} .= $ref_base;
			}
			$cds_pos++;
		}
	}

	if($self->{orientation} eq "-") {
		$self->{eco_coding} = reverse($self->{eco_coding});
		$self->{eco_coding} =~ tr/ACGTUacgtu/TGCAATGCAA/;
	}
	
	### Protein changes
	$self->{protein}{'alt'} = Bio::Perl::translate_as_string($self->{eco_coding});

	for(my $j = 1; $j <= length($self->{protein}{'alt'}); $j++) {
		my $ref_aa = substr($self->{protein}{'ref'}, $j - 1, 1);
		my $new_aa = substr($self->{protein}{'alt'}, $j - 1, 1);

		if( exists $self->{coding_SNP}{(3 * $j) - 2} ) {
			my $genome_position = $self->{coding_SNP}{(3 * $j) - 2};
			$self->{SNP_lists}{snps}{$genome_position}{codon_pos} = 1;
			$self->{SNP_lists}{snps}{$genome_position}{ref_aa} = $ref_aa;
			$self->{SNP_lists}{snps}{$genome_position}{new_aa} = $new_aa;

			if( $ref_aa ne $new_aa ) {
				$self->{AA_changes}{$j} = 1;
				$self->{NS_changes}{(3 * $j) - 2} = $genome_position;
				$self->{SNP_lists}{snps}{$genome_position}{ns_change} = 1;
				if   ($new_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{new_stop} = 1;}
				elsif($ref_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{lost_stop} = 1;}
			}
		}
				
		if( exists $self->{coding_SNP}{(3 * $j) - 1} ) { 
			my $genome_position = $self->{coding_SNP}{(3 * $j) - 1};
			$self->{SNP_lists}{snps}{$genome_position}{codon_pos} = 2;
			$self->{SNP_lists}{snps}{$genome_position}{ref_aa} = $ref_aa;
			$self->{SNP_lists}{snps}{$genome_position}{new_aa} = $new_aa;

			if( $ref_aa ne $new_aa ) {
				$self->{AA_changes}{$j} = 1;
				$self->{NS_changes}{(3 * $j) - 1} = $genome_position;
				$self->{SNP_lists}{snps}{$genome_position}{ns_change} = 1;
				if   ($new_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{new_stop} = 1;}
				elsif($ref_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{lost_stop} = 1;}
			}
		}
				
		if( exists $self->{coding_SNP}{3 * $j} ) {
			my $genome_position = $self->{coding_SNP}{3 * $j};
			$self->{SNP_lists}{snps}{$genome_position}{codon_pos} = 3;
			$self->{SNP_lists}{snps}{$genome_position}{ref_aa} = $ref_aa;
			$self->{SNP_lists}{snps}{$genome_position}{new_aa} = $new_aa;

			if( $ref_aa ne $new_aa ) {
				$self->{AA_changes}{$j} = 1;
				$self->{NS_changes}{3 * $j} = $genome_position;
				$self->{SNP_lists}{snps}{$genome_position}{ns_change} = 1;
				if   ($new_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{new_stop} = 1;}
				elsif($ref_aa eq '*') {$self->{SNP_lists}{snps}{$genome_position}{lost_stop} = 1;}
			}
		}
	}
	
	return(1);
}

1;
