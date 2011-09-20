#!/usr/bin/perl -w

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
use lib $FindBin::Bin;

use SNPlist;
use IndelList;
use GeneSNPlist;


### Variables
my $refseq_file = "";
my $refseq_error_file = "";
my $gff = "";
my $snp_file = "";
my $del_file = "";
my $ins_file = "";
my $dist_file = "";
my $chromosome;
my $start;
my $end;


### Get command line options
my %CMD;
GetCom();


### Get SNP lists
my $snps = new SNPlist();
$snps->get($chromosome, $start, $end, $snp_file);

### Get deletion list
my $deletions = new IndelList();
if($del_file ne "") { $deletions->get($chromosome, $start, $end, $del_file); }

### Get insertion list
my $insertions = new IndelList();
if($ins_file ne "") { $insertions->get($chromosome, $start, $end, $ins_file); }


### Parse marker distribution file to get peak position and heigth
my $peak_height = 0;
my $peak_pos = 0;
my %dist = ();

if($dist_file ne "") {
	open DIST, $dist_file or die "Cannot open marker distribution file.\n";
	while( <DIST> ) {
		my ($current_chr, $current_pos, $ref_count, $alt_count, $ratio) = split(" ", $_);

		if( ($current_chr eq $chromosome) && ($current_pos >= $start) && ($current_pos <= $end) ) {
			if(abs($ratio) > $peak_height) {
				$peak_height = abs($ratio);
				$peak_pos = $current_pos;
			}

			$dist{$current_pos} = $ratio;
		}
	}
}
else {
	$peak_pos = 1;
}
print STDERR "\n\nMax peak in interval: ". $peak_pos ."\n\n";


### Parse known reference sequence errors
my %ref_err = ();
if($refseq_error_file ne "") {
	open REFERR, $refseq_error_file or die "Cannot open refseq error file\n";
	while( <REFERR> ) {
		chomp;
		my ($chr, $pos) = split("\t", $_);
		if(! exists $ref_err{$chr}) {
			my %chr_hash = ();
			$ref_err{$chr} = \%chr_hash;
		}
		$ref_err{$chr}{$pos} = 1;
	}
}

### Prioritize SNPs
foreach my $snp_pos ( sort {$a <=> $b} keys %{$snps->{snps}} ) {
	$snps->{snps}{$snp_pos}{peak_distance} = abs($snp_pos - $peak_pos);
}


### Prioritize Deletions
foreach my $del_pos ( sort {$a <=> $b} keys %{$deletions->{indels}} ) {
	$deletions->{indels}{$del_pos}{peak_distance} = abs($del_pos - $peak_pos);
}

# Prioritize Insertions
foreach my $ins_pos ( sort {$a <=> $b} keys %{$insertions->{indels}} ) {
	$insertions->{indels}{$ins_pos}{peak_distance} = abs($ins_pos - $peak_pos);
}


### Functional analysis of SNPs and indels
my %coding_ann = ();
my %gene_ann = ();
my %seq_type = ();

if( ($gff ne "") && ($refseq_file ne "") ) {

	### Get gene annotation from gff
	open GFF, $gff or die "Cannot open gff file\n";
	while( <GFF> ) {

		# GFF format: chr, source, seq_type, start, end, score, orientation, frame, description
		my @columns = split(" ", $_);

		# Quick fix for A.thaliana TAIR8 annotation, disable for other species
		if($columns[0] =~ /Chr/) { $columns[0] = substr($columns[0], 3); }

		# Check if gene is in range
		if( ($columns[0] eq $chromosome) && ($columns[4] > $start) && ($columns[3] < $end) ) {

			# Split decription column
			my @desc = split(";", $columns[8]);

			# Get gene locus
			if( $columns[2] eq "gene" ) {
				if( substr($desc[0], 0, 2) eq "ID") {
					my ($junk, $gene_name)  = split("=", $desc[0]);
					if(! exists $coding_ann{$gene_name}) {
						#my %gene = ();
						#$coding_ann{$gene_name} = \%gene;

						my @gene_locus = ($columns[3], $columns[4], $columns[6], "");
						$gene_ann{$gene_name} = \@gene_locus;

						my %type = ();
						$seq_type{$gene_name} = \%type;
						$seq_type{$gene_name}{"gene"} = $columns[3] . "-" . $columns[4];
					}
				}
			}

			# Get mRNA locus
			elsif( $columns[2] eq "mRNA" ) {
				# TODO, currently not needed
			}

			# Get coding sequence
			elsif( $columns[2] eq "CDS" ) {
				my ($parent, $protein_name) = split(",", $desc[0]);
				my ($junk, $gene_id) = split("=", $parent);
				my $isoform = 1;
				my $gene_name = $gene_id;

				if( $gene_id =~ /\./ ) {
					($gene_name, $isoform) = split(/\./, $gene_id);
				}

				if($isoform == 1) {
					if(! exists $coding_ann{$gene_name}) {
						my %gene = ();
						$coding_ann{$gene_name} = \%gene;
					}

					my @cds = ($columns[3], $columns[4], $columns[6], "");
					$coding_ann{$gene_name}{$columns[3]} = \@cds;

					$seq_type{$gene_name}{"CDS"} = $columns[3] . "-" . $columns[4];
				}
			}

			# Get other sequence types
			elsif( ($columns[2] eq "exon") || ($columns[2] eq "five_prime_UTR") || ($columns[2] eq "three_prime_UTR") ) {
				my ($junk, $gene_id) = split("=", $desc[0]);
				my $isoform = 1;
				my $gene_name = $gene_id;
				if( $gene_id =~ /\./ ) {
					($gene_name, $isoform) = split(/\./, $gene_id);
				}

				if($isoform == 1) {
					$seq_type{$gene_name}{$columns[2]} = $columns[3] . "-" . $columns[4];
				}
			}
		}
	}


	### Load chromosome
	my $chr_seq = "";
	open GENOME, $refseq_file or die "Cannot open reference sequence file\n";
	while( <GENOME> ) {
		chomp;
		if(substr($_, 0, 1)  eq ">") {
			my $current_chr = substr($_, 1);
			$current_chr =~ s/\s.*//g;

			if($current_chr eq $chromosome) {
				while(<GENOME>) {
					chomp;
					if(substr($_, 0, 1)  eq ">") { 
						last; 
					}
					else {
						$chr_seq .= $_;
					}
				}
			}
		}

		if($chr_seq ne "" ) { last; }
	}

	### Extract gene seq
	foreach my $gene_name ( sort keys %coding_ann ) {
		$gene_ann{$gene_name}[3] = substr($chr_seq, $gene_ann{$gene_name}[0] - 1, $gene_ann{$gene_name}[1] - $gene_ann{$gene_name}[0] + 1);

		foreach my $start ( sort {$a <=> $b} keys %{$coding_ann{$gene_name}} ) {
			$coding_ann{$gene_name}{$start}[3] = substr($chr_seq, $start - 1, $coding_ann{$gene_name}{$start}[1] - $start + 1);
		}
	}

	### Add gene SNP and Indel annotation
	foreach my $gene_name ( sort keys %gene_ann ) {
		
		# Coding SNP/Indel
		if( exists $coding_ann{$gene_name} ) {
			my $gene = new GeneSNPlist($snps);
			my $results = $gene->get_gene_snps( $chromosome, $gene_ann{$gene_name}[0], $gene_ann{$gene_name}[1],
							$gene_ann{$gene_name}[2], $gene_name, 1, $gene_ann{$gene_name}[3], %{$coding_ann{$gene_name}});
			$gene->get_protein_changes();

		}

		# Noncoding SNP/Indels
		foreach my $seq_type ( keys %{$seq_type{$gene_name}} ) {
			my ($seq_type_start, $seq_type_end) = split( "-", $seq_type{$gene_name}{$seq_type} );

			# SNPs
			foreach my $snp_pos ( sort {$a <=> $b} keys %{$snps->{snps}} ) {

				if( ($snp_pos >= $seq_type_start) && ($snp_pos <= $seq_type_end) ) {

					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$snps->{snps}{$snp_pos}{stype} = $seq_type;
					}

					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $snps->{snps}{$snp_pos}{stype} ne "CDS") {
							$snps->{snps}{$snp_pos}{stype} = $seq_type;
							$snps->{snps}{$snp_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic/noncoding
					elsif( $seq_type eq "gene" ) {
						if(	($snps->{snps}{$snp_pos}{stype} ne "CDS") && 
							($snps->{snps}{$snp_pos}{stype} ne "five_prime_UTR") && 
							($snps->{snps}{$snp_pos}{stype} ne "three_prime_UTR")
						){
							$snps->{snps}{$snp_pos}{stype} = "intronic/noncoding";
							$snps->{snps}{$snp_pos}{gene_id} = $gene_name;
						}
					}
				}
			}

			# Deletions
			foreach my $del_pos ( sort {$a <=> $b} keys %{$deletions->{indels}} ) {
				if( ($del_pos >= $seq_type_start) && ($del_pos <= $seq_type_end) ) {
					
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$deletions->{indels}{$del_pos}{stype} = $seq_type;
					}

					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $deletions->{indels}{$del_pos}{stype} ne "CDS") {
							$deletions->{indels}{$del_pos}{stype} = $seq_type;
							$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic/noncoding
					elsif( $seq_type eq "gene" ) {
						if(     ($deletions->{indels}{$del_pos}{stype} ne "CDS") &&
							($deletions->{indels}{$del_pos}{stype} ne "five_prime_UTR") &&
							($deletions->{indels}{$del_pos}{stype} ne "three_prime_UTR")
						){
							$deletions->{indels}{$del_pos}{stype} = "intronic/noncoding";
							$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
						}
					}
				}
			}

			# Insertions
			foreach my $ins_pos ( sort {$a <=> $b} keys %{$insertions->{indels}} ) {
				if( ($ins_pos >= $seq_type_start) && ($ins_pos <= $seq_type_end) ) {
					
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$insertions->{indels}{$ins_pos}{stype} = $seq_type;
					}


					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $insertions->{indels}{$ins_pos}{stype} ne "CDS") {
							$insertions->{indels}{$ins_pos}{stype} = $seq_type;
							$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic/noncoding
					elsif( $seq_type eq "gene" ) {
						if(	($insertions->{indels}{$ins_pos}{stype} ne "CDS") &&
							($insertions->{indels}{$ins_pos}{stype} ne "five_prime_UTR") &&
							($insertions->{indels}{$ins_pos}{stype} ne "three_prime_UTR")
						){
							$insertions->{indels}{$ins_pos}{stype} = "intronic/noncoding";
							$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
						}
					}
				}
			}
		}
	}
}


### Print SNP results
open SNPOUT, ">snp.priority.txt" or die "Cannot open SNP output file\n";
foreach my $pos (sort {$snps->{snps}{$a}{peak_distance} <=> $snps->{snps}{$b}{peak_distance}} keys %{$snps->{snps}} ) {
	print SNPOUT "$chromosome\t$pos\t" . 
		$snps->{snps}{$pos}{ref_base} . "\t" . $snps->{snps}{$pos}{new_base} . "\t" . 
		$snps->{snps}{$pos}{peak_distance} . "\t" . $snps->{snps}{$pos}{support} . "\t" . 
		sprintf("%.2f", $snps->{snps}{$pos}{concordance}) ."\t" . $snps->{snps}{$pos}{quality};

	if(exists $ref_err{$chromosome}{$pos}) { print SNPOUT "\tREFERR"; }
	else { print SNPOUT "\tNEWSNP"; }
	
	if( ($gff ne "") && ($refseq_file ne "") ) {
		print SNPOUT "\t" . $snps->{snps}{$pos}{stype};

		if($snps->{snps}{$pos}{gene_id} ne "NA") {
			print SNPOUT "\t" . $snps->{snps}{$pos}{gene_id} . "\t1";
		}

		if($snps->{snps}{$pos}{cds_pos} != 0) {
			my $syn_nonsyn = "Syn";
			if($snps->{snps}{$pos}{ns_change} == 1) { $syn_nonsyn = "Nonsyn"; }

			print SNPOUT "\t" . $snps->{snps}{$pos}{cds_pos} . "\t" . $snps->{snps}{$pos}{codon_pos} . "\t$syn_nonsyn\t" . 
				$snps->{snps}{$pos}{ref_aa} . "\t" . $snps->{snps}{$pos}{new_aa};
		}
	}

	print SNPOUT "\n";
}
close SNPOUT;


### Print deletion results
if($del_file ne "") {
	open DELOUT, ">deletion.priority.txt" or die "Cannot open deletion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$deletions->{indels}} ) {
		
		print DELOUT "$chromosome\t$start\t" . $deletions->{indels}{$start}{end} ."\t". $deletions->{indels}{$start}{seq} ."\t". 
				$deletions->{indels}{$start}{peak_distance} ."\t". $deletions->{indels}{$start}{support} ."\t". 
				sprintf("%.2f", $deletions->{indels}{$start}{concordance}) ."\t". $deletions->{indels}{$start}{quality};
				
		if( ($gff ne "") && ($refseq_file ne "") ) {	
			print DELOUT "\t". $deletions->{indels}{$start}{stype};

			if( $deletions->{indels}{$start}{gene_id} ne "NA") {
				print DELOUT "\t". $deletions->{indels}{$start}{gene_id}. "\t1";
			}
		}
	
		print DELOUT "\n";
	}

	close DELOUT;
}


### Print insertion results
if($ins_file ne "") {
	open INSOUT, ">insertion.priority.txt" or die "Cannot open insertion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$insertions->{indels}} ) {

		print INSOUT "$chromosome\t$start\t" . $insertions->{indels}{$start}{end} ."\t". $insertions->{indels}{$start}{seq} ."\t".
				$insertions->{indels}{$start}{peak_distance} ."\t". $insertions->{indels}{$start}{support} ."\t". 
				sprintf("%.2f", $insertions->{indels}{$start}{concordance}) ."\t". $insertions->{indels}{$start}{quality};
				
		if( ($gff ne "") && ($refseq_file ne "") ) {
			print INSOUT "\t". $insertions->{indels}{$start}{stype};

			if( $insertions->{indels}{$start}{gene_id} ne "NA") {
				print INSOUT "\t". $insertions->{indels}{$start}{gene_id}. "\t1";
			}
		}

		print INSOUT "\n";
	}
	
	close INSOUT;
}

exit(0);


### Read command line parameters
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory:
--snp      STRING     SNP file in Shore format (\"quality_variant.txt\")
--chrom    STRING     Chromosome of target region
--start    INT        Start of target region
--end      INT        End of target region

Rank by peak distance (optional):
--dist     STRING     Output file of either SHOREmap_interval.pl or SHOREmap_denovo.pl
                      If file is not specified peak position is set to 1.

Functional SNP and indel annotation (optional, only used if both are specified):
--genome   STRING     Reference sequence file (chromosome names have to be equal to SNP file)
--gff      STRING     Gene annotation in GFF format

Correct for known reference sequence errors (optional):
--referr   STRING     Reference sequence errors (chr | position)

Small indel annotation (optional):
--del      STRING     Deletion file in Shore format
--ins      STRING     Insertion file in Shore format
\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "snp=s", "dist=s", "chrom=s", "start=s", "end=s", "genome=s", "gff=s", "referr=s", "del=s", "ins=s");

	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify chromosome of target region\n") unless defined($CMD{chrom});
	die("Please specify start of target region\n") unless defined($CMD{start});
	die("Please specify end of target region\n") unless defined($CMD{end});

	$snp_file    = $CMD{snp};
	$chromosome  = $CMD{chrom};
	$start       = $CMD{start};
	$end         = $CMD{end};

	$dist_file = "";
	$gff = "";
	$refseq_file = "";
	$refseq_error_file = "";
	$del_file = "";
	$ins_file = "";

	if(defined $CMD{dist}) { $dist_file = $CMD{dist}; }
	if(defined $CMD{gff}) { $gff = $CMD{gff}; }
	if(defined $CMD{genome}) { $refseq_file = $CMD{genome}; }
	if(defined $CMD{referr}) { $refseq_error_file = $CMD{referr}; }
	if(defined $CMD{del}) { $del_file = $CMD{del}; }
	if(defined $CMD{ins}) { $ins_file = $CMD{ins}; }
}

