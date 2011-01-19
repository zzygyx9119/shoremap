#!/usr/bin/perl

# --------------------------------------------------------------------
#  ShoreMap extension to SHORE:
#  Identification of causal mutations using IlluminaGA2 sequencing data
# 
#  Written by Stephan Ossowski and Korbinian Schneeberger
#  --------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $referror = "";
my $marker;
my $marker2 = "";
my $region_chr;
my $consensus;
my $region_begin;
my $region_end;
my $parent1 = "";
my $parent2 = "";

my %CMD;
GetCom();



### Read in Ref Errors ------------------------------------------------------------------------
my %REFERROR = ();

if ($referror ne "") {
	open REFERRORFILE, $referror;
	while (my $line = <REFERRORFILE>) {
		my @a = split " ", $line;
		$REFERROR{$a[0]."#".$a[1]} = 1;
	}
	close REFERRORFILE;
	print STDERR "Finished reading referror file\n";
}


### Read in Marker ----------------------------------------------------------------------------
my %MARKER = ();

open MARKERFILE, $marker or die "Cannot open marker file\n";
while (my $line=<MARKERFILE>) {
	my @a = split " ", $line;
	if (not defined($REFERROR{$a[1]."#".$a[2]})) {
		$MARKER{$a[1]."#".$a[2]} = $a[4];
	}
}
close MARKERFILE or die "Marker file won't close\n";
print STDERR "Finished reading marker file\n";


### Read in Marker for second parent ----------------------------------------------------------
my %MARKER2 = ();

if ($marker2 ne "") {
open MARKERFILE, $marker2 or die "Cannot open marker file\n";
while (my $line=<MARKERFILE>) {
        my @a = split " ", $line;
        if (not defined($REFERROR{$a[1]."#".$a[2]})) {
                $MARKER2{$a[1]."#".$a[2]} = $a[4];
        }
}
close MARKERFILE or die "Marker file won't close\n";
print STDERR "Finished reading marker file 2\n";
}



### Get the counts at the marker positions ----------------------------------------------------
my @SNPCHR = ();
my @SNPPOS = ();
my @REFFREQ = ();
my @ALLELEFREQ = ();

#my $curr_chr = "NA";

open CONS, $consensus or die "Cannot open consensus file\n";
my $region_passed_flag = 0;

PARSING: while(my $line = <CONS>) {
	my @a = split " ", $line;
	my $chromosome = $a[0];
	my $position   = $a[1];	
	my $coverage   = $a[3];
	my $refbase    = uc($a[44]);

	if ($position%1000000 == 0) {
		print STDERR "Reading consensus at chromosome: ", $chromosome, " position: ", $position, "\n";
	}

	if ( (defined($MARKER{$chromosome."#".$position}) or defined($MARKER2{$chromosome."#".$position})) &&
	    ($refbase eq "A" || $refbase eq "C" || $refbase eq "G" || $refbase eq "T")) {

		if ($region_chr == $chromosome and $region_begin <= $position and $region_end >= $position) {
			$region_passed_flag = 1;

			my $snpbase = "";

        	        if (defined($MARKER{$chromosome."#".$position})) {
                	        $snpbase = uc($MARKER{$chromosome."#".$position});
	                }
        	        else {
                	        # First parent is col like
                        	$snpbase = $refbase;
	                }

			# refcount serves as count for the second parent
			if (defined($MARKER2{$chromosome."#".$position})) {
				$refbase = uc($MARKER2{$chromosome."#".$position});
			}
			
			my $snp_count = 0;
			my $ref_count = 0;
	
			if ($coverage != 0) {
	
				if ($snpbase eq "A") {
					$snp_count = $a[4];
				}
				if ($snpbase eq "C") {
					$snp_count = $a[5];
				}
				if ($snpbase eq "G") {
					$snp_count = $a[6];
				}
				if ($snpbase eq "T") {
					$snp_count = $a[7];
				}
	
				if ($refbase eq "A") {
					$ref_count = $a[4];
				}
				if ($refbase eq "C") {
					$ref_count = $a[5];
				}
				if ($refbase eq "G") {
					$ref_count = $a[6];
				}
				if ($refbase eq "T") {
					$ref_count = $a[7];
				}
			}
	
			push @SNPCHR, $chromosome;
			push @SNPPOS, $position;
			push @REFFREQ, $ref_count;
			push @ALLELEFREQ, $snp_count;
		
			#print $chromosome, "\t", $position, "\n";

		}
		else {
			if ($region_passed_flag == 1) {
				last PARSING;
			}
		}

	}

}
close CONS or die "Consensus file won't close\n";;
print STDERR "Finished reading consensus file\n";



### Plot allele frequencies --------------------------------------------------------------------

### Output the whole thing
my $outputfile = "SHOREmap.allele_freq.txt";
open OUT, "> ".$outputfile or die "Cannot open outputfile\n";
for (my $i = 0; $i < @SNPCHR; $i++) {
	print OUT $SNPCHR[$i], "\t", $SNPPOS[$i], "\t", $REFFREQ[$i], "\t", $ALLELEFREQ[$i], "\n";
}
close OUT;

### Call R for plotting
my $pdffile = "SHOREmap.allele_freq.$region_chr.$region_begin.$region_end.pdf";
my $cmd = "R --slave --vanilla --args $pdffile $outputfile $parent1 $parent2 < ".$FindBin::Bin."/SHOREmap_allele.R";
print STDERR $cmd, "\n";
system($cmd); 


sub min {
	my ($a, $b) = @_;
	return $a if $a < $b;
	return $b;
}

sub max {
	my ($a, $b) = @_;
	return $a if $a > $b;
	return $b;
}


### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

SHOREmap v2.0 (Marker file format changed!)

Mandatory:
--consensus    STRING   SHORE consensus file
--marker       STRING   Marker file
--chromosome   STRING   Chromosome ID
--begin        INT      Start position
--end          INT      Start position

Optional:
--marker2      STRING   Marker file (if none of the parents is the reference strain)
--referrors    STRING   Known reference sequence errors file

Plotting:
--parent1      STRING   Name of first parent (default \"\")
--parent2      STRING   Name of second parent (default \"\")

See documentation for file formats.

SHOREmap is written by Korbinian Schneeberger and Stephan Ossowski.
Max Planck Institute for Developmental Biology, Tübingen, 2009.

\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "parent1=s", "parent2=s", "chromosome=s", "begin=s", "end=s", "marker=s","marker2=s","consensus=s", "referrors=s");

        die("Please specify marker file\n") unless defined($CMD{marker});
        die("Please specify consensus file\n") unless defined($CMD{consensus});
        die("Please specify chromosome\n") unless defined($CMD{chromosome});
        die("Please specify begin\n") unless defined($CMD{begin});
        die("Please specify end\n") unless defined($CMD{end});

        $marker = $CMD{marker};
    	$region_chr = $CMD{chromosome};
    	$region_begin = $CMD{begin};
    	$region_end = $CMD{end};
        $consensus = $CMD{consensus};

	if (defined($CMD{parent1}) or defined($CMD{parent2})) {
		if (not defined($CMD{parent1}) or not defined($CMD{parent2})) {
			die("Please specify both parents\n");
		}
		$parent1 = $CMD{parent2};
		$parent2 = $CMD{parent1};
	}

	if (defined($CMD{referrors})) {
		$referror = $CMD{referrors};
	}

	if (defined($CMD{marker2})) {
                $marker2 = $CMD{marker2};
        }
}

