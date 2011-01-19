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
my $method = "boost";
my $consensus;
my $chrsizes;
my $window_step;
my $window_size_string;
my $expect = 1.0;
my $boost_single = 1;
my $min_allele_window = 100;
my $parent1 = "";
my $parent2 = "";
my $no_alleles = 0;
my $old_format = 0;
my $verbose = 0;

my $reg_chromosome = -1;
my $reg_begin = -1;
my $reg_end = -1;
my $reg_freq_min = -1;
my $reg_freq_max = -1;


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
	print STDERR "Finished reading referror file\n" if $verbose == 1;
}


### Read in Marker ----------------------------------------------------------------------------
my %MARKER = ();
my %REF = ();

open MARKERFILE, $marker or die "Cannot open marker file\n";
while (my $line=<MARKERFILE>) {
	my @a = split " ", $line;
	if ($old_format == 0){
		if (not defined($REFERROR{$a[1]."#".$a[2]})) {
			$MARKER{$a[1]."#".$a[2]} = $a[4];
			$REF{$a[1]."#".$a[2]} = $a[3];
		}
	}
	else {
	        if (not defined($REFERROR{$a[1]."#".$a[2]})) {
			$MARKER{$a[0]."#".$a[1]} = $a[3];
			$REF{$a[0]."#".$a[1]} = $a[2];
		}
	}
}
close MARKERFILE or die "Marker file won't close\n";
print STDERR "Finished reading marker file\n" if $verbose == 1;



### Get the counts at the marker positions ----------------------------------------------------
my @SNPCHR = ();
my @SNPPOS = ();
my @REFFREQ = ();
my @ALLELEFREQ = ();

#my $curr_chr = "NA";

open CONS, $consensus or die "Cannot open consensus file\n";
while(my $line = <CONS>) {
	my @a = split " ", $line;
	my $chromosome = $a[0];
	my $position   = $a[1];	
	my $coverage   = $a[3];

	if ($position%1000000 == 0) {
		print STDERR "Reading consensus at chromosome: ", $chromosome, " position: ", $position, "\n" if $verbose == 1;
	}


	if (defined($MARKER{$chromosome."#".$position})) {

		my $snpbase = uc($MARKER{$chromosome."#".$position});
		my $refbase = uc($REF{$chromosome."#".$position});

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
}
close CONS or die "Consensus file won't close\n";
print STDERR "Finished reading consensus file\n" if $verbose == 1;



### Create sliding windows --------------------------------------------------------------------

# Set windowsizes
my @window_sizes = split ",", $window_size_string;

for (my $w = 0; $w < @window_sizes; $w++) {
	my $window_size = $window_sizes[$w];

	print STDERR "Analysing window size: ", $window_size, "\n" if $verbose == 1;

	my $output_string = "";

	my $last_pos = 0;
	my $marker_count = 0;
	my $ref_count = 0;
	my $allele_count = 0;

	my $j = 0;
	for (my $i = 0; $i < @SNPCHR; $i++) {

#print STDERR $i, "\t", $j, "\t\t";

        	my $chr = $SNPCHR[$i];
	        my $pos = $SNPPOS[$i];
        	my $ref = $REFFREQ[$i];
	        my $mut = $ALLELEFREQ[$i];

		# delete old values
		# todo finish last window?
		if ($SNPCHR[$j] ne $chr) {
			$j = $i;
	                $marker_count = 0;
	                $ref_count = 0;
        	        $allele_count = 0;
			$last_pos = 0; 
       		}

        	while($SNPPOS[$j] < $pos - $window_size + 1) {
                	$marker_count--;
			if ($marker_count < 0) {
				$marker_count = 0;
			}
	                $ref_count -= $REFFREQ[$j];
			if ($ref_count < 0) {
				$ref_count = 0;
			}
        	        $allele_count -= $ALLELEFREQ[$j];
			if ($allele_count < 0) {
				$allele_count = 0;
			}
                	$j++;
        	}

		# add new value
		$marker_count++;
	        $ref_count += $REFFREQ[$i];
        	$allele_count += $ALLELEFREQ[$i];

		

		# check if the regions is in the selected interval - if zooming was switched on
		if ($reg_chromosome == -1 or ($reg_chromosome == $chr and $reg_begin <= int($pos-(($pos-$SNPPOS[$j])/2)) and $reg_end >= int($pos-(($pos-$SNPPOS[$j])/2)))) {

		# Report if more than
		if ($pos > $last_pos + $window_step - 1) {
                	if ($marker_count!= 0) {

				$output_string .= $chr."\t";
				$output_string .= int($pos-(($pos-$SNPPOS[$j])/2))."\t";
				$output_string .= $ref_count."\t".$allele_count."\t";

				$last_pos = $pos;
#print STDERR "printed\n";

				my $val = 0;
				my $power = 1;
				my $pi = 3.141593;
				my $max_val = 1000;


				###############################################################################
				# allele:
				###############################################################################

				if ($method eq "afreq") {

					if ($ref_count + $allele_count != 0) {
						$val = $ref_count / ($ref_count + $allele_count);
					}
					else {
						$val = 0.5;
					}

                                }


				###############################################################################
				# abs:
				###############################################################################

				if ($method eq "abs") {

					$val = $allele_count."\t".$ref_count;

                                }

				###############################################################################
			        # r-value:
			        ###############################################################################

				if ($method eq "r") {

					if ($ref_count > $allele_count) {
						if ($allele_count != 0) {
							$val = $ref_count / $allele_count;
						}
						else {
							$val = $ref_count;	
						}
					}
					else {
						if ($ref_count == $allele_count) {
							$val = 0;
						}
						else {
							if ($ref_count != 0) {
								$val = (-1) * ($allele_count / $ref_count);
							}
							else {
								$val = (-1) * $allele_count;
							}
						}
					}

				}
			       

				###############################################################################
			        # Cos:
			        ###############################################################################

				if ($method eq "cos") {

					my $total = $allele_count + $ref_count;
					my $larger = max($allele_count, $ref_count);
					if ($total == 0) {
						$val = 0;
					}
					else {
						my $perc = $larger / $total;
						$val = (cos($perc*$pi*2)/2)+0.5;
					}

					if ($ref_count < $allele_count) {
						$val = $val * (-1);
					}

				}

				###############################################################################
				# Log:
				###############################################################################
		
				if ($method eq "log") {
		
					my $total = $allele_count + $ref_count;
					if ($total == 0) {
						$val = 0;
					}
					else {
						my $val1 = 0;
						if ($allele_count < $ref_count) {
							$val1 = $ref_count/$total;
						}
						else {
                                	                $val1 = $allele_count/$total;
						}
					
						if ($val1 >= 0.45 and $val1 <= 0.56) {
							$val = 0;
						}
						else {
							$val = log($val1-0.55)+4.60517;
						}
					}

					if ($ref_count < $allele_count) {
						$val = $val * (-1);
					}

				}

				###############################################################################
				# Hyperbel:
				###############################################################################

				if ($method eq "hyperbola") {

					my $total = $allele_count + $ref_count;
        	                        if ($total == 0) {
                	                        $val = 0;
                        	        }
                                	else {
                                        	my $val1 = 0;
	                                        if ($allele_count < $ref_count) {
        	                                        $val1 = $ref_count/$total;
                	                        }
                        	                else {
                                	                $val1 = $allele_count/$total;
                                        	}

	                                        $val = $val1**6;
					}

					if ($ref_count < $allele_count) {
						$val = $val * (-1);
					}

                                }				


				###############################################################################
				# Boost:
				###############################################################################
				
				if ($method eq "boost") {

					my $total = $ref_count+$allele_count;
					if ($total < $min_allele_window) {
						$val = 0;
					}
					else {
						my $obs;
						if ($ref_count > $allele_count) {
							$obs = $ref_count/$total;
						}
						else {
							$obs = $allele_count/$total;
						}
	
						if ($obs == 0) {
							$val = 0;
						}
						elsif ($obs >= $expect and $boost_single == 0) {
							$val = $max_val;
						}
						else {
							my $div = 1 - (($expect) / $obs)**$power;
							if ($div == 0) {
								$val = $max_val;
							}
							else {
								$val = min($max_val, abs(1 / $div));
							}
						}
					}
	
					if ($ref_count < $allele_count and $no_alleles == 0) {
						$val = $val * (-1);
					}

				}

				###############################################################################


				$output_string .= "\t".$val."\n";

			}
                }}
        }


	### Output the whole thing
	my $outputfile = "SHOREmap.".$method.".winstep".$window_step.".winsize".$window_size.".txt";
	open OUT, "> ".$outputfile or die "Cannot open outputfile\n";
	print OUT $output_string;
	close OUT;


	my $pdffile = "SHOREmap.".$method.".winstep".$window_step.".winsize".$window_size.".pdf";
	my $cmd = "";
	if ($method eq "abs") {
                $cmd = "R --slave --vanilla --args $chrsizes $pdffile $outputfile $window_step $window_size $parent1 $parent2 < ".$FindBin::Bin."/SHOREmap_abs.R 2> /dev/null";
	}
	elsif ($method eq "afreq") {
		if ($reg_chromosome == -1) {
			$cmd = "R --slave --vanilla --args $chrsizes $pdffile $outputfile $window_step $window_size $parent1 $parent2 2> /dev/null < ".$FindBin::Bin."/SHOREmap_afreq.R";
		}
		else {
			#$cmd = "R --slave --vanilla --args $chrsizes $pdffile $outputfile $window_step $window_size $reg_chromosome $reg_begin $reg_end $reg_freq_min $reg_freq_max 2> /dev/null < ".$FindBin::Bin."/SHOREmap_afreq_reg.R";
			$cmd = "R --slave --vanilla --args $chrsizes $pdffile $outputfile $window_step $window_size $reg_chromosome $reg_begin $reg_end $reg_freq_min $reg_freq_max < ".$FindBin::Bin."/SHOREmap_afreq_reg.R";
		}
	}
	else {
		$cmd = "R --slave --vanilla --args $chrsizes $pdffile $outputfile $window_step $window_size $parent1 $parent2 < ".$FindBin::Bin."/SHOREmap_2.R 2> /dev/null";
	}
        print STDERR $cmd, "\n" if $verbose == 1;
        system($cmd);

}

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

SHOREmap v2.0 
(-Marker file format changed)
(-Two non-reference parental lines possible)
(-Allele frequency visualization)
(-Zoom level)

Mandatory:
--consensus    STRING   SHORE consensus file
--marker       STRING   Marker file (ref & parent or parent1 & parent2)
--chrsizes     STRING   Chromosome sizes file
--agg          STRING   Agglomeration: \"boost\", \"r\", \"afreq\" (default: boost)
--exp          DOUBLE   For \"boost\" only: Expected allele ratio for causal region (default 1.0)

Optional:
--referrors    STRING   Known reference sequence errors file
-verbose                Be talkative
-old-format             Old format of the variation files

Plotting:
--parent1      STRING   Name of first parent (default \"\")
--parent2      STRING   Name of second parent (default \"\")
--min          INT      Minimum number of alleles within window (100)
-no_alleles             Do not distinguish between parents (e.g. with boost and 50% allele freq)

--chromosome   INT	Zoom to one chromosome ...
--begin        INT      ... from here ...
--end          INT      ... to here with a ...
--minfreq      INT      ... minimal to ...
--maxfreq      INT      ... maximal allele frequency.

--windowstep   INT      Number of bp between datapoints. Default is 10000.
--windowsize   STRING   Comma separated windowsizes. 
               Default is: 50000,100000,150000,200000,250000,300000,350000,400000,450000,500000

See documentation for file formats.

SHOREmap was developed by 
Korbinian Schneeberger and Stephan Ossowski.

\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "chromosome=s", "begin=s", "end=s", "minfreq=s", "maxfreq=s", "verbose", "parent1=s", "parent2=s", "agg=s", "min=f", "exp=f", "marker=s","consensus=s", "chrsizes=s", "windowsize=s", "windowstep=s", "referrors=s", "no_alleles","old-format");

        die("Please specify marker file\n") unless defined($CMD{marker});
        die("Please specify consensus file\n") unless defined($CMD{consensus});
        die("Please specify chromosome sizes file\n") unless defined($CMD{chrsizes});

	if (defined($CMD{agg})) {
		$method = $CMD{agg};
	}
        $marker = $CMD{marker};
        $consensus = $CMD{consensus};
	$chrsizes = $CMD{chrsizes};

	if ($method ne "hyperbola" and $method ne "log" and $method ne "cos" and $method ne "boost" and $method ne "r" and $method ne "afreq" and $method ne "abs") {
		die("Agglormeration function must be one of \"boost\", \"hyperbola\", \"log\", \"cos\", \"r\", \"afreq\", \"abs\".\n");
	}

	if (defined($CMD{parent1}) or defined($CMD{parent2})) {
		if (not defined($CMD{parent1}) or not defined($CMD{parent2})) {
			die("Please specify both parents\n");
		}
		$parent1 = $CMD{parent2};
		$parent2 = $CMD{parent1};
	}

	if (defined($CMD{verbose})) {
                $verbose = 1;
        }

	if (defined($CMD{no_alleles})) {
                $no_alleles = 1;
        }

	if (defined($CMD{old_format})) {
		$old_format = 1;
	}

	if (defined($CMD{min})) {
                $min_allele_window = $CMD{min};
        }

	if (defined($CMD{exp})) {
                $expect = $CMD{exp};
        }

	if (defined($CMD{referrors})) {
		$referror = $CMD{referrors};
	}

	if (defined($CMD{windowstep})) {
		$window_step = $CMD{windowstep};
	}
	else {
		$window_step = 10000;
	}

	if (defined($CMD{windowsize})) {
                $window_size_string = $CMD{windowsize};
        }
        else {
                $window_size_string = "50000,100000,150000,200000,250000,300000,350000,400000,450000,500000";
        }

	if (defined($CMD{chromosome}) or defined($CMD{begin}) or defined($CMD{end})) {
		if (not(defined($CMD{chromosome}) and defined($CMD{begin}) and defined($CMD{end}))) {
			die("You need to define chromosome, begin and end to zoom in.\n");
		}
		$reg_chromosome = $CMD{chromosome};
		$reg_begin = $CMD{begin};
		$reg_end = $CMD{end};
		$reg_freq_min = $CMD{minfreq};
		$reg_freq_max = $CMD{maxfreq};
	}

}

