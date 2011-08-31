#!/usr/bin/perl

##############################################################################################
#  ShoreMap:
#  Identification of causal mutations through Next Gen Sequencing
# 
#  Written by 
#  Stephan Ossowski
#  Korbinian Schneeberger
#  Karl Nordstroem
##############################################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin;

### Command line parameters ##################################################################
my $expect;
my $confidence;
my $chrsizes;
my $out_folder;
my $marker;
my $marker_format;
my $consensus;
my $consensus_format;
my $window_size;
my $window_step;

my $peak_window_size;
my $peak_window_step;

my $filter_min_marker;
my $filter_min_coverage;

my $outlier_window_size;
my $outlier_pvalue;

my $misphenotyped;

my $reg_chromosome;
my $reg_begin;
my $reg_end;
my $reg_freq_min;
my $reg_freq_max;

my $referror;
my $background2;
my $verbose;

my $runid;
my $r_max;
my $plot_r;
my $boost_max;
my $plot_boost;

### Additional global variables ############################################################## 

my %CHR2SIZE = ();
my %REFERROR = ();
my %ALLELE1 = ();
my %ALLELE2 = ();
my %CMD;

init();

### Collect marker counts and sliding window values

my %CHR2POS2ALLELE1_COUNT = ();
my %CHR2POS2ALLELE2_COUNT = ();
my %CHR2POS2ERROR_COUNT = ();

read_allele_counts(); 


### Print allele counts

my $outputfile = $out_folder."/SHOREmap.winsize1.txt";
open OUT, "> ".$outputfile or die "Cannot open outputfile\n";

foreach my $chr (sort {$a cmp $b} keys %CHR2POS2ALLELE1_COUNT) {
	foreach my $pos (sort {$a <=> $b} keys %{$CHR2POS2ALLELE1_COUNT{$chr}}) {
		print OUT $chr, "\t", $pos, "\t", $CHR2POS2ALLELE1_COUNT{$chr}{$pos}, "\t", $CHR2POS2ALLELE2_COUNT{$chr}{$pos}, "\t", $CHR2POS2ERROR_COUNT{$chr}{$pos}, "\n";
	}
}



### Create sliding windows ###################################################################

# Set windowsizes

#my %WS2POS2FREQ = ();

#for (my $w = 0; $w < @window_sizes; $w++) {
# 	$window_size = $window_sizes[$w];
# 	$window_step = $window_size;
# 
# 	print STDERR "Parsing marker-wise AFE.\n" if $verbose == 1;
# 	my $outputfile = $out_folder."/SHOREmap.winsize1.txt";
#         open OUT, "> ".$outputfile or die "Cannot open outputfile\n";
# 
# 	foreach my $chr (sort {$a cmp $b} keys %CHR2POS2ALLELE1_COUNT) {
# 	
# 		my $last_report = 0;
# 
# 		# Window summary
# 		my $marker_sum = 0;
# 		my $allele1_sum = 0;
#                 my $allele2_sum = 0;
#                 my $error_sum = 0;
# 
# 		# Single marker data of current window
# 		my @marker_pos = ();
# 		my @allele1_counts = ();
# 		my @allele2_counts = ();
# 		my @error_counts = ();
# 
# 		foreach my $pos (sort {$a <=> $b} keys %{$CHR2POS2ALLELE1_COUNT{$chr}}) {
# 
# 			## remove markers that are not in the current window anymore
# 	        	while(defined($marker_pos[0]) and $marker_pos[0] < $pos - $window_size + 1) {
# 
#         	        	$marker_sum--;
# 				if ($marker_sum < 0) {
# 					$marker_sum = 0;
# 				}
# 		                $allele1_sum -= $allele1_counts[0];
# 				if ($allele1_sum < 0) {
# 					$allele1_sum = 0;
# 				}
#         		        $allele2_sum -= $allele2_counts[0];
# 				if ($allele2_sum < 0) {
# 					$allele2_sum = 0;
# 				}
# 				$error_sum -= $error_counts[0];
#                                 if ($error_sum < 0) {
#                                         $error_sum = 0;
#                                 }			
# 
# 				shift @marker_pos;
# 				shift @allele1_counts;
# 				shift @allele2_counts;
# 				shift @error_counts;
#         		}
# 
# 			# add new value
# 			$marker_sum++;
# 	        	$allele1_sum += $CHR2POS2ALLELE1_COUNT{$chr}{$pos};
# 	        	$allele2_sum += $CHR2POS2ALLELE2_COUNT{$chr}{$pos};
# 			$error_sum += $CHR2POS2ERROR_COUNT{$chr}{$pos};
# 
# 			# store new values
# 			push @marker_pos, $pos;
# 			push @allele1_counts, $CHR2POS2ALLELE1_COUNT{$chr}{$pos};
# 			push @allele2_counts, $CHR2POS2ALLELE2_COUNT{$chr}{$pos};
# 			push @error_counts, $CHR2POS2ERROR_COUNT{$chr}{$pos};
# 
# 			if ($marker_sum!= 0) {
# 
# 				# Calc where to report current window
# 				my $report_pos = int($pos- (($pos-$marker_pos[0]) / 2));
# 	
# 				# If current window distant enough to the last report (step size)
# 				# and if there are enough marker in the current window
# 
# 				if (($report_pos > $last_report + $window_step - 1) and
#                 	    	    ($filter_min_marker <= $marker_sum or $window_size == 1)) 
# 				{
# 
# 					#my $boost = get_boost($allele1_sum, $allele2_sum);
# 					print OUT $chr, "\t", $report_pos, "\t", $allele1_sum, "\t", $allele2_sum, "\t", $error_sum, "\n"; #, $boost, "\n";
# 
# 					$last_report = $report_pos;
# 
# 				}
#                 	}
# 		}
#         }
# 
# 	close OUT;
#}


my $pdffile = "$out_folder/SHOREmap.pdf";
my $cmd = "R --slave --vanilla --args $expect $chrsizes $pdffile $out_folder/SHOREmap.zoom_region.txt $window_size $window_step $FindBin::Bin $out_folder $outlier_window_size $outlier_pvalue $confidence $misphenotyped $filter_min_marker $filter_min_coverage $r_max $plot_r $boost_max $plot_boost $peak_window_size $peak_window_step $runid < ".$FindBin::Bin."/SHOREmap_plot.R"; # 2> /dev/null";
print STDERR $cmd, "\n" if $verbose == 1;
$cmd .= " 2> /dev/null" if $verbose == 0;
system($cmd);





##############################################################################################
# Subroutines

# Calc boost value
sub get_boost {
	my ($allele1_sum, $allele2_sum) = @_;

	my $boost = 0;
	my $sum = $allele1_sum + $allele2_sum;

	if ($sum == 0) {
		return $boost;#
	}

	my $obs = $allele1_sum/$sum;
                
	if ($obs != 0) {
		my $div = 1 - ($expect / $obs);
                if ($div == 0) {
	                $boost = 1000;
                }
                else {
                	$boost = min(1000, abs(1 / $div));
                }
                
	}

	return $boost;
}



# Read in marker alleles
sub read_allele_counts {

	open FILE, $consensus or die "Cannot open consensus file\n";
	while(my $line = <FILE>) {
	        my @a = split " ", $line;

        	my $chromosome;
	        my $position; 
		if ($consensus_format eq "shore") {
			$chromosome = $a[0];
			$position = $a[1];
		}
		elsif ($consensus_format eq "tab") {
			$chromosome = $a[0]; # DUMMY
			$position = $a[1];
		}

		if (substr($chromosome, 0, 3) eq "Chr") {
			$chromosome =~ s/Chr//g;
		}

		my $id = $chromosome."#".$position;  

		if (defined($ALLELE1{$id})) {

                        my $allele1 = uc($ALLELE1{$id});
                        my $allele2 = uc($ALLELE2{$id});

			my $count_allele1 = 0;
			my $count_allele2 = 0;
			my $count_error   = 0;

			my $coverage = 0;
			my $count_lower_allele = 0;
			my $count_higher_allele = 0;

			if ($consensus_format eq "shore") {

				$count_allele1 = $a[4] if ($allele1 eq "A"); 
				$count_allele1 = $a[5] if ($allele1 eq "C"); 
				$count_allele1 = $a[6] if ($allele1 eq "G"); 
				$count_allele1 = $a[7] if ($allele1 eq "T"); 

				$count_allele2 = $a[4] if ($allele2 eq "A");
                                $count_allele2 = $a[5] if ($allele2 eq "C");
                                $count_allele2 = $a[6] if ($allele2 eq "G");
                                $count_allele2 = $a[7] if ($allele2 eq "T");

                                $count_error += $a[4] if ($allele1 ne "A" and $allele2 ne "A");
                                $count_error += $a[5] if ($allele1 ne "C" and $allele2 ne "C");
                                $count_error += $a[6] if ($allele1 ne "G" and $allele2 ne "G");
                                $count_error += $a[7] if ($allele1 ne "T" and $allele2 ne "T");

				$coverage = $count_allele1 + $count_allele2 + $count_error;
				$count_lower_allele = $count_allele1 > $count_allele2 ? $count_allele2 : $count_allele1;
				$count_higher_allele = $count_allele1 > $count_allele2 ? $count_allele1 : $count_allele2;

			}
			elsif ($consensus_format eq "tab") {
                	        $chromosome = $a[0];
                        	$position = $a[1];
                        	$count_allele1 = $a[2];
                        	$count_allele2 = $a[3];
                        	$count_error = $a[4];
				
			}

			#if (	$coverage >= $filter_min_coverage and
			#	($filter_errors == 0 or $count_error < $count_lower_allele or $count_error < 0.2 * $count_higher_allele)
			#) {
	        	        $CHR2POS2ALLELE1_COUNT{$chromosome}{$position} = $count_allele1;
        	        	$CHR2POS2ALLELE2_COUNT{$chromosome}{$position} = $count_allele2;
	        	        $CHR2POS2ERROR_COUNT{$chromosome}{$position} = $count_error;
			#}
        	}
	}

}

sub read_referror {

        open FILE, $referror;
       	while (my $line = <FILE>) {
               	my @a = split " ", $line;
		my $chr = $a[0];
		my $pos = $a[1];

		# Check consistency
                if (not defined($CHR2SIZE{$chr})) {
			print STDERR "RefError file. Ignoring entry. \"$chr\" does not match entries in chromosome sizes file.\n"
		}
		else {

	                if ($pos =~ m/[^0-9.]/ ) { die("Referror position not numeric ($pos).\n"); }

			if ($CHR2SIZE{$chr} < $pos) { die("Referror position out of chromosome bounds ($chr:$pos).\n");}

                	$REFERROR{$a[0]."#".$a[1]} = 1;
		}
       	}
        close FILE;
       	print STDERR "Finished reading referror file\n" if $verbose == 1;

	
}

sub read_marker {

	open FILE, $marker or die "Cannot open file. ($marker)\n";
	while (my $line=<FILE>) {
        	my @a = split " ", $line;
		my $chr;
		my $pos;
		my $allele1;
		my $allele2;
	        if ($marker_format eq "shore"){
			$chr = $a[1]; 
			$pos = $a[2];
			$allele1 = $a[3];
			$allele2 = $a[4];
        	}
        	elsif ($marker_format eq "vcf") {
			$chr = $a[0];
                        $pos = $a[1]; ## TODO check indels are those correctly doing what I think?
                        $allele1 = $a[3];
                        $allele2 = $a[4];
        	}

		# Check for consistency
		if (not defined($CHR2SIZE{$chr})) {
			#die("Marker reading: Wrong chromosome annoation in marker file: \"$chr\" Marker format is $marker_format.\nDoes not match entries in chromosome sizes file.\n");
			print STDERR "Marker file. Ignoring entry. \"$chr\" does not match entries in chromosome sizes file.\n"
		}
		else {

			if ($pos =~ m/[^0-9.]/ ) { die("Marker position not numeric ($pos).\n"); }

			if ($CHR2SIZE{$chr} < $pos) { die("Marker position out of chromosome bounds ($chr:$pos).\n");}

			# Store
			if (not defined($REFERROR{$chr."#".$pos})) {
				if ($background2 == 0) {
					$ALLELE1{$chr."#".$pos} = $allele1;
        	                	$ALLELE2{$chr."#".$pos} = $allele2;
				}
				else {
					$ALLELE1{$chr."#".$pos} = $allele2;
                        	        $ALLELE2{$chr."#".$pos} = $allele1;
				}
			}
                }
	
	}
	close FILE or die "Marker file won't close. ($marker)\n";

	print STDERR "Finished reading marker file\n" if $verbose == 1;

}

# Read in chromosome sizes
sub read_chromosomes {

	open FILE, "$chrsizes" or die "Cannot open $chrsizes\n";
	while (<FILE>) {
		chomp();
		my @a = split " ";
		if ($a[1] =~ m/[^0-9.]/ ) { die("Chromosome size not numeric ($_).\n");}
		$CHR2SIZE{$a[0]} = $a[1];
	}

}

sub write_log {

	open FILE, ">$out_folder/SHOREmap.log";

	print FILE "# perl $0 --target $expect --chrsizes $chrsizes --folder $out_folder --marker $marker --marker-format $marker_format --consen $consensus --consen-format $consensus_format --window-size $window_size --window-step $window_step --mis-phenotyped $misphenotyped --min-marker $filter_min_marker --min-coverage $filter_min_coverage --outlier-window-size $outlier_window_size --outlier-pvalue $outlier_pvalue ";

	if ($confidence == 2) {
		print FILE "-no-interval ";
	}
	else {
		print FILE "--conf $confidence ";
	}

	if ($reg_chromosome ne "") {
		print FILE "--chromosome $reg_chromosome --begin $reg_begin --end $reg_end --minfreq $reg_freq_min --maxfreq $reg_freq_max ";
	}
	if ($referror ne "") {
		print FILE "--referrors $referror ";
	}
	if ($background2 != 0)  {
		print FILE "-background2 ";
	}
	if ($verbose != 0) {
		print FILE "-verbose ";
	}
	if ($plot_boost != 0) {
		print FILE "-plot_boost ";
	}
	if ($plot_r != 0) {
                print FILE "-plot_r ";
        }

	print FILE "\n";
		
	close FILE;
}

sub write_zoom_region {
	if ($reg_chromosome ne "") {
	        open FILE, ">$out_folder/SHOREmap.zoom_region.txt" or die "Cannot open file $out_folder/SHOREmap.zoom_region.txt\n";
        	print FILE $reg_chromosome, "\t", $reg_begin, "\t", $reg_end, "\t", $reg_freq_min, "\t", $reg_freq_max, "\n";
	        close FILE;
	}
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


### Read command line parameters
sub init {

	my @usage = ("$0

SHOREmap 2.0

Mandatory:
--chrsizes              STRING   Chromosome sizes file
--folder                STRING   Output folder
--marker                STRING   Marker file
--marker-format         STRING   Marker file format, \"shore\" or 
                                 \"vcf\" (default: \"shore\")
--consen                STRING   Consensus file 
--consen-format         STRING   Consensus file format, \"shore\" 
                                 or \"tab\" (default: \"shore\")

Confidence interval:
--target                DOUBLE   Target allele frequency 
                                 (default: 1.0)
--conf                  DOUBLE   Confidence level
                                 (default: 0.99)

--peak-window-size      INT      (default: 50000)
--peak-window-step      INT      (default: 10000)
                                 Used for initial seed finding

Visulization:
--window-size           INT      (default: 50000)
--window-step           INT      (default: 10000)
                                 Used for smoothed visulization
                                 and \"boost\"-value calculation

Filter:
--min-marker            INT      Filter windows with low numbers 
                                 of markers (default: 10)
--min-coverage          INT      Filter single marker with low 
                                 average coverage (default: 0)
--outlier-window-size   INT      Window size to assess local
                                 allele frequency used for 
                                 outlier removal 
                                 (default: 200000)
--outlier-pvalue        DOUBLE   p-value used for outlier
                                 removal (default: 0.05)

Phenotyping:
--mis-phenotyped        DOUBLE   Degree of putatively
                                 mis-scored plants
                                 (default: 0.00)

Zooming:
--chromosome            INT      Zoom to chromosome ..
--begin                 INT      .. from here ..
--end                   INT      .. to here with a ..
--minfreq               INT      .. minimal to ..
--maxfreq               INT      .. maximal frequency.


Optional:
--referrors             STRING   Reference errors file
-plot-r                          Plot frequency calculation (\"r\")
-plot-boost                      Plot frequency calculation 
                                 (\"boost\")
-background2                     Mutation is in second parent
-no-interval                     Switch off confidence interval
                                 calculation
-verbose                         Be talkative

See documentation for file formats.
");

## HIDDEN FLAG
# boost max, r max, runid

	$expect = 1.0;
	$confidence = 0.99;
	$chrsizes = "";
	$out_folder = "";
	$marker = "";
	$marker_format = "shore";
	$consensus = "";
	$consensus_format = "shore";
	$window_size = 50000;
	$window_step = 10000;

	$peak_window_size = 50000;
	$peak_window_step = 10000;

	$filter_min_marker = 0;
	$filter_min_coverage = 0;

	$outlier_window_size = 200000;
	$outlier_pvalue = 0.05;
	
	$misphenotyped = 0.00;

	$reg_chromosome = "";
	$reg_begin = -1;
	$reg_end = -1;
	$reg_freq_min = -1;
	$reg_freq_max = -1;

	$referror = "";
	$background2 = 0;
	$verbose = 0;	

	$boost_max = 10000;
	$plot_boost = 0;
	$r_max = 10000;
	$plot_r = 0;

	$runid = 1;

	if (@ARGV+0 == 0) {
		print @usage, "\n";
		exit(0);
	}

        GetOptions(\%CMD, "target=f", "conf=f", "chrsizes=s", "folder=s", "marker=s", "marker-format=s", "consen=s", "consen-format=s", "window-size=i", "window-step=i", "min-marker=i", "min-coverage=i", "outlier-window-size=i", "outlier-pvalue=f", "mis-phenotyped=f", "chromosome=i", "begin=i", "end=i", "minfreq=f", "maxfreq=f", "verbose", "background2", "referrors=s", "no-interval", "runid=i", "boost-max=i", "plot-boost", "r-max=i", "plot-r"); 


        die("Please specify chromosome sizes file\n") unless defined($CMD{chrsizes});
        die("Please specify output folder\n") unless defined($CMD{folder});
        die("Please specify marker file\n") unless defined($CMD{marker});
        die("Please specify consensus file\n") unless defined($CMD{consen});

	####################################################################################
	## Store values and check for consistency

        if (defined($CMD{background2})) {
                $background2 = 1;
        }

	if (defined($CMD{target})) {
		$expect = $CMD{target};
	}
	if (!($expect >= 0 and $expect <= 1.0)) {
		die("target must be between 0 and 1.\n");
	}

	if (defined($CMD{conf})) {
                $confidence = $CMD{"conf"};
		if ($confidence <= 0 or $confidence > 1) {
			die ("Confidence level needs to be larger 0 and smaller 1.\n");
		}
        }

	if (defined($CMD{"no-interval"})) {
		$confidence = 2;
	} 

	$chrsizes = $CMD{chrsizes};
	if (!-e $chrsizes) {
                die("File not found: $chrsizes\n");
        }
	read_chromosomes();

        if (defined($CMD{referrors})) {
                $referror = $CMD{referrors};
                if (!-e $referror) {
                        die("File not found: $referror\n");
                }
		read_referror();
        }
	
	$out_folder = $CMD{folder};
	if (-e $out_folder) {
		die("Outfolder exists already.\n");
	}
	mkdir($out_folder);

        if (defined($CMD{"marker-format"})) {
                $marker_format = $CMD{"marker-format"};
                if ($marker_format ne "shore" and $marker_format ne "vcf") {
                        die("marker-format has to be \"shore\" or \"vcf\"\n");
                }
        }

        $marker = $CMD{marker};
	if (!-e $marker) {
                die("File not found: $marker\n");
        }
	read_marker();

        $consensus = $CMD{consen};
	if (!-e $consensus) {
                die("File not found: $consensus\n");
        }

	if (defined($CMD{"consen-format"})) {
                $consensus_format = $CMD{"consen-format"};
                if ($marker_format ne "shore" and $marker_format ne "tab") {
                        die("marker-format has to be \"shore\" or \"tab\"\n");
                }
        }

	if (defined($CMD{"window-size"})) {
                $window_size = $CMD{"window-size"};
	}
	if ($window_size =~ m/[^0-9.]/ ) { die("Window size is not numeric ($window_size).\n");}
        if ($window_size <= 1) { die("Window size must be larger than 1 ($window_size).\n");}

	if (defined($CMD{"window-step"})) {
		$window_step = $CMD{"window-step"};
		if ($window_step =~ m/[^0-9.]/ ) { die("Window step size is not numeric ($window_step).\n");}
		if ($window_step < 1) { die("Window step size smaller than 1 not valid ($window_step).\n");}
	}

	if (defined($CMD{"min-marker"})) {
                $filter_min_marker = $CMD{"min-marker"};
                if ($filter_min_marker =~ m/[^0-9.]/ ) { die("Minimal number of markers is not numeric ($filter_min_marker).\n");}
        }

	if (defined($CMD{"min-coverage"})) {
                $filter_min_coverage = $CMD{"min-coverage"};
                if ($filter_min_coverage =~ m/[^0-9.]/ ) { die("Minimal coverage is not numeric ($filter_min_coverage).\n");}
        }

	if (defined($CMD{"outlier-window-size"})) {
                $outlier_window_size = $CMD{"outlier-window-size"};
                if ($outlier_window_size =~ m/[^0-9.]/ ) { die("Outlier window size is not numeric ($outlier_window_size).\n");}
        }

	if (defined($CMD{"outlier-pvalue"})) {
                $outlier_pvalue = $CMD{"outlier-pvalue"};
                if ($outlier_pvalue =~ m/[^0-9.]/ ) { die("Outlier p-value is not numeric ($outlier_pvalue).\n");}
        }

	if (defined($CMD{"mis-phenotyped"})) {
		$misphenotyped = $CMD{"mis-phenotyped"};
		if ($misphenotyped =~ m/[^0-9.]/ ) { die("Misphenotyping is not numeric ($misphenotyped).\n");}
		if ($misphenotyped < 0 or $misphenotyped > 1) { die("Misphenotyping needs to be between 0 and 1.\n");}
	}


	if (defined($CMD{chromosome}) or defined($CMD{begin}) or defined($CMD{end}) or defined($CMD{minfreq}) or defined($CMD{maxfreq})) {
                if (not(defined($CMD{chromosome}) and defined($CMD{begin}) and defined($CMD{end}) and defined($CMD{minfreq}) and defined($CMD{maxfreq}))) {
                        die("You need to define chromosome, begin, end, minfreq and maxfreq to zoom in.\n");
                }
                $reg_chromosome = $CMD{chromosome};
		if (not defined($CHR2SIZE{$reg_chromosome})) {
                        die("Wrong chromosome annoation for zoom.\n Does not match entries in chromosome sizes file.\n");
                }
		
                $reg_begin = $CMD{begin};
		if ($reg_begin =~ m/[^0-9.]/ ) { die("Begin of zoom region not numeric ($reg_begin).\n");}
                if ($reg_begin < 1 or $CHR2SIZE{$reg_chromosome} < $reg_begin) { die("Zoom begin out of chromosome bounds ($reg_chromosome:$reg_begin).\n");}

                $reg_end = $CMD{end};
		if ($reg_end =~ m/[^0-9.]/ ) { die("End of zoom region not numeric ($reg_end).\n");}
                if ($reg_end < 1 or $CHR2SIZE{$reg_chromosome} < $reg_end) { die("Zoom end out of chromosome bounds ($reg_chromosome:$reg_end).\n");}
                if ($reg_end < $reg_begin) { die("Zoom begin larger than zoom end ($reg_begin:$reg_end).\n");}

                $reg_freq_min = $CMD{minfreq};
		if ($reg_freq_min =~ m/[^0-9.]/ ) { die("Min. freq. of zoom region not numeric ($reg_freq_min).\n");}
                if ($reg_freq_min < 0 or $reg_freq_min > 1) { die("Zoom min. freq. out of bounds ($reg_freq_min).\n");}

                $reg_freq_max = $CMD{maxfreq};
		if ($reg_freq_max =~ m/[^0-9.]/ ) { die("Max. freq. of zoom region not numeric ($reg_freq_max).\n");}
                if ($reg_freq_max < 0 or $reg_freq_max > 1) { die("Zoom max. freq. out of bounds ($reg_freq_max).\n");}
		if ($reg_freq_min > $reg_freq_max) { die("Zoom min. freq. larger than zoom max. freq. ($reg_freq_min:$reg_freq_max).\n");}

        }

	if (defined($CMD{verbose})) {
                $verbose = 1;
        }

        if (defined($CMD{"runid"})) {
                $runid = $CMD{"runid"};
        }

	if (defined($CMD{"boost-max"})) {
                $boost_max = $CMD{"boost-max"};
        }

	if (defined($CMD{"plot-boost"})) {
                $plot_boost = 1;
		if (defined($CMD{"plot-r"})) {
			die("Select either \"r\" or \"boost\".\n");
		}
        }

	if (defined($CMD{"r-max"})) {
                $r_max = $CMD{"r-max"};
        }

	if (defined($CMD{"plot-r"})) {
                $plot_r = 1;
		if (defined($CMD{"plot-boost"})) {
                        die("Select either \"r\" or \"boost\".\n");
                }
        }

	####################################################################################

	mkdir($out_folder);
	write_log();
	write_zoom_region();
	

}

