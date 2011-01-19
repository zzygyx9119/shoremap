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
my $variant;
my $window_step;
my $window_size;
my $type = "freq";
my %CMD;
GetCom();



### Get the  frequencies at called variant position --------------------------------------------
my $sum = 0;
my $count = 1;
my $chr = "NA";
my $win_start = 1;
my $win_end = $window_size;

open VAR, $variant or die "Cannot open variant file\n";

while(my $line = <VAR>) {
	my @a = split " ", $line;

	### New chromosome reached
	if( $chr ne $a[1] ) {
		if( $chr ne "NA" ) {
			my $avg_freq = $sum/$count;
			if( $type eq "freq") {
				print "$chr\t$win_start\t$win_end\t$avg_freq\n";
			}
			else {
				print "$chr\t$win_start\t$win_end\t$sum\n";
			}
		}
		$chr = $a[1];
		$win_start = 1;
		$win_end = $window_size;
	}

	### Sliding window end reached
	while( $a[2] > $win_end) {
		my $avg_freq = $sum/$count;
		if( $type eq "freq") {
			print "$chr\t$win_start\t$win_end\t$avg_freq\n";
		}
		else {
			print "$chr\t$win_start\t$win_end\t$sum\n";
		}

		$sum = 0;
		$count = 1;
		$win_start += $window_step;
		$win_end += $window_step;
	}
	
	$sum += $a[7];
	$count++;
}

close VAR or die "Variation file won't close\n";;
print STDERR "Finished reading variation file\n";


### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

SHOREmap v2.0 (Marker file format changed!)

Mandatory:
--variant      STRING   SHORE consensus file
--winsize      INT      Sliding window size for plot
--winstep      INT      Sliding window step size for plot
--type         STRING   <freq | sum>

See documentation for file formats.

SHOREmap is written by Korbinian Schneeberger and Stephan Ossowski.
Max Planck Institute for Developmental Biology, Tübingen, 2009.

\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "variant=s","winsize=s", "winstep=s", "type=s");

        die("Please specify variant file\n") unless defined($CMD{variant});
        die("Please specify winsize\n") unless defined($CMD{winsize});
        die("Please specify winstep\n") unless defined($CMD{winstep});

        $variant     = $CMD{variant};
    	$window_size = $CMD{winsize};
    	$window_step = $CMD{winstep};

	if(defined($CMD{type})) {
		$type = $CMD{type};
	}
}

