#! /usr/bin/perl
use strict;

my $usage= "\n$0 consensus_summary.txt markerfile\n\n" ;
my $file = shift or die $usage;
my $subset = shift or die $usage;

#print STDERR "\n... will chop \"Chr\" of the chromosome identifier of consensus_summary\n";

my %POS = ();

open FILE, $subset or die $usage;
while (my $line = <FILE>) {
        my @a = split " ", $line;
	my $id = $a[1]."#".$a[2];
	$POS{$id} = 1;
	#print STDERR $id, "\n";
}
close FILE;

open FILE, $file or die $usage;
while (my $line = <FILE>) {
        my @a = split " ", $line;
	#if (substr($a[0], 0, 3) eq "Chr") {
	#	$a[0] = substr($a[0], 3);
	#}
	my $id = $a[0]."#".$a[1];
	if (defined($POS{$id})) {
		print $line;
	}
}

