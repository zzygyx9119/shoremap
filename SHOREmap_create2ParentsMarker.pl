#! /usr/bin/perl
use strict;
use warnings;


my $usage = "\n$0 variants1 reference1 variants2 reference2\n\n";

my $var1 = shift or die $usage;
my $ref1 = shift or die $usage;
my $var2 = shift or die $usage;
my $ref2 = shift or die $usage;

my %REF1 = ();
my %REF2 = ();
my %VAR1 = ();
my %VAR2 = ();

my $MIN_QUAL = 25;

open FILE, $var1 or die $usage;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($MIN_QUAL <= $a[5]) {
		my $id = $a[1] * 100000000 + $a[2]; 
		$VAR1{$id} = $line;
	}
}
close FILE;

open FILE, $var2 or die $usage;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($MIN_QUAL <= $a[5]) {
		my $id = $a[1] * 100000000 + $a[2]; 
		$VAR2{$id} = $line;
	}
}
close FILE;

open FILE, $ref1 or die $usage;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($MIN_QUAL <= $a[5]) {
		my $id = $a[1] * 100000000 + $a[2]; 
		#if (defined($VAR2{$id})) {
			$REF1{$id} = $line;
		#}
	}
}
close FILE;

open FILE, $ref2 or die $usage;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($MIN_QUAL <= $a[5]) {
		my $id = $a[1] * 100000000 + $a[2]; 
		#if (defined($VAR1{$id})) {
			$REF2{$id} = $line;
		#}
	}
}
close FILE;

open COMBINED_195, ">markers.combined.c195.txt.not_sorted";
open COMBINED_190, ">markers.combined.c190.txt.not_sorted";
open COMBINED_175, ">markers.combined.c175.txt.not_sorted";


foreach my $m (sort {$a <=> $b} keys %VAR1) {
	my $cref = 0;
	for (my $i = $m-100; $i<=$m+100; $i++) {
		if (defined($REF1{$i})) {
			$cref++;
		}
	}
	# change alleles
	my @a = split " ", $VAR1{$m};
	my $new_line = $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[4]."\t".$a[3]."\t".$a[5]."\t".$a[6]."\t".$a[7]."\t".$a[8]."\n"; 
	if ($cref >= 195) {
                if (defined($REF2{$m}) and not defined($VAR2{$m})) {
                        print COMBINED_195 $new_line;
                }
        }
	if ($cref >= 190) {
		if (defined($REF2{$m}) and not defined($VAR2{$m})) {
			print COMBINED_190 $new_line;
		}
	}
	if ($cref >= 175) {
                if (defined($REF2{$m}) and not defined($VAR2{$m})) {
                        print COMBINED_175 $new_line;
                }
        }	
}

foreach my $m (sort {$a <=> $b} keys %VAR2) {
	my $cref = 0;
        for (my $i = $m-100; $i<=$m+100; $i++) {
                if (defined($REF2{$i})) {
                        $cref++;
                }
        }
	if ($cref >= 195) {
                if (defined($REF1{$m}) and not defined($VAR1{$m})) {
                        print COMBINED_195 $VAR2{$m};
                }
        }
	if ($cref >= 190) {
		if (defined($REF1{$m}) and not defined($VAR1{$m})) {
			print COMBINED_190 $VAR2{$m};
		}
	}
	if ($cref >= 175) {
                if (defined($REF1{$m}) and not defined($VAR1{$m})) {
                        print COMBINED_175 $VAR2{$m};
                }
        }
}

system("sort -k2n -k3n markers.combined.c195.txt.not_sorted > markers.combined.c195.txt");
system("rm markers.combined.c195.txt.not_sorted");

system("sort -k2n -k3n markers.combined.c190.txt.not_sorted > markers.combined.c190.txt");
system("rm markers.combined.c190.txt.not_sorted");

system("sort -k2n -k3n markers.combined.c175.txt.not_sorted > markers.combined.c175.txt");
system("rm markers.combined.c175.txt.not_sorted");


		






