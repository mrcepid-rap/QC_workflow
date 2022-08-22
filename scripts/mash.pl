#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open (SAMP, "wes_samples.470k.txt") || die "Cannot open file: $!";

my %samps;

foreach (<SAMP>) {

	chomp $_;
	$samps{$_} = {'PTV' => 0, 'SYN' => 0, 'MISSENSE' => 0};

}

open (LST, "per_indv_list.txt") || die "Cannot open file: $!";

foreach (<LST>) {

	chomp $_;
	open (CURR, "$_") || die "Cannot open file: $!";

	foreach my $line (<CURR>) {

		chomp $line;
		my @data = split("\t", $line);

		if ($data[3] eq "0/1") {
			$samps{$data[0]}{$data[1]} += 1;
		} else {
			$samps{$data[0]}{$data[1]} += 2;
		}
		
	}

	close CURR;

}
#print Dumper(\%samps);
print "eid\tPTV\tMISSENSE\tSYN\n";
foreach (keys %samps) {

	print "$_\t$samps{$_}{PTV}\t$samps{$_}{MISSENSE}\t$samps{$_}{SYN}\n";

}
