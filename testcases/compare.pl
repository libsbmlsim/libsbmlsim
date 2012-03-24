#!/usr/bin/perl

use strict;
use warnings;

my $HOME = $ENV{"HOME"};
my $BaseDir="./cases/semantic";

my $verbose = 0;
my $red     = "\x1b[1;31m";
my $green   = "\x1b[1;32m";
my $default = "\x1b[0m";

my $arg = $ARGV[0];
if ($arg eq "-v") {
	$verbose = 1;
	shift(@ARGV);
}

foreach my $num (@ARGV) {
	my (%hash1, %hash2);
	my $msg;

### Open settings
	my $Ta = 0.0;
	my $Tr = 0.0;
	open(SET, "$BaseDir/$num/$num-settings.txt") or die;
	while(<SET>) {
		chomp;
		if (/absolute: (\S+)/) {
			$Ta = $1;
		}
		if (/relative: (\S+)/) {
			$Tr = $1;
		}
	}
	close(SET);

### My simulation result
	my $linenum = 0;
	if (!open(IN1, "$num-results.csv")) {
		print "Model $num ... ". "[No result]\n";
		next;
	}
	while(<IN1>) {
		chomp;
		my @array1 = split(/,/);
		push(@{$hash1{$linenum}},@array1);  # add to %dataHash $key array
		$linenum++;
	}
	close(IN1);

### Result from SBML test case
	$linenum = 0;
	open(IN2, "$BaseDir/$num/$num-results.csv") or die "$!\n";
	while(<IN2>) {
		chomp;
		my @array2 = split(/,/);
		push(@{$hash2{$linenum}},@array2);  # add to %dataHashの$key array
		$linenum++;
	}
	close(IN2);

# print column name first.
	if ($verbose) {
		for (my $i = 0; $i < @{$hash1{0}}; $i++) {
			print "$hash1{0}[$i]\t";
		}
		print "\n";
	}

# Ta stand for the absolute tolerance for this test case,
# Tr stand for the relative tolerance for this test case,
# Cij stand for the expected correct value for row i, column j, of the result data set
# Uij stand for the the user's uploaded result value.
# 
# |Cij − Uij| ≤ ( Ta + Tr × |Cij| )
	foreach my $index (sort {$a <=> $b} keys %hash1){
		if ($index != 0) {
			print "$hash1{$index}[0]\t" if $verbose;
			my $num_error = 0;
			for (my $i = 1; $i < @{$hash1{$index}}; $i++) {
				my $delta = $hash2{$index}[$i] - $hash1{$index}[$i];
				if ( abs($delta) > $Ta + $Tr*abs($hash2{$index}[$i]) ) {
					$num_error++;
					print "$red". "$hash1{$index}[$i]". "$default\t" if $verbose;
				} else {
					print "$hash1{$index}[$i]\t" if $verbose;
				}
			}
			$msg = "$green"."OK";
			if ($num_error > 0) {
				$msg = "$red" . "NG";
			}
			print "[$msg". "$default]\n" if $verbose;
		}
	}
	print "Model $num ... ". "[$msg". "$default]\n";
}
