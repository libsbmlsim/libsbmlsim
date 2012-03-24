#!/usr/bin/perl
# ./genresult.pl csv_file list_of_valuables(in comma separated list)
# (ex.) ./genresult.pl out.csv S1,S2,k1

use strict;
use warnings;

my $vars = "time," . $ARGV[1];
my $steps = $ARGV[2];
print $vars . "\n";
pop @ARGV;
pop @ARGV;

my %hash = ();
my @nameArray;
my $linenum = 0;

while(<>) {
  $linenum++;
  s/\s+//g;
  if (/time/) {  # if it is first line
    @nameArray = split(/,/);
  } else {
    if ($linenum > ($steps+2) ) {
      exit;
    }
    my $n = 0;
    foreach my $i (split(/,/)) {
      $i =~ s/1.#INF/inf/;
      $i =~ s/-1.#INF/-inf/;
      $i =~ s/1.#QNAN/nan/;
      $hash{$nameArray[$n]} = $i;
      $n++;
    }
    my $m = 0;
    foreach my $j (split(/,/, $vars)) {
      if ($m > 0) {
        print ",";
      }
      print "$hash{$j}";
      $m++;
    }
    print "\n";
  }
}
