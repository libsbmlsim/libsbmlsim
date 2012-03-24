use strict;
use warnings;
use utf8;
use 5.010;

use lib '.';
use libsbmlsim;

my $xml;
open my $fh, "<../sample.xml";
for (<$fh>) {
  $xml .= $_;
}
close $fh;

my $result = libsbmlsim::simulateSBMLFromString($xml, 4000.0, 0.1, 100, 1, 41, 0);
say $result->getNumOfRows;
