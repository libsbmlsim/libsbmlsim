use strict;
use warnings;
use utf8;
use 5.010;

use blib './blib';
use libsbmlsim;
use Data::Dumper;

my $xml;
open my $fh, "<src/MAPK.xml";
for (<$fh>) {
  $xml .= $_;
}
close $fh;

my $result = libsbmlsim::simulateSBMLFromString($xml, 4000.0, 0.1, 100, 1, 41, 0);
say $result->getNumOfRows;
