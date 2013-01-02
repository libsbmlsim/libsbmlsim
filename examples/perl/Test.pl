#
#  <!--------------------------------------------------------------------------
#  This file is part of libSBMLSim.  Please visit
#  http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
#  information about libSBMLSim and its latest version.
# 
#  Copyright (C) 2011-2013 by the Keio University, Yokohama, Japan
# 
#  This library is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation.  A copy of the license agreement is provided
#  in the file named "LICENSE.txt" included with this software distribution.
#  ---------------------------------------------------------------------- -->
use strict;
use warnings;
use utf8;
use 5.010;

use lib '.';
use libsbmlsim;

my $xml;
my $resfile = "test.dat";
open my $fh, "<../sample.xml";
for (<$fh>) {
  $xml .= $_;
}
close $fh;

my $result = libsbmlsim::simulateSBMLFromString($xml, 25.0, 0.01, 10, 1, libsbmlsim::MTHD_RUNGE_KUTTA, 0);
#say $result->getNumOfRows;
libsbmlsim::write_result($result, $resfile);
say "Simulation result is written to $resfile.";
