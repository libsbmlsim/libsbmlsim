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
$:.unshift File.join(File.dirname(__FILE__), ".")
require 'libsbmlsim'

#res = Libsbmlsim::simulateSBMLFromFile('MAPK.xml', 4000.0, 0.1, 100, 1, Libsbmlsim::MTHD_RUNGE_KUTTA, 0)
#puts res.getErrorMessage if res.isError
res = Libsbmlsim::simulateSBMLFromFile('../sample.xml', 25.0, 0.01, 100, 1, Libsbmlsim::MTHD_RUNGE_KUTTA, 0)

numOfRows = res.getNumOfRows
puts "numOfRows: #{numOfRows}"

numOfSp = res.getNumOfSpecies
puts "numOfSpeies: #{numOfSp}"

numOfParam = res.getNumOfParameters
puts "numOfParameters: #{numOfParam}"

numOfComp = res.getNumOfCompartments
puts "numOfCompartments: #{numOfComp}"

timeName = res.getTimeName
puts "timeName: #{timeName}"

puts "Species Name"
numOfSp.times do |n|
  sname = res.getSpeciesNameAtIndex(n)
  puts "  #{sname}"
end

puts "Parameter Name"
numOfParam.times do |n|
  pname = res.getParameterNameAtIndex(n)
  puts "  #{pname}"
end

puts "Compartment Name"
numOfComp.times do |n|
  cname = res.getCompartmentNameAtIndex(n)
  puts "  #{cname}"
end

puts "Values"
numOfRows.times do |i|
  t = res.getTimeValueAtIndex(i)
  print "#{t}"
  numOfSp.times do |j|
    sname = res.getSpeciesNameAtIndex(j)
    v = res.getSpeciesValueAtIndex(sname, i)
    print " #{v}"
  end
  print "\n"
  break if i == 25;
end

