#!/usr/bin/python
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

from libsbmlsim import *

def main():
## readSBMLFromString example
#  f = open('./src/MAPK.xml')
#  xml = f.read()
#  f.close()
#  result = simulateSBMLFromString(xml, 4000.0, 0.1, 100, 1, 41, 0)

## error handling
#  result = simulateSBMLFromFile('hogeMAPK.xml', 4000.0, 0.1, 100, 1, MTHD_RUNGE_KUTTA, 0)
#  if result.isError():
#    print result.getErrorMessage()
  result = simulateSBMLFromFile('../sample.xml', 25.0, 0.01, 100, 1, MTHD_RUNGE_KUTTA, 0)

  numOfRows = result.getNumOfRows()
  print "numOfRows: " + str(numOfRows)

  numOfSp = result.getNumOfSpecies()
  print "numOfSpecies: " + str(numOfSp)

  numOfParam = result.getNumOfParameters()
  print "numOfParameters: " + str(numOfParam)

  numOfComp = result.getNumOfCompartments()
  print "numOfCompartments: " + str(numOfComp)

  timeName = result.getTimeName()
  print "timeName: " + timeName

  print "Species Name"
  for i in range(numOfSp):
    sname = result.getSpeciesNameAtIndex(i)
    print "  " + sname

  print "Parameters Name"
  for i in range(numOfParam):
    pname = result.getParameterNameAtIndex(i)
    print "  " + pname

  print "Compartment Name"
  for i in range(numOfComp):
    cname = result.getCompartmentNameAtIndex(i)
    print "  " + cname

  print "Values"
  for i in range(numOfRows):
    t = result.getTimeValueAtIndex(i)
    print t,
    for j in range(numOfSp):
      sname = result.getSpeciesNameAtIndex(j)
      v = result.getSpeciesValueAtIndex(sname, i)
      print str(v),
    print
    if i == 25:
      break

if __name__ == "__main__":
  main()

