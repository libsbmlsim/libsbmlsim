#!/usr/bin/python

from libsbmlsim import *

def main():
#  f = open('./src/MAPK.xml')
#  xml = f.read()
#  f.close()

#  result = simulateSBMLFromString(xml, 4000.0, 0.1, 100, 1, 41, 0)
  result = simulateSBMLFromFile('../../MAPK.xml', 4000.0, 0.1, 100, 1, MTHD_RUNGE_KUTTA, 0)

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
    if i == 10:
      break

if __name__ == "__main__":
  main()

