#!/usr/bin/env python
#  <!--------------------------------------------------------------------------
#  This file is part of libSBMLSim.  Please visit
#  http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
#  information about libSBMLSim and its latest version.
# 
#  Copyright (C) 2011-2017 by the Keio University, Yokohama, Japan
# 
#  This library is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation.  A copy of the license agreement is provided
#  in the file named "LICENSE.txt" included with this software distribution.
#  ---------------------------------------------------------------------- -->

#
# This script will run a grid search for a given global parameter and
# calculate error from a given experimental data.
#
from libsbml import *
from libsbmlsim import *
from numpy import genfromtxt
from numpy import linspace

def main():
    # Experimental Data
    datafile = './data.csv'
    data = importData(datafile)
    # SBML model
    modelfile = './simple.xml'
    d = readSBML(modelfile)
    # Simulation
    simulation_time = 10
    dt = 0.1
    # 1D Grid Search
    param_id = 'k'
    min_value = 0
    max_value = 1.0
    grid_spacing = 0.05
    grid_points = (max_value - min_value) / grid_spacing + 1
    for param_value in linspace(min_value, max_value, num=grid_points):
        updateParameter(d, param_id, param_value)
        result = integrate(d, simulation_time, dt)
        error = calcError(result, data)
        print "Parameter", param_id, "=", param_value, "\tError =", error

def updateParameter(sbml, parameter_id, new_value):
    m = sbml.getModel()
    p = m.getParameter(parameter_id)
    p.setValue(new_value)

def calcError(result, data):
    error = 0.0
    numOfRows = result.getNumOfRows()
    numOfSp = result.getNumOfSpecies()
    for i in range(numOfRows):
        t = result.getTimeValueAtIndex(i)
        for di in data:
            if di['time'] == t:
                for j in range(numOfSp):
                    sname = result.getSpeciesNameAtIndex(j)
                    value_sim = result.getSpeciesValueAtIndex(sname, i)
                    error += (value_sim - di[sname])**2

    return error

def integrate(sbml, simulation_time, dt):
    result = simulateSBMLFromString(sbml.toSBML(), simulation_time, dt, 1, 0, MTHD_RUNGE_KUTTA, 0)
    if result.isError():
        print result.getErrorMessage()

    return result

def debugPrint(result):
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

def importData(datafile):
    data = genfromtxt(datafile, delimiter=',', names=True)
    print "Experimental Data"
    print data.dtype.names
    print data
    return data

if __name__ == "__main__":
  main()

