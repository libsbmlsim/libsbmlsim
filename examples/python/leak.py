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
from libsbmlsim import *
from numpy import linspace
import gc

def main():
    # SBML model
    modelfile = './mapk.xml'
    #d = readSBML(modelfile)
    # Simulation
    simulation_time = 4000
    dt = 1
    # 1D Grid Search
    min_value = 0.1
    max_value = 1000
    grid_spacing = 0.05
    grid_points = (max_value - min_value) / grid_spacing + 1
    param_id = 'K'
    for param_value in linspace(min_value, max_value, num=grid_points):
        result = integrateF(modelfile, simulation_time, dt)
        print "Parameter", param_id, "=", param_value
        #result.__swig_destroy__(result)  # free myResult object
        result.thisown = 0
        del result
        gc.collect()

def integrate(sbml, simulation_time, dt):
    #result = simulateSBMLFromString(sbml.toSBML(), simulation_time, dt, 1, 0, MTHD_RUNGE_KUTTA, 0)
    result = simulateSBMLFromString(sbml.toSBML(), simulation_time, dt, 1, 0, MTHD_EULER, 0)
    if result.isError():
        print result.getErrorMessage()

    return result

def integrateF(sbmlfile, simulation_time, dt):
    result = simulateSBMLFromFile(sbmlfile, simulation_time, dt, 1, 0, MTHD_RUNGE_KUTTA, 0)
    result.thisown = 0
    if result.isError():
        print result.getErrorMessage()

    return result

if __name__ == "__main__":
  main()

