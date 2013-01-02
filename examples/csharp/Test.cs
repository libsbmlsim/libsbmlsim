/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2013 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
using System;

public class Test
{
  static void Main()
  {
    myResult result;
    result = libsbmlsim.simulateSBMLFromFile("../sample.xml", 25.0, 0.01, 100, 1, libsbmlsim.MTHD_RUNGE_KUTTA, 0);
    int numOfRows = result.getNumOfRows();
    Console.WriteLine("numOfRows: " + numOfRows);

    int numOfSp = result.getNumOfSpecies();
    Console.WriteLine("numOfSpecies: " + numOfSp);

    int numOfParam = result.getNumOfParameters();
    Console.WriteLine("numOfParameters: " + numOfParam);

    int numOfComp = result.getNumOfCompartments();
    Console.WriteLine("numOfCompartments: " + numOfComp);

    String timeName = result.getTimeName();
    Console.WriteLine("TimeName: " + timeName);

    Console.WriteLine("Species Name:");
    for (int i = 0; i < numOfSp; i++) {
      String sname = result.getSpeciesNameAtIndex(i);
      Console.WriteLine("  " + sname);
    }

    Console.WriteLine("Parameter Name:");
    for (int i = 0; i < numOfParam; i++) {
      String pname = result.getParameterNameAtIndex(i);
      Console.WriteLine("  " + pname);
    }

    Console.WriteLine("Compartment Name:");
    for (int i = 0; i < numOfComp; i++) {
      String cname = result.getCompartmentNameAtIndex(i);
      Console.WriteLine("  " + cname);
    }

    Console.WriteLine("Values:");
    for (int i = 0; i < numOfRows; i++) {
      double t = result.getTimeValueAtIndex(i);
      Console.Write(t);
      for (int j = 0; j < numOfSp; j++) {
        String sname = result.getSpeciesNameAtIndex(j);
        double val = result.getSpeciesValueAtIndex(sname, i);
        Console.Write(" " + val);
      }
      Console.WriteLine();
      if (i == 25)
        break;
    }

    libsbmlsim.write_result(result, "test.dat");
    Console.WriteLine("Simulation result is written to test.dat.");
  }
}
