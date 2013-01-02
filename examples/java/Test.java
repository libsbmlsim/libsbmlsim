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
//package jp.ac.keio.bio.fun.libsbmlsim;
import jp.ac.keio.bio.fun.libsbmlsim.*;

public class Test {
  public static void main(String[] args) {
    System.loadLibrary("sbmlsimj");

    /*
       SBMLDocument d = libsbml.readSBML("./src/MAPK.xml");
       String docstr = d.toSBML();
       myResult result = libsbmlsim.simulateSBMLFromString(docstr, 4000.0, 0.1, 100, 1, 41, 0);
    */
    System.out.println("Simulate SBML model from File");
    /*
       System.out.println("Intentionally loading file not exists");
       myResult result_with_error = libsbmlsim.simulateSBMLFromFile("does_not_Exists.xml", 4000.0, 0.1, 100, 1, libsbmlsim.MTHD_RUNGE_KUTTA, 0);
       if (result_with_error.isError()) {
         System.out.println(result_with_error.getErrorMessage());
       }
    */
    System.out.println("Loading existing and valid SBML file");
    myResult result = libsbmlsim.simulateSBMLFromFile("../sample.xml", 25.0, 0.01, 100, 1, libsbmlsim.MTHD_RUNGE_KUTTA, 0);

    //libsbmlsim.print_result(result);
    //libsbmlsim.write_csv(result, "result.csv"); // Output result to a csv file
    //libsbmlsim.write_result(result, "result.dat"); // Output result to a file

    int numOfRows = result.getNumOfRows();
    System.out.println("numOfRows: " + numOfRows);

    int numOfSp = result.getNumOfSpecies();
    System.out.println("numOfSpecies: " + numOfSp);

    int numOfParam = result.getNumOfParameters();
    System.out.println("numOfParameters: " + numOfParam);

    int numOfComp = result.getNumOfCompartments();
    System.out.println("numOfCompartments: " + numOfComp);

    String timeName = result.getTimeName();
    System.out.println("TimeName: " + timeName);

    System.out.println("Species Name:");
    for (int i = 0; i < numOfSp; i++) {
      String sname = result.getSpeciesNameAtIndex(i);
      System.out.println("  " + sname);
    }

    System.out.println("Parameter Name:");
    for (int i = 0; i < numOfParam; i++) {
      String pname = result.getParameterNameAtIndex(i);
      System.out.println("  " + pname);
    }

    System.out.println("Compartment Name:");
    for (int i = 0; i < numOfComp; i++) {
      String cname = result.getCompartmentNameAtIndex(i);
      System.out.println("  " + cname);
    }

    System.out.println("Values:");
    for (int i = 0; i < numOfRows; i++) {
      double t = result.getTimeValueAtIndex(i);
      System.out.print(t);
      for (int j = 0; j < numOfSp; j++) {
        String sname = result.getSpeciesNameAtIndex(j);
        double val = result.getSpeciesValueAtIndex(sname, i);
        System.out.print(" " + val);
      }
      System.out.println();
      if (i == 25)
        break;
    }
  }
}
