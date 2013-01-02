           LibSBMLSim: The library for simulating SBML models

                     LibSBMLSim development team
             http://fun.bio.keio.ac.jp/software/libsbmlsim/
                  mailto:sbmlsim@fun.bio.keio.ac.jp

-- Last modified: Wed, 02 Jan 2013 23:03:21 +0900

* Overview
  LibSBMLSim is a library for simulating an SBML model which contains
  Ordinary Differential Equations (ODEs). LibSBMLSim provides simple
  command-line tool and several APIs to load an SBML model, perform
  numerical integration (simulate) and export its results.
  Both explicit and implicit methods are supported on libSBMLSim.
  LibSBMLSim is confirmed to pass all SBML Level-2 Version 4 and Level-3
  Version 1 test cases (sbml-test-cases-2.0.2.zip, available from
  http://sourceforget.net/projects/sbml/files/test-suite/2.0.2/).
  The libSBMLSim code is portable. It is written in C programming language
  (ANSI C89) and it does not depend on other third-party libraries
  except libSBML(*1).
  The library should build and work without serious troubles on Unix
  based operating systems (Linux, MacOSX and FreeBSD) and on Windows
  (with Visual C++).
  LibSBMLSim also provides several language bindings like Java, Python,
  C# and Ruby. Perl binding is already included in the source tree, but
  is not able to create through the single build process (see the
  description below).

  (*1 libSBML: http://sbml.org/Software/libSBML)

  LibSBMLSim can be used to create your own SBML capable simulator,
  plug-in, web based application and web services. The API is quite
  straight forward. You can run a simulation and generate a result
  file in Comma Separated Values (CSV) with a few lines of codes.
  === Python ============================
    from libsbmlsim import *
    r = simulateSBMLFromFile('sbml.xml', 20.0, 0.1, 10, 0, MTHD_RUNGE_KUTTA, 0)
    write_csv(r, 'result.csv')
  =======================================

  Please see the 'API.txt' and 'examples' directory for further information
  on libSBMLSim APIs.
  
* Installation
- Dependencies
  LibSBMLSim requires libSBML to be installed on your system.
  Please follow the instruction on (*1) and install libSBML.
  
- Required software packages to compile libSBMLSim
  CMake(*2) is required to compile libSBMLSim. Please download
  and install CMake-2.8.7 or above from (*2) before building libSBMLSim.
  If you want to use language bindings of libSBMLSim, please
  download and install SWIG-2.0.4 or above from (*3).
  (Note: If you installed SWIG from MacPorts, please install
         swig-java, swig-python, swig-ruby, swig-csharp which are required
         to compile language bindings for libSBMLSim.)

  (*2 CMake: http://cmake.org/)
  (*3 SWIG:  http://swig.org/)

- How to build libSBMLSim
  1. Extract the archive file
   % tar xvzf libsbmlsim-1.1.0.tar.gz (for tar ball)
   % unzip libsbmlsim-1.1.0.zip       (for zip archive)
  2. Compile
   % cd libsbmlsim-1.1.0/build
   % cmake ..
   % ccmake .
     CUI from cmake will be launched. Please confirm that
     cmake have automatically detected the installed location of
     libSBML. You can check the installed location from the
     following values:
     (ex. on MacOSX)
       LIBSBML_INCLUDE_DIR            /usr/local/include
       LIBSBML_LIBRARY                /usr/local/lib/libsbml.dylib

     If libSBML is not detected automatically, you can manually
     specify the installed location through this menu.

     If you want to build language bindings, please turn on the
     corresponding compile option.
       WITH_JAVA     ... build with Java bindings
       WITH_PYTHON   ... build with Python bindings
       WITH_RUBY     ... build with Ruby bindings
       WITH_CSHARP   ... build with C# bindings

     Once you press [c] key, cmake will run the configure procedure
     and tries to detect SWIG, Java, Python, C# and Ruby (depending on
     which language bindings you enabled). Once the configuration
     is done, press [g] key and cmake will generate Makefile.
     After Makefile is generated, just run

   % make
   % sudo make install
     which will compile and install the library, command-line tool
     and header files on your system. Default prefix (install
     directory) is
       - /usr/local                  ... on Linux and MacOSX
       - C:\Program Files\libsbmlsim ... on Windows
     (Note: You can change the prefix from the UI of ccmake)

- Installed files
  Following files are installed on your system.
  = Unix based systems (Linux, MacOSX, etc.)
    $prefix/bin/simulateSBML         ... SBML simulator
    $prefix/lib/libsbmlsim-static.a  ... Static library 
               /libsbmlsim.dylib     ... Dynamic library (on MacOSX)
               /libsbmlsim.so        ... Dynamic library (on Linux)
    $prefix/include/libsbmlsim       ... Header files
    $prefix/share/libsbmlsim/        ... Sample files (SBML, results)
                            /c       ... Sample C code
                            /cpp     ... Sample C++ code
                            /csharp  ... Sample C# code and language bindings
                            /java    ... Sample Java code and language bingings
                            /python  ... Sample Python code and language bingings
                            /ruby    ... Sample Ruby code and language bingings

  = Windows
    $prefix\bin\simulateSBML.exe     ... SBML simulator
    $prefix\bin\sbmlsim.dll          ... Dynamic library
    $prefix\lib\sbmlsim-static.lib   ... Static library
               \sbmlsim.lib          ... Lib file for dynamic library
    $prefix\include\libsbmlsim       ... Header files
    $prefix\share\libsbmlsim\        ... Sample files (SBML, results)
                            \c       ... Sample C code
                            \cpp     ... Sample C++ code
                            \cpp     ... Sample C++ code
                            \csharp  ... Sample C# code and language bindings
                            \java    ... Sample Java code and language bingings
                            \python  ... Sample Python code and language bingings
    (Note: Ruby binding is not supported on Windows)

* Usage
- simulateSBML
  simulateSBML is a simple SBML simulator which accept SBML file as
  an input, and then output "out.csv" as a result.
  Usage: simulateSBML [option] filename(SBML)
    -t # : specify simulation time (ex. -t 100 )
    -s # : specify simulation step (ex. -s 100 )
    -d # : specify simulation delta (ex. -d 0.01 [default:1/4096])
           dt is calculated in (delta)*(time)/(step)
    -l   : use lazy method for integration
    -n   : do not use lazy method
    -a   : print Species Value in Amount
    -A # : specify absolute tolerance for variable stepsize (ex. -A 1e-03 [default:1e-04])
    -R # : specify relative tolerance for variable stepsize (ex. -R 0.1 [default:1e-04])
    -M # : specify the max change rate of stepsize (ex. -M 1.5 [default:2.0])
    -m # : specify numerical integration algorithm (ex. -m 3 )
        1: Runge-Kutta
        2: AM1 & BD1 (implicit Euler)
        3: AM2 (Crank Nicolson)
        4: AM3
        5: AM4
        6: BD2
        7: BD3
        8: BD4
        9: AB1 (explicit Euler)
       10: AB2
       11: AB3
       12: AB4
       13: Runge-Kutta-Fehlberg
       14: Cash-Karp
       (AM: Adams-Moulton, BD: Backward-Difference, AB: Adams-Bashforth.
        Number after synonim specifies the order of integration.
        For example, AM2 is "2nd order Adams-Moulton" method)

- Scripts for "SBML test cases"
  LibSBMLSim provides scripts to easily run SBML test cases (*4)
  and compare the results with it. Generated results are compatible
  with Online SBML Test Suite (*4), so you can run all tests with
  this scripts and upload the results to Online SBML Test Suite.
  The scripts are not installed, you will find them under "testcases"
  directory in the extracted source directory (libsbmlsim/testcases/). 

    libsbmlsim/testcases/simulateSBML  ... SBML simulator
                        /runall.sh     ... Script which will run all tests
                        /compare.pl
                        /genresult.pl

  "simulateSBML" simulates SBML model and generates simulation result
  as a CSV file, which is identical with the one installed under
  $prefix/bin . "runall.sh" will call simulateSBML for all SBML test
  cases, and compare the result with the one from SBML test cases.
  "compare.pl" and "genresult.pl" are scripts which will support some
  functions called from runall.sh.
  The SBML test cases are not included in this distribution, so please
  download them from (*5). After downloading sbml-test-cases-X.Y.Z.zip,
  unzip the archive and move (or copy) "cases/" directory to 
  libsbmlsim/testcases directory. The directory structure will be:

    libsbmlsim/testcases/simulateSBML  ... SBML simulator
                        /runall.sh     ... Script which will run all tests
                        /compare.pl
                        /genresult.pl
                        /cases/semantic/00001 ... Test case 1
                        /cases/semantic/00002 ... Test case 2
                        /cases/semantic/00003 ... Test case 3
                        /cases/semantic/...

  Following command will test all 980 tests and print out the
  results, whether the simulation result matches with the result
  with the one from SBML test cases.

    % ./runall.sh
    00001: 5 : 50 : [S1,S2] : [S1,S2] : [S1,S2]
      print amount
      time:5 step:50 dt:0.000024
      Model 00001 ... [OK]
    00002: 5.0 : 50 : [S1,S2] : [S1,S2] : [S1,S2]
      print amount
      time:5 step:50 dt:0.000024
      Model 00002 ... [OK]
    00003: 5.0 : 50 : [S1,S2] : [S1,S2] : [S1,S2]
      print amount
      time:5 step:50 dt:0.000024
      Model 00003 ... [OK]

  If the simulation result doesn't match with the one from
  SBML test cases, then the result will be marked as "[NG]".

  (*4 Online SBML Test Suite: http://sbml.org/Facilities/Online_SBML_Test_Suite)
  (*5 SBML-test-cases-2.0.2 : http://sourceforge.net/projects/sbml/files/test-suite/2.0.2/sbml-test-cases-2.0.2.zip/download)

* LibSBMLSim API and its language bindings
  Example usage of libSBMLSim APIs are as follows.
  Please see the 'API.txt' and 'examples' directory for further information
  on libSBMLSim APIs.

- C, C++ API
  === C code ============================
  #include "libsbmlsim/libsbmlsim.h"
  ...
  /*
   * Simulate sbml.xml to time=20 with dt=0.1, print_interval=10
   * by 4th-order Runge-Kutta Method.
   */
  myResult *r = simulateSBMLFromFile("sbml.xml", 20, 0.1, 10, 0, MTHD_RUNGE_KUTTA, 0);
  /*
   * Export simulation result as CSV file
   */
  write_csv(r, "result.csv"); 
  /*
   * Free Result object
   */
  free_myResult(r);
  =====================================

- Java, Python, Ruby bindings
  LibSBMLSim API is also provided for several language bindings.
  === Python ============================
  from libsbmlsim import *
  r = simulateSBMLFromFile('sbml.xml', 20.0, 0.1, 10, 0, MTHD_RUNGE_KUTTA, 0)
  write_csv(r, 'result.csv')
  =======================================
  
  === Java ==============================
  import jp.ac.keio.bio.fun.libsbmlsim.*;
  ...
  System.loadLibrary("sbmlsimj");
  myResult r = libsbmlsim.simulateSBMLFromFile("sbml.xml", 20.0, 0.1, 10, 0, libsbmlsim.MTHD_RUNGE_KUTTA, 0);
  libsbmlsim.write_csv(r, "result.csv");
  =======================================
  
  === Ruby ==============================
  require 'libsbmlsim'
  r = Libsbmlsim::simulateSBMLFromFile('sbml.xml', 20.0, 0.1, 10, 0, Libsbmlsim::MTHD_RUNGE_KUTTA, 0)
  Libsbmlsim::write_csv(r, 'result.csv')
  =======================================
  
  === C# ================================
  using System;
  public class Test
  {
    static void Main()
      {
        myResult result = libsbmlsim.simulateSBMLFromFile("sbml.xml", 20.0, 0.1, 10, 0, libsbmlsim.MTHD_RUNGE_KUTTA, 0);
        libsbmlsim.write_csv(result, "test.csv");
      }
  }
  =======================================
  
  Please see the 'API.txt' and 'examples' directory for further information.
  The 'examples' directory contains sample code for test application
  in several programming languages (C, C++, Java, Python, Ruby, C# and Perl).

Have fun!
-- 
LibSBMLSim development team <sbmlsim@fun.bio.keio.ac.jp>
