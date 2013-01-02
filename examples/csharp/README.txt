* How to compile & run a test program:
  (1) If you are in a example directory under libsbml distirbution:
      % cp ../../build/src/bindings/csharp/CSharpBinaries/libsbmlsimcs.bundle .
      % cp ../../build/src/bindings/csharp/CSharpBinaries/libsbmlsimcsP.dll .
      % gmcs -reference:libsbmlsimcsP.dll Test.cs
      % mono Test.exe

  (2) If you are in installed directory (ex. /usr/local/share/libsbml/java)
      % gmcs -reference:libsbmlsimcsP.dll Test.cs
      % mono Test.exe

* Output from this test program will be:
    numOfRows: 26
    numOfSpecies: 6
    numOfParameters: 12
    numOfCompartments: 1
    TimeName: time
    Species Name:
      CDK1
      PLK1
      APC
      inact_CDK1
      inact_PLK1
      inact_APC
    Parameter Name:
      a1
      a2
      ...
    Compartment Name:
      default
    Values:
    0 0 0 0 1 1 1
    1 0.1 7.75156113389025E-07 1.36667397929948E-48 1 1 1
    2 0.2 0.000362912432714945 6.28088494355585E-27 1 1 1
    3 0.3 0.0126348384543206 2.04870868660958E-14 1 1 1
    4 0.4 0.132391832294817 4.49639413243631E-06 1 1 1
    ...
    Simulation result is written to test.dat.
