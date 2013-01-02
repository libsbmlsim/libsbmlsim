* How to compile & run a test program:
  (1) If you are in a example directory under libsbml distirbution:
      % cp ../../build/src/bindings/java/libsbmlsimj.jnilib .
      % javac -cp ../../build/src/bindings/java/libsbmlsimj.jar   Test.java
      % java  -cp ../../build/src/bindings/java/libsbmlsimj.jar:. Test

  (2) If you are in installed directory (ex. /usr/local/share/libsbml/java)
      % javac -cp libsbmlsimj.jar   Test.java
      % java  -cp libsbmlsimj.jar:. Test

* Output from this test program will be:
    Simulate SBML model from File
    Loading existing and valid SBML file
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
    0.0 0.0 0.0 0.0 1.0 1.0 1.0
    1.0 0.10000000000000007 7.751561133890252E-7 1.3666739792994853E-48 1.0 1.0 1.0
    2.0 0.20000000000000015 3.6291243271494456E-4 6.280884943555849E-27 1.0 1.0 1.0
    3.0 0.3000000000000002 0.012634838454320584 2.0487086866095784E-14 1.0 1.0 1.0
    4.0 0.4000000000000003 0.13239183229481652 4.496394132436305E-6 1.0 1.0 1.0
    ...
