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
    Intentionally loading file not exists
    File Not Found
    Loading existing and valid SBML file
    numOfRows: 26
    numOfSpecies: 6
    ...
