* How to run a test script:
  (1) If you are in a example directory under libsbml distirbution:
      % cp ../../build/src/bindings/python/libsbmlsim.py .
      % cp ../../build/src/bindings/python/_libsbmlsim.so .
      % python Test.py

  (2) If you are in installed directory (ex. /usr/local/share/libsbml/python)
      % python Test.py

* Output from this test program will be:
    numOfRows: 26
    numOfSpecies: 6
    numOfParameters: 12
    numOfCompartments: 1
    timeName: time
    Species Name
      CDK1
      PLK1
      APC
      inact_CDK1
      inact_PLK1
      inact_APC
    Parameters Name
      a1
      a2
      ...
