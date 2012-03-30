* How to run a test script:
  (1) If you are in a example directory under libsbml distirbution:
      % cp ../../build/src/bindings/ruby/libsbmlsim.bundle .
      % ruby Test.rb

  (2) If you are in installed directory (ex. /usr/local/share/libsbml/ruby)
      % ruby Test.rb

* Output from this test program will be:
    numOfRows: 26
    numOfSpeies: 6
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
    Parameter Name
      a1
      a2
      ...
