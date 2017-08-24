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

* How to run a grid search script:
  At first, you have to install python bindings of libSBML.
  Please refer to libSBML installation guide[1] to install it.
  (1) If you are in a example directory under libsbml distirbution:
      % cp ../../build/src/bindings/python/libsbmlsim.py .
      % cp ../../build/src/bindings/python/_libsbmlsim.so .
      % python gridSearch.py

  (2) If you are in installed directory (ex. /usr/local/share/libsbml/python)
      % python gridSearch.py

  [1] http://sbml.org/Software/libSBML/Downloading_libSBML

* Output from this test program will be:
    Experimental Data
    ('time', 'A', 'B', 'k', 'default')
    [(  0.,   1.00000000e+00,  0.        ,  0.8,  1.)
     (  1.,   4.49329095e-01,  0.5506709 ,  0.8,  1.)
     (  2.,   2.01896636e-01,  0.79810336,  0.8,  1.)
     (  3.,   9.07180327e-02,  0.90928197,  0.8,  1.)
     (  4.,   4.07622516e-02,  0.95923775,  0.8,  1.)
     (  5.,   1.83156656e-02,  0.98168433,  0.8,  1.)
     (  6.,   8.22976146e-03,  0.99177024,  0.8,  1.)
     (  7.,   3.69787127e-03,  0.99630213,  0.8,  1.)
     (  8.,   1.66156115e-03,  0.99833844,  0.8,  1.)
     (  9.,   7.46587770e-04,  0.99925341,  0.8,  1.)
     ( 10.,   3.35463607e-04,  0.99966454,  0.8,  1.)]
    Parameter k = 0.0       Error = 17.2431693098
    Parameter k = 0.05      Error = 9.5415096976
    Parameter k = 0.1       Error = 5.57659227977
    ...
    Parameter k = 0.8       Error = 4.83392664417e-35
    Parameter k = 0.85      Error = 0.00217235929845
    Parameter k = 0.9       Error = 0.00794285269809
    Parameter k = 0.95      Error = 0.0163981861688
    Parameter k = 1.0       Error = 0.0268413252801
