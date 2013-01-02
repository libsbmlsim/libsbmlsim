* How to compile & run a test program:
  [Windows]
    1. Copy following files under
         "C:\Program Files\libsbmlsim-1.1\share\libsbmlsim\csharp"       (32bit)
         "C:\Program Files (x86)\libsbmlsim-1.1\share\libsbmlsim\csharp" (64bit)
       to your temporary directory (say, your Desktop\temp).
        - Test.cs
        - libsbmlsimcs.dll
        - libsbmlsimcsP.dll
    2. Copy following .dll files under
         "C:\Program Files\libsbmlsim-1.1\bin"       (32bit)
         "C:\Program Files (x86)\libsbmlsim-1.1\bin" (64bit)
       to your temporary directory (Desktop\temp).
       These DLLs are required to be placed in the same directory of above files.
        - libbz2.dll
        - libiconv.dll
        - libsbml.dll
        - libxml2.dll
        - zlib1.dll
    3. Copy sample SBML model (sample.xml) files under
         "C:\Program Files\libsbmlsim-1.1\share\libsbmlsim"       (32bit)
         "C:\Program Files (x86)\libsbmlsim-1.1\share\libsbmlsim" (64bit)
       to your Desktop (not Desktop\temp).
    4. Open cmd.exe and navigate to Desktop\temp like:
       cd Desktop\temp
    5. Compile Test.cs with following commandline:
       csc.exe /reference:libsbmlsimcsP.dll Test.cs
       (csc.exe is usually located under
        "C:\Windows\Microsoft.NET\Framework\v2.0.50727" )
    6. Run "Test.exe".
       Test.exe
    Note: If you have VisualStudio installed, you can compile the test code
    on VisualStudio. Create a new Project and add Test.cs to your project as
    a source code, add libsbmlsimP.dll as a reference file, and then "Build".
    To run your compiled C# binary, above 7 DLLs should have to be placed on
    the same directory of Test.exe.

  [Unix (Linux and MacOSX) systems]
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
