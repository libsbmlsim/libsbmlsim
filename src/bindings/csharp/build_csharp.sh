#!/bin/sh

swig -includeall -csharp -namespace "libsbmlsimcs" -dllimport "libsbmlsimcs" libsbmlsim.i
gcc -O2 -I../ -fPIC -c libsbmlsim_wrap.c
gcc -shared *.o ../../../build/src/CMakeFiles/sbmlsim.dir/*.o -L/usr/local/lib -lsbml -o libsbmlsimcs.bundle
gmcs -target:library -out:libsbmlsimcsP.dll libsbmlsim.cs myResult.cs LibsbmlsimErrorCode.cs SWIGTYPE_p_p_char.cs SWIGTYPE_p_double.cs SWIGTYPE_p_BOOLEAN.cs libsbmlsimPINVOKE.cs

# Test compile
#cp ../../../examples/csharp/Test.cs .
#cp ../../../examples/sample.xml .
#gmcs -reference:libsbmlsimcsP.dll Test.cs
#mono Test.exe
