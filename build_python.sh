#!/bin/sh

swig -python libsbmlsim.i && \
gcc -O2 -fPIC -c libsbmlsim_wrap.c -I/usr/include/python2.6 && \
gcc -shared *.o src/CMakeFiles/sbmlsim.dir/*.o -lpython -L/usr/local/lib -lsbml -o _libsbmlsim.so
