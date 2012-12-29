#!/bin/sh

swig -const -perl libsbmlsim.i && \
/opt/local/bin/perl Makefile.PL && \
make
#gcc -fPIC -c libsbmlsim_wrap.c -I`perl -e 'use Config; print $Config{archlib};'`/CORE && \
#gcc -shared *.o src/CMakeFiles/sbmlsim.dir/*.o -L/usr/local/lib -lsbml -o libsbmlsim.so
