#!/bin/sh

swig -java libsbmlsim.i && \
gcc -c -Wall -I/System/Library/Frameworks/JavaVM.framework/Headers -I./src libsbmlsim_wrap.c && \
g++ -shared *.o src/CMakeFiles/sbmlsim.dir/*.o -o libsbmlsim.jnilib -Lsrc -L/usr/local/lib -lsbml && \
mv Test.java Test.java.org && \
LC_ALL=C javac *.java && \
rm sbmlsim.jar && \
LC_ALL=C jar cvf sbmlsim.jar *.class && \
mv Test.java.org Test.java && \
javac -cp .:sbmlsim.jar:libsbmlj.jar Test.java && \
java -cp .:sbmlsim.jar:libsbmlj.jar Test
