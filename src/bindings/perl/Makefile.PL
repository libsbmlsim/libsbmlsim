use ExtUtils::MakeMaker;
WriteMakefile(
    CC => 'gcc',
    INC => '-I/usr/local/include -I../../libsbmlsim',
    NAME => 'libsbmlsim',
    LIBS => ['-L/usr/local/lib -lsbml'],
    OBJECT => '*.o ../../../build/src/CMakeFiles/sbmlsim.dir/*.o'
);
