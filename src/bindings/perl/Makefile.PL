use ExtUtils::MakeMaker;
WriteMakefile(
    CC => 'gcc-mp-4.4',
    INC => '-I/usr/local/include',
    NAME => 'libsbmlsim',
    LIBS => ['-L/usr/local/lib -lsbml'],
    OBJECT => '*.o src/CMakeFiles/sbmlsim.dir/*.o'
);
