use ExtUtils::MakeMaker;
WriteMakefile(
  NAME => 'libsbmlsim',
  LIBS => ['-lsbml'],
  OBJECT => '*.o src/CMakeFiles/sbmlsim.dir/*.o'
);
