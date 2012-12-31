using System;

public class Test
{
  static void Main()
  {
    myResult result;
    result = libsbmlsim.simulateSBMLFromFile("../sample.xml", 25.0, 0.01, 100, 1, libsbmlsim.MTHD_RUNGE_KUTTA, 0);
    libsbmlsim.write_result(result, "test.dat");
    Console.WriteLine("Simulation result is written to test.dat.");
  }
}
