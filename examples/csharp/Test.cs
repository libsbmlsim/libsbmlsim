using System;

public class Test
{
  static void Main()
  {
    myResult result;
    result = Libsbmlsim.simulateSBMLFromFile("../sample.xml", 25.0, 0.01, 100, 1, Libsbmlsim.MTHD_RUNGE_KUTTA, 0);
    Libsbmlsim.write_result(result, "test.dat");
    Console.WriteLine("Simulation result is written to test.dat.");
  }
}
