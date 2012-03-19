import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.SBMLReader;
import org.sbml.libsbml.libsbml;

public class Test {
  public static void main(String[] args) {
    System.loadLibrary("sbmlsim");
    System.loadLibrary("sbmlj");

    /*
    SBMLDocument d = libsbml.readSBML("./src/MAPK.xml");
    String docstr = d.toSBML();
    myResult result = libsbmlsim.simulateSBMLFromString(docstr, 4000.0, 0.1, 100, 1, 41, 0);
    */
    System.out.println("Simulate SBML model from File");
    myResult result = libsbmlsim.simulateSBMLFromFile("./src/MAPK.xml", 4000.0, 0.1, 100, 1, 41, 0);

    //libsbmlsim.print_result(result);
    //libsbmlsim.write_csv(result, "result.csv");
    //libsbmlsim.write_result(result, "result.dat");

    int numOfRows = result.getNumOfRows();
    System.out.println("numOfRows: " + numOfRows);

    int numOfSp = result.getNumOfSpecies();
    System.out.println("numOfSpecies: " + numOfSp);
    
    int numOfParam = result.getNumOfParameters();
    System.out.println("numOfParameters: " + numOfParam);

    int numOfComp = result.getNumOfCompartments();
    System.out.println("numOfCompartments: " + numOfComp);

    String timeName = result.getTimeName();
    System.out.println("TimeName: " + timeName);

    System.out.println("Species Name:");
    for (int i = 0; i < numOfSp; i++) {
      String sname = result.getSpeciesNameAtIndex(i);
      System.out.println("  " + sname);
    }
        
    System.out.println("Parameter Name:");
    for (int i = 0; i < numOfParam; i++) {
      String pname = result.getParameterNameAtIndex(i);
      System.out.println("  " + pname);
    }

    System.out.println("Compartment Name:");
    for (int i = 0; i < numOfComp; i++) {
      String cname = result.getCompartmentNameAtIndex(i);
      System.out.println("  " + cname);
    }

    System.out.println("Values:");
    for (int i = 0; i < numOfRows; i++) {
      double t = result.getTimeValueAtIndex(i);
      System.out.print(t);
      for (int j = 0; j < numOfSp; j++) {
        String sname = result.getSpeciesNameAtIndex(j);
        double val = result.getSpeciesValueAtIndex(sname, i);
        System.out.print(" " + val);
      }
      System.out.println();
      if (i == 10)
        break;
    }
  }
}
