import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.SBMLReader;
import org.sbml.libsbml.libsbml;

public class Test {
  public static void main(String[] args) {
    System.loadLibrary("sbmlsim");
    System.loadLibrary("sbmlj");

    SBMLDocument d = libsbml.readSBML("./src/MAPK.xml");
    String docstr = d.toSBML();
    myResult result = libsbmlsim.simulateSBMLFromString(docstr, 4000.0, 0.1, 100, 1, 41, 0);

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
        int index = i * numOfSp + j;
        double sp = result.getSpeciesValueAtIndex(index);
        System.out.print(" " + sp);
      }
      System.out.println();

      if (i == 10)
        break;
    }

    /*
    SWIGTYPE_p_p_char pp = result.getColumn_name_sp();
    int numOfSp = result.getNum_of_columns_sp();
    for (int i = 0; i < numOfSp; i++) {
      String ssp = libsbmlsim.stringArray_getitem(pp, i);
      System.out.println(ssp);
    }
    System.out.println(result.getColumn_name_time());
    */
  }
}
