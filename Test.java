import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.SBMLReader;
import org.sbml.libsbml.libsbml;

public class Test {
  public static void main(String[] args) {
    _myResult result = new _myResult();
    System.loadLibrary("sbmlsim");
    System.loadLibrary("sbmlj");
    SBMLDocument d = libsbml.readSBML("./src/MAPK.xml");
    result = libsbmlsim.simulateSBMLFromString(d.toSBML(), 4000.0, 0.1, 100, 1, 41, 0);
  }
}
