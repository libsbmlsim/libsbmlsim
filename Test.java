import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.SBMLReader;
import org.sbml.libsbml.libsbml;

public class Test {
  public static void main(String[] args) {
    System.loadLibrary("sbmlsim");
    System.loadLibrary("sbmlj");

    _myResult result = new _myResult();
    SBMLDocument d = libsbml.readSBML("./src/MAPK.xml");
    String docstr = d.toSBML();
    //System.out.println(docstr);
    result = libsbmlsim.simulateSBMLFromString(docstr, 4000.0, 0.1, 100, 1, 41, 0);
    int numOfRows = result.getNum_of_rows();
    System.out.println("numOfRows: " + numOfRows);
    doubleArray darr = doubleArray.frompointer(result.getValues_time());
    for (int i = 0; i < numOfRows; i++) {
      double da = darr.getitem(i);
      System.out.println(da);
    }
        
        
    /*
    SWIGTYPE_p_double dd = result.getValues_time();
    for (int i = 0; i < numOfRows; i++) {
      double val = libsbmlsim.doubleArray_getitem(dd, i);
      System.out.println(val);
    }
    */
    /*
    double[] dd = result.getValues_time();
    System.out.println(dd.length);
    */

    /*
    String[] ss = result.getColumn_name_sp();
    System.out.println(ss.length);
    System.out.println(ss[0]);
    */
    SWIGTYPE_p_p_char pp = result.getColumn_name_sp();
    int numOfSp = result.getNum_of_columns_sp();
    for (int i = 0; i < numOfSp; i++) {
      String ssp = libsbmlsim.stringArray_getitem(pp, i);
    System.out.println(ssp);
    }
    System.out.println(result.getColumn_name_time());
  }
}
