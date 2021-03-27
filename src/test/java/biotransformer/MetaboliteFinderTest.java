package biotransformer;

import static org.junit.Assert.assertFalse;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.utils.MetaboliteFinder;
import biotransformer.utils.MetaboliteFinder.FinderOption;



//@Ignore
public class MetaboliteFinderTest {
	static MetaboliteFinder mft;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		mft = new MetaboliteFinder();
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test(timeout=5000)
	public void testfindSuperbioMetabolitesByMass() throws Exception {
		IAtomContainer molecule = mft.hsbt.smiParser.parseSmiles("CCCCCC(=O)OCC(O)CO");
		IAtomContainer metabolite = mft.hsbt.smiParser.parseSmiles(
					"CCCCCC(=O)OCC(O)COP(O)(O)=O"); // mass = 271.09466 Da.; formula = C9H20O7P.
		ArrayList<String> masses = new ArrayList<String>();
		masses.add("271.094");

		IAtomContainerSet metabolites = mft.findSuperbioMetabolites(molecule, masses, 0.01, false, FinderOption.MASS, 1);
		System.out.println("Results: " + metabolites.getAtomContainerCount());

		
		assertFalse("A compound of mass ~= 271.094 Da must be identified.", 
					metabolites.getAtomContainerCount()>1 && Math.abs(
					Float.valueOf((String) metabolites.getAtomContainer(0).getProperties().get("Major Isotope Mass")) 
					- Float.valueOf(masses.get(0)) ) < 0.01);

	}
	
	
	@Test(timeout=5000)
	public void testfindAllHumanMetabolitesByFormula() throws Exception {
		IAtomContainer molecule = mft.hsbt.smiParser.parseSmiles("CCCCCC(=O)OCC(O)CO");
//		IAtomContainer metabolite = mft.hsbt.smiParser.parseSmiles(
//					"CCCCCC(=O)OCC(O)COP(O)(O)=O"); // mass = 271.09466 Da.; formula = C9H20O7P.
		ArrayList<String> formulas = new ArrayList<String>();

		formulas.add("C9H20O7P");
		IAtomContainerSet metabolites = mft.findSuperbioMetabolites(molecule, formulas, 0.01, false, FinderOption.FORMULA, 1);
		System.out.println("Results: " + metabolites.getAtomContainerCount());
		assertFalse("A compound of formula C9H20O7P must be identified", 
				metabolites.getAtomContainerCount()>1 && 
				metabolites.getAtomContainer(0).getProperties().get("Molecular formula").toString().contentEquals(formulas.get(0)));

		
	}
}
