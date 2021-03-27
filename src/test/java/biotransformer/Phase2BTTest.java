package biotransformer;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

public class Phase2BTTest {
	static Phase2BTransformer p2bt;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		p2bt = new Phase2BTransformer(BioSystemName.HUMAN);
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

	@Test
	public void test() throws Exception {
		IAtomContainer molecule 	= p2bt.getSmiParser().parseSmiles("Oc1ccc(cc1O)[C@H]3Oc2cc(O)cc(O)c2C[C@@H]3O");
		IAtomContainer metabolite 	= p2bt.getSmiParser().
				parseSmiles("O[C@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(OC2OC(C(O)C(O)C2O)C(O)=O)C(O)=C1");
		
		ArrayList<Biotransformation> p2bbios = p2bt.applyPhase2Transformations(
				molecule, true, true, true, 0.5);
		IAtomContainerSet p2bmets = p2bt.extractProductsFromBiotransformations(p2bbios);
		
		assertTrue("There must be at least 1 biotranformation", (p2bbios != null && p2bbios.size()>0));
		
		assertTrue("4'-O-glucuronide must be a metabolite.", 
				ChemStructureExplorer.atomContainerInclusionHolds(p2bmets, metabolite));
		
		
	}

}
