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

import biotransformer.btransformers.HGutBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

public class HGutBTTest {
	static HGutBTransformer hgbt;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		hgbt = new HGutBTransformer();
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
		IAtomContainer molecule 	= hgbt.getSmiParser().parseSmiles("Oc1ccc(cc1O)[C@H]3Oc2cc(O)cc(O)c2C[C@@H]3O");
		IAtomContainer metabolite 	= hgbt.getSmiParser().parseSmiles("OC1=CC=C(C=C1O)CC2OC(=O)CC2");
		
		ArrayList<Biotransformation> hgbbios = hgbt.simulateGutMicrobialMetabolism(
				molecule, true, true, 1, 0.5);
		IAtomContainerSet hgbmets = hgbt.extractProductsFromBiotransformations(hgbbios);
		
		assertTrue("There must be at least 1 biotranformation", (hgbbios != null && hgbbios.size()>0));
		
		assertTrue("5-[(3,4-dihydroxyphenyl)methyl]oxolan-2-one must be a metabolite.", 
				ChemStructureExplorer.atomContainerInclusionHolds(hgbmets, metabolite));
		
	}

}
