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
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

public class ECBasedBTTest {
	
	static ECBasedBTransformer ecbt;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ecbt = new ECBasedBTransformer(BioSystemName.HUMAN);
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
		IAtomContainer molecule = ecbt.getSmiParser().parseSmiles("CCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CC\\C=C\\CC");
		IAtomContainer choline = ecbt.getSmiParser().parseSmiles("C[N+](C)(C)CCO");
		ArrayList<Biotransformation> ecbbios = ecbt.simulateECBasedMetabolismChain(
				molecule, true, true, 1, 0.5);
		IAtomContainerSet ecbmets = ecbt.extractProductsFromBiotransformations(ecbbios);		

		assertTrue("There must be at least 1 biotranformation", (ecbbios != null && ecbbios.size()>0));
		
		assertTrue("Choline must be a metabolite.", 
				ChemStructureExplorer.atomContainerInclusionHolds(ecbmets, choline));

	}

}
