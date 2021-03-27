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

import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

//@Ignore("This test class would be ignore but can later be run if the user includes EnviPath data.")
public class EnvMicroBTTest {
	static EnvMicroBTransformer envbt;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		envbt = new EnvMicroBTransformer();
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
		IAtomContainer molecule 	= envbt.getSmiParser().parseSmiles("CC(=O)Nc1ccc(cc1)C([O-])=O");
		IAtomContainer acetamide 	= envbt.getSmiParser().parseSmiles("CC(N)=O");
		
		ArrayList<Biotransformation> envbios = envbt.applyEnvMicrobialTransformationsChain(
				molecule, true, true, 1, 0.5);
		IAtomContainerSet envmets = envbt.extractProductsFromBiotransformations(envbios);
		
		assertTrue("There must be at least 1 biotranformation", (envbios != null && envbios.size()>0));
		
		assertTrue("acetamide must be an environmental microbial metabolite.", 
				ChemStructureExplorer.atomContainerInclusionHolds(envmets, acetamide));		
	}

}
