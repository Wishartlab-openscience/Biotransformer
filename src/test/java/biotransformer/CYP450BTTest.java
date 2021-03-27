package biotransformer;

import static org.junit.Assert.assertFalse;
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
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

public class CYP450BTTest {
	static Cyp450BTransformer hCyp450;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		hCyp450 = new Cyp450BTransformer(BioSystemName.HUMAN);
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
	public void test() throws Exception{
		IAtomContainer molecule = hCyp450.getSmiParser().parseSmiles("CC(=O)Nc1ccccc1");
		IAtomContainer acetaminophen = hCyp450.getSmiParser().parseSmiles("CC(=O)NC1=CC=C(C=C1)O");
		IAtomContainer cyp450metabo_mode_2 = hCyp450.getSmiParser().parseSmiles("CC(=O)NC12C(C=CC=C1)O2"); 
		IAtomContainer _2_Hydroxy_N_phenylacetamide = hCyp450.getSmiParser().parseSmiles("C(C(=O)NC1=CC=CC=C1)O");
		
//		ArrayList<Biotransformation> cyp450bios = hCyp450.predictCyp450BiotransformationChain(
//				molecule, true, true, 1, 0.5);
		ArrayList<Biotransformation> cyp450bios_1 = hCyp450.predictCyp450BiotransformationsByMode(molecule, 1, true, true, 0);
		ArrayList<Biotransformation> cyp450bios_2 = hCyp450.predictCyp450BiotransformationsByMode(molecule, 2, true, true, 0);
		ArrayList<Biotransformation> cyp450bios_3 = hCyp450.predictCyp450BiotransformationsByMode(molecule, 3, true, true, 0);
		IAtomContainerSet cyp450mets_1 = hCyp450.extractProductsFromBiotransformations(cyp450bios_1);
		IAtomContainerSet cyp450mets_2 = hCyp450.extractProductsFromBiotransformations(cyp450bios_2);
		IAtomContainerSet cyp450mets_3 = hCyp450.extractProductsFromBiotransformations(cyp450bios_3);
		
//		System.out.println("cyp450mets (mode 1): " + cyp450mets_1.getAtomContainerCount());
//		System.out.println("cyp450mets (mode 2): " + cyp450mets_2.getAtomContainerCount());
//		System.out.println("cyp450mets (mode 3): " + cyp450mets_3.getAtomContainerCount());
		
		
		assertTrue("There must five unique metabolites predicted in mode 1: ", (cyp450mets_1 != null && cyp450mets_1.getAtomContainerCount()==5));
		assertTrue("There must four unique metabolites predicted in mode 2: ", (cyp450mets_2 != null && cyp450mets_2.getAtomContainerCount()==4));
		assertTrue("There must seven unique metabolites predicted in mode 3: ", (cyp450mets_3 != null && cyp450mets_3.getAtomContainerCount()==7));

		assertTrue("Acetaminophen must be a metabolite predicted in mode 1", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_1, acetaminophen));
		assertTrue("Acetaminophen must be a metabolite predicted in mode 2", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_2, acetaminophen));
		assertTrue("Acetaminophen must be a metabolite predicted in mode 3", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_3, acetaminophen));

		assertTrue("_2_Hydroxy_N_phenylacetamide (CC(=O)NC12C(C=CC=C1)O2) must be predicted in mode 1", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_1, _2_Hydroxy_N_phenylacetamide));
		assertFalse("_2_Hydroxy_N_phenylacetamide (CC(=O)NC12C(C=CC=C1)O2) must not be a metabolite predicted in mode 2", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_2, _2_Hydroxy_N_phenylacetamide));
		
		assertFalse("cyp450metabo_mode_2 (CC(=O)NC12C(C=CC=C1)O2) must not be predicted in mode 1", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_1, cyp450metabo_mode_2));
		assertTrue("cyp450metabo_mode_2 (CC(=O)NC12C(C=CC=C1)O2) must be a metabolite predicted in mode 2", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets_2, cyp450metabo_mode_2));
	
		
		
	}

}
