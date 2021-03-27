package biotransformer;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.utils.ChemStructureExplorer;

public class ChemStructureTest {


	protected static IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
	protected static SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
	protected static SmilesParser	smiParser		= new SmilesParser(builder);
	public static IAtomContainerSet molecules		= DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	public static InChIGeneratorFactory inchiGenFactory;

	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		inchiGenFactory = InChIGeneratorFactory.getInstance();
		IAtomContainer molecule_1 	= smiParser.parseSmiles("CC(O)=Nc1ccc(O)cc1"); //4301564
		IAtomContainer molecule_2 	= smiParser.parseSmiles("CC(=O)Nc1ccc(cc1)C([O-])=O");
		IAtomContainer molecule_3 	= smiParser.parseSmiles("CNCC(=O)OC");
		IAtomContainer molecule_4 	= smiParser.parseSmiles("CNCC(OC)=O");
		IAtomContainer molecule_5 	= smiParser.parseSmiles("Oc1ccc(cc1O)[C@H]3Oc2cc(O)cc(O)c2C[C@@H]3O");
		IAtomContainer molecule_6 	= smiParser.parseSmiles("CCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CC\\C=C\\CC");
		
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_1).get("atomContainer"));
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_2).get("atomContainer"));
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_3).get("atomContainer"));
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_4).get("atomContainer"));
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_5).get("atomContainer"));
		molecules.addAtomContainer((IAtomContainer) ChemStructureExplorer.annotateWithProperties(molecule_6).get("atomContainer"));
		
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

//		System.out.println(molecules.getAtomContainer(1).getProperty("InChIKey").toString());
//		System.out.println(molecules.getAtomContainer(2).getProperty("InChIKey").toString());
		
		assertEquals("Formula must be C8H9NO2.", "C8H9NO2", 
				molecules.getAtomContainer(0).getProperty("Molecular formula").toString());
		
		assertEquals("Molecules 2 and 3 must have the same InChiKey.", true, ChemStructureExplorer.inchikeyEqualityHolds(molecules.getAtomContainer(2), 
				molecules.getAtomContainer(3)));

		assertEquals("atomContainerInclusionHolds() failed.", true, ChemStructureExplorer.atomContainerInclusionHolds(molecules,
				molecules.getAtomContainer(1)));
	
		assertEquals("This compound is a metabolizable polyphenol derivative.", true, 
				ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(molecules.getAtomContainer(4)));

		assertEquals("There must be 5 unique molecules.", 5, 
				ChemStructureExplorer.uniquefy(molecules).getAtomContainerCount());

	}

}
