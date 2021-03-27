/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.IOException;

/**
 * A class that implements functions to analyze the structure of a molecule
 * (e.g. through structure search). It also implements functions on collections
 * of molecules, such as the removal of duplicates.
 */

/**
 * @author Yannick Djoumbou
 *
 */

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
//import org.openscience.cdk.aromaticity.Aromaticity;
//import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
//import org.openscience.cdk.exception.InvalidSmilesException;
//import org.openscience.cdk.fragment.ExhaustiveFragmenter;
import org.openscience.cdk.graph.ConnectivityChecker;
//import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.isomorphism.Pattern;
//import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
//import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
//import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
//import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import com.google.common.base.Objects;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;
import biotransformer.fingerprint.ChemStructureFingerprinter;
//import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;

// import ambit2.structure2name.IUPACNameGenerator;


public class ChemStructureExplorer {

	protected static IChemObjectBuilder 	builder 		= SilentChemObjectBuilder.getInstance();
	protected static SmilesParser			smiParser		= new SmilesParser(builder);
	public static SmilesGenerator 			smiGen			= SmilesGenerator.isomeric();
	public static InChIGeneratorFactory inchiGenFactory;
	public static HBondAcceptorCountDescriptor hbaDCountDescriptor = new HBondAcceptorCountDescriptor();
	public static HBondDonorCountDescriptor hbdCountDescriptor 	= new HBondDonorCountDescriptor();
	public static RotatableBondsCountDescriptor rbCountDescriptor		= new RotatableBondsCountDescriptor();
	public static XLogPDescriptor xLogpDescriptor = new XLogPDescriptor();

	
	
	public ChemStructureExplorer() throws CDKException{
		ChemStructureExplorer.inchiGenFactory = InChIGeneratorFactory.getInstance();

	}
	
	public static String generateIsomericSmiles(IAtomContainer atc){	
		try{
			return smiGen.create(atc);
		}catch(Exception e){
			return null;
		}	
	}
	
	public static int ringCount(IAtomContainer atc){
		int ringCount = 0;
		RingSearch ringSearch = new RingSearch(atc);
		if(ringSearch != null){
			ringCount = ringSearch.numRings();
		}
		
		return ringCount;
	}
	
	public static LinkedHashMap<String, Object> createAtomContainerFromSmiles(String smiles, boolean standardize) throws Exception{
//		return smiParser.parseSmiles(smiles);			
		LinkedHashMap<String, Object> results = new LinkedHashMap<String, Object>();
		results.put("atomContainer", null);
		results.put("errors", new ArrayList<String>());
				
		try {
			IAtomContainer atc_0 = smiParser.parseSmiles(smiles);
			IAtomContainer atc = atc_0;
			if (standardize == true){
				atc = ChemStructureManipulator.standardizeMoleculeWithCopy(atc_0);
			}
			results.put("atomContainer", atc);
						
		}catch (CDKException c){
			results.put("errors", c.getMessage());
		}
		return results;	
	}

	public static  LinkedHashMap<String, Object>  createAtomContainerFromSmilesWithProperties(String smiles, boolean standardize) throws Exception{
		LinkedHashMap<String, Object> results = new LinkedHashMap<String, Object>();
		results.put("atomContainer", null);
		results.put("errors", new ArrayList<String>());
				
		try {
			
			IAtomContainer atc_0 = smiParser.parseSmiles(smiles);
			IAtomContainer atc = atc_0;
			if (standardize == true){
				atc = ChemStructureManipulator.standardizeMoleculeWithCopy(atc_0);
			}

			LinkedHashMap<String, String> physchemprops = computePhysicoChemicalProperties(atc);			
//			InChIGenerator gen0 = inchiGenFactory.getInChIGenerator(atc);
//			atc.setProperty("InChI", gen0.getInchi());
//			atc.setProperty("InChIKey", gen0.getInchiKey());			
			atc.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(atc));
			
			for(String p : physchemprops.keySet()) {
				atc.setProperty(p, physchemprops.get(p));
			}
//			atc.setProperty("Major Isotope Mass", physchemprops.get("Major Isotope Mass"));
//			atc.setProperty("ALogP", physchemprops.get("ALogP"));

			results.put("atomContainer", atc);
			
			
		}catch (CDKException c){
			results.put("errors", c.getMessage());
		}
		
		return results;			
	}

	
	public static  LinkedHashMap<String, Object>  annotateWithProperties(IAtomContainer atc) throws CDKException{
		LinkedHashMap<String, Object> results = new LinkedHashMap<String, Object>();
		results.put("atomContainer", null);
		results.put("errors", new ArrayList<String>());
				
		try {
			LinkedHashMap<String, String> physchemprops = computePhysicoChemicalProperties(atc);			
//			InChIGenerator gen0 = inchiGenFactory.getInChIGenerator(atc);
//			atc.setProperty("InChI", gen0.getInchi());
//			atc.setProperty("InChIKey", gen0.getInchiKey());
			atc.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(atc));
			
//			atc.setProperties(physchemprops);
			
			
			for(String p : physchemprops.keySet()) {
				atc.setProperty(p, physchemprops.get(p));
//				System.out.println(p + " : " + physchemprops.get(p));
			}
//			atc.setProperty("Major Isotope Mass", physchemprops.get("Major Isotope Mass"));
//			atc.setProperty("ALogP", physchemprops.get("ALogP"));

			results.put("atomContainer", atc);
//			System.out.println(atc);
//			System.out.println(atc.getProperties());
//			System.out.println(atc.getProperties());
			
			
		}catch (CDKException c){
			results.put("errors", c.getMessage());
		}
		
		return results;			
	}
	
	
	
	public static boolean isAromatic(IAtomContainer molecule){
		boolean aro = false;
		
		for(IBond b : molecule.bonds()){
			if(b.isAromatic()){
				aro = true;
				break;
			}
		}
		
		return aro;
	}
	
	
	/**
	 * Given a molecule A and a SMARTS expression S, list all occurrences of S
	 * in A. Each occurrence is represented as an ArrayList of atom indexes.
	 * 
	 * @param smartsString
	 * @param molecule
	 * @return A list of occurrences of S in A represented as an ArrayList of
	 *         atom indexes
	 * @throws SMARTSException
	 */
	public static List<List<Integer>> findAllOccurences(String smartsString,
		IAtomContainer molecule) throws SMARTSException {
		// might want to use IsomorphismTester as in SMIRKSManager
		List<List<Integer>> matches = new ArrayList<List<Integer>>();
		SmartsPatternCDK smartsPattern = new SmartsPatternCDK(smartsString);
		if (smartsPattern.hasSMARTSPattern(molecule) > 0) {
			matches = smartsPattern.getUniqueMatchingAtoms(molecule);
			IAtomContainer structure = smartsPattern.getMatchingStructure(molecule);

//			for (int k = 0; k < structure.getAtomCount(); k++) {
//				System.out.println(structure.getAtom(k));
//			}
//			System.out.println(matches);
		}

		return matches;
	}

	/**
	 * Given a molecule A and a MetabolicReaction M, find whether A matches all
	 * the structural constraints imposed by M
	 * 
	 * @param reaction
	 * @param molecule
	 * @return A boolean value that specifies whether or not A matches all the
	 *         structural constraints imposed by M
	 * @throws SMARTSException
	 * @throws CDKException
	 * @throws IOException 
	 */
	public static boolean compoundMatchesReactionConstraints(MetabolicReaction reaction,
			IAtomContainer molecule) throws SMARTSException, CDKException, IOException {
		boolean match = true;

		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		if (match) {
//			System.out.println("reaction.getReactantSMARTS().size() : " + reaction.getReactantSMARTS().size() );
//			System.out.println("reaction.getReactantSMARTS() == null? : " + reaction.getReactantSMARTS() == null );
//			System.out.println(reaction.name);
			
			for (int j = 0; j < reaction.getReactantSMARTS().size(); j++) {
//				System.out.println(reaction.getReactantSMARTS().get(j));
				Pattern smp = SmartsPattern.create(reaction.getReactantSMARTS().get(j), bldr);				
				boolean status = smp.matches(molecule);
//				SmartsPatternCDK smartsPattern = new SmartsPatternCDK(reaction.getReactantSMARTS().get(j));
//				boolean status = (smartsPattern.hasSMARTSPattern(molecule)>0);
//				SMARTSQueryTool querytool = new SMARTSQueryTool(reaction.getReactantSMARTS().get(j), SilentChemObjectBuilder.getInstance());
//				boolean status = querytool.matches(molecule);
				
				if (!status) {
					match = false;
					break;
				}
			}
		}
		if(match && (!(reaction.getExcludedReactantsSMARTS() == null || reaction.getExcludedReactantsSMARTS().isEmpty()))){
			for (int i = 0; i < reaction.getExcludedReactantsSMARTS().size(); i++) {
	
				Pattern smp2 = SmartsPattern.create(reaction.getExcludedReactantsSMARTS().get(i), bldr);
				boolean status2 = smp2.matches(molecule);
				if (status2) {
					match = false;
					break;
				}
			}
		}

		return match;

	}

	/**
	 * Given a molecule M and a SMARTS expression S, return all the
	 * substructures of M that match the expression S
	 * 
	 * @param smartsString
	 * @param molecule
	 * @return return all the substructures of M that match the expression S
	 * @throws SMARTSException
	 */
	public IAtomContainer findAllMatchingStructures(String smartsString,
		IAtomContainer molecule) throws SMARTSException {
		// might want to use IsomorphismTester as in SMIRKSManager
		SmartsPatternCDK smartsPattern = new SmartsPatternCDK(smartsString);
		IAtomContainer structure = smartsPattern.getMatchingStructure(molecule);
		return structure;
	}

	/**
	 * Given two molecules, find out whether they are equal. This function uses
	 * CDK's Pattern class to check for equality.
	 * 
	 * @param mol1
	 * @param mol2
	 * @return a boolean value that specifies whether the two molecules are
	 *         equal.
	 */
	/**
	 * WARNING: THIS METHOD SEEMS NOT TO BE SYMMERTRIC. WHICH IS WHY I HAD TO
	 * TEST BOTH DIRECTION. MIGHT BE NICE TO CHECK THE INCHIKEY EQUALLITY RATHER
	 */
	public static boolean equalityHolds(IAtomContainer mol1, IAtomContainer mol2) {
		Pattern mol1Pattern = Pattern.findIdentical(mol1);
		Pattern mol2Pattern = Pattern.findIdentical(mol2);
		return (mol1Pattern.matches(mol2) || mol2Pattern.matches(mol1));
	}

	/**
	 * Given two molecules, find out whether they are equal. This function uses
	 * InChI strings generated by CDK's InChIGeneratorFactory class as a
	 * comparator.
	 * 
	 * @param mol1
	 * @param mol2
	 * @return a boolean value that specifies whether the two molecules are
	 *         equal.
	 * @throws Exception 
	 */
	public static boolean inchikeyEqualityHolds(IAtomContainer mol1, IAtomContainer mol2)
			throws Exception {
		boolean equal = false;

		// Generate factory - throws CDKException if native code does not load
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		// Get InChIGenerator
		String inchikey1 = mol1.getProperty("InChIKey");
		String inchikey2 = mol2.getProperty("InChIKey");
		
		
		if(inchikey1 ==null){
			try {
				InChIGenerator gen1 = factory.getInChIGenerator(mol1);
				inchikey1 = gen1.getInchiKey();
			}catch(CDKException e){
				try {
					InChIGenerator gen1 = factory.getInChIGenerator(ChemStructureManipulator.preprocessContainer(mol1));
					inchikey1 = gen1.getInchiKey();
				} catch(NullPointerException n){
						System.err.println("Could not generate inchikey for molecule " + mol1.getProperty(CDKConstants.TITLE));
						return false;
					}
			}
		}

		if(inchikey2 ==null){
			try {
				InChIGenerator gen2 = factory.getInChIGenerator(mol2);
				inchikey2 = gen2.getInchiKey();
			}catch(CDKException e){
				try {
				InChIGenerator gen2 = factory.getInChIGenerator(ChemStructureManipulator.preprocessContainer(mol2));
				inchikey2 = gen2.getInchiKey();
				}	catch(NullPointerException n){
					System.err.println("Could not generate inchikey for molecule " + mol2.getProperty(CDKConstants.TITLE));
					return false;
				}
			}		
		}
		
		
		return inchikey1.equals(inchikey2);
	}

	/**
	 * Given a molecule (or AtomContainer) M and an AtomContainerSet S, check
	 * whether M is included in the set S.
	 * 
	 * @param molecules
	 * @param mol
	 * @return
	 * @throws Exception 
	 */
	public static boolean atomContainerInclusionHolds(IAtomContainerSet molecules,
			IAtomContainer mol) throws Exception {
		boolean inclusion = false;
		for (int i = 0; i < molecules.getAtomContainerCount(); i++) {

			if (inchikeyEqualityHolds(mol, molecules.getAtomContainer(i))) {
				inclusion = true;
				break;
			}
		}
		return inclusion;
	}

	/**
	 * Given an AtomContainerSet S, return an AtomContainerSet Su with only
	 * unique compounds of S
	 * 
	 * @param molecules
	 * @return A subset of the input AtomContainerSet S with only unique
	 *         compounds
	 * @throws Exception
	 */
	public static IAtomContainerSet uniquefy(IAtomContainerSet molecules)
			throws Exception {
		if (molecules != null && (!molecules.isEmpty()) && molecules.getAtomContainerCount() > 1) {
			
			IAtomContainerSet uniqueContainer = DefaultChemObjectBuilder.getInstance()
					.newInstance(IAtomContainerSet.class);
			
			uniqueContainer.addAtomContainer(molecules.getAtomContainer(0));
			
			for (int i = 1; i < molecules.getAtomContainerCount(); i++) {
				if (! ( (molecules.getAtomContainer(i) == null) || atomContainerInclusionHolds(uniqueContainer,
						molecules.getAtomContainer(i) ))) {
					uniqueContainer.addAtomContainer(molecules.getAtomContainer(i));
				}
			}

			return uniqueContainer;
		}

		else
			return molecules;

	}

	/**
	 * This functions looks at a molecule and returns each of its components
	 * (moieties) as separate molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return A set of molecules
	 * @throws CDKException
	 */
	public static IAtomContainerSet checkConnectivity(IAtomContainer molecule)
			throws CDKException {
		IAtomContainerSet partitions = new AtomContainerSet();
		IAtomContainerSet ms = ConnectivityChecker.partitionIntoMolecules(molecule);
//		System.out.println("MS: " + ms.getAtomContainerCount() + " molecule(s).");
	
		for (int k = 0; k < ms.getAtomContainerCount(); k++) {
			IAtomContainer current_metabolite = ms.getAtomContainer(k);
			AtomContainerManipulator
					.percieveAtomTypesAndConfigureAtoms(current_metabolite);
			for (IAtom atom : current_metabolite.atoms())
				if (atom.getImplicitHydrogenCount() == null) {
					atom.setImplicitHydrogenCount(0);
				}
			AtomContainerManipulator.suppressHydrogens(current_metabolite);
			partitions.addAtomContainer(current_metabolite);
		}
		return partitions;
	}
	
	public static boolean containsCarbon(IAtomContainer molecule) {
		boolean carbon = false;
		for(IAtom at : molecule.atoms()){
			if(at.getAtomicNumber() == 6){
				carbon = true;
				break;
			}
		}	
		return carbon;
	}
	
	
	public static boolean isPolyphenolOrDerivative(IAtomContainer molecule) throws SMARTSException, CDKException, CloneNotSupportedException{
		boolean polyphenol =  false;
		
		if(isMetabolizablePolyphenolOrDerivative(molecule)){
			// This must be removed and the definition of isPolyphenolOrDerivative should be improved. I noticed that
			// Urolithin B (OC1=CC2=C(OC(=O)C3=CC=CC=C23)C=C1) is not a polyphenolOr derivative but isMetabolizablePolyphenolOrDerivative.
			// Which is not very logical !!!
			polyphenol = true;
		} else {
			String constraintsValid = "["
					+ "$(O=[#6;R0](-[#6;R0]-,=[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // (dihydro)chalcone
					+ "$(O=[#6](-[#6;R0]-,=[#6;R0]-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)," // (dihydro)chalcone
					+ "$([H][#8;R0]-[#6;R0](-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#8][H])," // 
					+ "$([H][#8]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // 
					+ "$([#6;R0](=[#6;R0]/[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)," // stilbene
					+ "$([#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]-[#6;R1]=,:1[#6;R1](-[$(C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])=C([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // phenyl- and benzoic acids and conjugates
					+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[#1,#8])[#6;R1](-[#1,#8])=,:[#6;R1](-[#1,#8])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2])," // anthocyanidin
					+ "$([H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([H])=,:[#6;R1]1[C;R1]-,:1([H])-,:[#8][#6]=,:2[#6]=,:[#6](!@-[#8;X2])[#6]=,:[#6](!@-[#8;X2])[#6]=,:2-,:[C;R1]([H])([H])-,:[C;R1]-,:1([H])[#8])," // flavan-3-ol
					//+ "$([H][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]([H])=[#6;R1]-1[C;R1]1([H])[#8]-[#6]-2=[#6](-[#6](!@-[#8;X2])=[#6]-[#6](!@-[#8;X2])=[#6]-2)[C;R1]([H])([H])[C;R1]1([H])[#8])," // flavan-3-ol
					+ "$([H][#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1]([H])=,:[#6;R1]1[C;R1]1([H])[#8;R1]-[#6;R2]=,:2[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R2]=,:2-[#6;R1](=[O;X1])[C;R1]1([H])[#1,OX2H1,OX1-,$([#8]-[#6])])," // flavanone/flavanonol
					+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," //flavone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," // flavonol
					+ "$([O;X1]=[#6;R1]1[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#8;R1][#6;R1]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavone
					+ "$([H][C;R1]1([H])[#8;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6;R1](=[O;X1])[C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavanone
					+ "$([O;R0]=[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1](-[#6;R0]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-[#8;R1]-1)," // phenylvalerolactone
					+ "$([H][#6;R1]1=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]1[H])," // pyrrogallol or conjugates
					+ "$([H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:1[H])," // catechol and conjugates
					+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1-[#8;R0][H])," // catechol
					+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0][H])=,:[#6;R1]1-[#8;R0][H])," // pyrrogallol
					+ "$([H][#8;R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]=,:1[H])," // phloroglucinol
					+ "$([#8]-[#6;X3](=[O;R0])-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)," // benzoic acid
					+ "$([H][#6](=O)-[#6;R1]=,:1[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // benzaldehyde derivative and conjugates
					+ "$([#8;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R0]-,=[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-,=[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // curcominoid
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C2([H])[#8]-[#6](=O)-[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[#8][H])=,:[#6;R2]3-[#6;R2]3=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3-[#6](=O)-[#8]C12[H])-[#8][H])-[#8][H])," // ellagitannin pattern1
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#8])C([H])([#8]C([H])([#8][H])C1([H])[#8]-[#6]-2=O)C([H])([H])[#8])-[#8][H])," // ellagitannin pattern2
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C([H])([#8])C1([H])[#8])-[#8][H])," // ellagitannin pattern3
					+ "$([H][#8]!@-[#6]1=,:[#6][#6]2=,:[#6]3[#6]([#8][#6](=O)[#6]4=,:[#6][#6](!@-[#8][H])=,:[#6](!@-[#8][H])[#6]([#8][#6]2=O)=,:[#6]34)=,:[#6]1!@-[#8][H])," // ellagic acid
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:12)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]1)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2](=,:[#6;R1]1)[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:1[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]1=,:2)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]=,:2[#6;R1]=,:1)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]([#6;R1]=,:1)[#8;R1][#6;R1](=[O;R0])[#6;R2]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]21)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R1](=[O;R0])[#8;R1][#6;R2]1=,:2)," // Urolithin backbone
					+ "$([H][#8]-[#6;R0](=[O;R0])-[#6;R1]=,:1[#6;R2]=,:[#6;R2][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // fused carboxylated benzene
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2][#6;R2]=,:1)," // fused hydroxylated benzene
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1),"  // fused hydroxylated benzene
					+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylvaleric acid
					+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylbutyric acid
									
					// benzo[c]chromen‐6‐one (Urolithins must be hydroxylated at one or more positions)
					// Episin, J.C. (2013); Biological Significance of Urolithins, the Gut Microbial Ellagic Acid-Derived Metabolites: The Evidence So Far; 
					// Evidence-Based Complementary and Alternative Medicine; Volume 2013, Article ID 270418; 15 pages; http://dx.doi.org/10.1155/2013/270418
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([H])[#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:3[#6;R2]([#6;R2]=,:2[#6;R1]([H])=,:[#6;R1]1[H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]=,:3[H])"
	
					+ "]";
					// add depside
			
			SmartsPatternCDK smartsPatternValid = new SmartsPatternCDK(constraintsValid);
			String constraintsInvalid_1 ="[H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1";
			String constraintsInvalid_2 ="[H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1";		
			// Combination of 1 and 2
			String constraintsInvalid_3 ="[$([H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),$([H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)]";
					
			String athocyanidin = "[$([#8]-[#6;R1]-1=[#6;R1](-[#8;R1]-[#6]-2=[#6]-[#6](=O)-[#6]=[#6]-[#6]-2=[#6;R1]-1)-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8;A;H1X2][#6]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)[C;R1]1([#8;A;H1X2])[#8;R1]-[#6]-2=[#6]-[#6]([#8;A;H1X2])=[#6]-[#6]=[#6]-2-[#6;R1]=[#6;R1]1-[#8]),"
					+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1](-"
					+ "[R0;#1,$([O][H]),$([O]-[C]([H])[H])])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])]";
	
			SmartsPatternCDK flavonoidPattern 		=  new SmartsPatternCDK("[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6]-1-,=[#6;R1]-[#6;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#8;R1]-1),"
					+ "$(O=[#6]1[#6]=,:[#6]([#8][#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]12)-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)"
					+ "]");	                      
			
	
			SmartsPatternCDK otherFlavonoidsPattern = new SmartsPatternCDK("["
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]1=,:[#6][#6][#6]2=,:[#6]([#8]1)[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#8]-[#6]=,:2[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6](-[#8])[#6]=,:2-[#6]-[#6]1[H]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]=,:1[#6][#6]2=,:[#6]([#8][#6]=,:1[H])[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#6]-[#6]=,:2[#6](-[#8])=,:[#6]([H])[#6](-[#8])=,:[#6]([H])[#6]=,:2-[#8]C1([H])[H])"
			                         + "]");
			SmartsPatternCDK anthocyanidinPattern = new SmartsPatternCDK(athocyanidin);
			SmartsPatternCDK isoflavonoidPattern 	=  new SmartsPatternCDK("[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6;R1]-1-,=[#6]-[#8;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#6;R1]-1");
			SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
			
			if( smartsPatternValid.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || 
				 flavonoidPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0 ||
				 isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonePattern.hasSMARTSPattern(molecule)>0 || 			
					ChemStructureExplorer.findAllOccurences(constraintsInvalid_3,molecule).size()>0) {
					polyphenol = true;
			}
			
			
		}
			
		return polyphenol;
	}
	
	
	public static boolean isMetabolizablePolyphenolOrDerivative(IAtomContainer mol) throws SMARTSException, CDKException, CloneNotSupportedException{
		IAtomContainer molecule = mol.clone();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		boolean polyphenol =  false;
		String constraintsValid = "["
				+ "$(O=[#6;R0](-[#6;R0]-,=[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // (dihydro)chalcone
				+ "$(O=[#6](-[#6;R0]-,=[#6;R0]-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)," // (dihydro)chalcone
				+ "$([H][#8;R0]-[#6;R0](-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#8][H])," // 
				+ "$([H][#8]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // 
				+ "$([#6;R0](=[#6;R0]/[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)," // stilbene
				+ "$([#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]-[#6;R1]=,:1[#6;R1](-[$(C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])=C([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // phenyl- and benzoic acids and conjugates
				+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[#1,#8])[#6;R1](-[#1,#8])=,:[#6;R1](-[#1,#8])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2])," // anthocyanidin
				+ "$([H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([H])=,:[#6;R1]1[C;R1]-,:1([H])-,:[#8][#6]=,:2[#6]=,:[#6](!@-[#8;X2])[#6]=,:[#6](!@-[#8;X2])[#6]=,:2-,:[C;R1]([H])([H])-,:[C;R1]-,:1([H])[#8])," // flavan-3-ol
				//+ "$([H][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]([H])=[#6;R1]-1[C;R1]1([H])[#8]-[#6]-2=[#6](-[#6](!@-[#8;X2])=[#6]-[#6](!@-[#8;X2])=[#6]-2)[C;R1]([H])([H])[C;R1]1([H])[#8])," // flavan-3-ol
				+ "$([H][#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1]([H])=,:[#6;R1]1[C;R1]1([H])[#8;R1]-[#6;R2]=,:2[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R2]=,:2-[#6;R1](=[O;X1])[C;R1]1([H])[#1,OX2H1,OX1-,$([#8]-[#6])])," // flavanone/flavanonol
				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," //flavone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," // flavonol
				+ "$([O;X1]=[#6;R1]1[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#8;R1][#6;R1]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavone
				+ "$([H][C;R1]1([H])[#8;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6;R1](=[O;X1])[C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavanone
				+ "$([O;R0]=[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1](-[#6;R0]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-[#8;R1]-1)," // phenylvalerolactone
				+ "$([H][#6;R1]1=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]1[H])," // pyrrogallol or conjugates
				+ "$([H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:1[H])," // catechol and conjugates
				+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1-[#8;R0][H])," // catechol
				+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0][H])=,:[#6;R1]1-[#8;R0][H])," // pyrrogallol
				+ "$([H][#8;R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]=,:1[H])," // phloroglucinol
				+ "$([#8]-[#6;X3](=[O;R0])-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)," // benzoic acid
				+ "$([H][#6](=O)-[#6;R1]=,:1[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // benzaldehyde derivative and conjugates
				+ "$([#8;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R0]-,=[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-,=[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // curcominoid
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C2([H])[#8]-[#6](=O)-[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[#8][H])=,:[#6;R2]3-[#6;R2]3=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3-[#6](=O)-[#8]C12[H])-[#8][H])-[#8][H])," // ellagitannin pattern1
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#8])C([H])([#8]C([H])([#8][H])C1([H])[#8]-[#6]-2=O)C([H])([H])[#8])-[#8][H])," // ellagitannin pattern2
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C([H])([#8])C1([H])[#8])-[#8][H])," // ellagitannin pattern3
				+ "$([H][#8]!@-[#6]1=,:[#6][#6]2=,:[#6]3[#6]([#8][#6](=O)[#6]4=,:[#6][#6](!@-[#8][H])=,:[#6](!@-[#8][H])[#6]([#8][#6]2=O)=,:[#6]34)=,:[#6]1!@-[#8][H])," // ellagic acid
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:12)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]1)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2](=,:[#6;R1]1)[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:1[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]1=,:2)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]=,:2[#6;R1]=,:1)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]([#6;R1]=,:1)[#8;R1][#6;R1](=[O;R0])[#6;R2]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]21)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R1](=[O;R0])[#8;R1][#6;R2]1=,:2)," // Urolithin backbone
				+ "$([H][#8]-[#6;R0](=[O;R0])-[#6;R1]=,:1[#6;R2]=,:[#6;R2][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // fused carboxylated benzene
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2][#6;R2]=,:1)," // fused hydroxylated benzene
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1),"  // fused hydroxylated benzene
				+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylvaleric acid
				+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylbutyric acid
				
				// benzo[c]chromen‐6‐one (Urolithins must be hydroxylated at one or more positions)
				// Episin, J.C. (2013); Biological Significance of Urolithins, the Gut Microbial Ellagic Acid-Derived Metabolites: The Evidence So Far; 
				// Evidence-Based Complementary and Alternative Medicine; Volume 2013, Article ID 270418; 15 pages; http://dx.doi.org/10.1155/2013/270418
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([H])[#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:3[#6;R2]([#6;R2]=,:2[#6;R1]([H])=,:[#6;R1]1[H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]=,:3[H])"
				+ "]";
				// add depside
		
		SmartsPatternCDK smartsPatternValid = new SmartsPatternCDK(constraintsValid);
		
		String sp_phenols = "[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]-[#6;R1]=,:1[#6;R1](-[$(C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])=C([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]";
//		System.out.println("sp_phenols: " + ChemStructureExplorer.findAllOccurences(sp_phenols,molecule).size());
//		System.out.println("constraintsValid: " + ChemStructureExplorer.findAllOccurences(constraintsValid,molecule).size());
		
		/**
		 * (R1) Marín, L. et al. (2015); Bioavailability of Dietary Polyphenols and Gut Microbiota Metabolism: Antimicrobial Properties; Biomed Res Int. 2015; 2015: 905215.; doi:  10.1155/2015/905215
		 * (R2) Deprez, S. et al. (2001); “Transport of proanthocyanidin dimer, trimer, and polymer across monolayers of human intestinal epithelial Caco-2 cells,” Antioxidants and Redox Signaling, vol. 3, no. 6, pp. 957–967
		 * (R3) Monagas, M. et al. (2010); “Insights into the metabolism and microbial biotransformation of dietary  avan-3-ols and the bioactivity of their metabolites,” Food and Function, vol. 1, no. 3, pp. 233–253.
		 * 
		 * Oligomers with a degree of polymerization >3 are not absorbed in the small intestine, and therefore they are metabolized in the colon (R1-R3).
		 */
		String constraintsInvalid_1 ="[H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1";
		String constraintsInvalid_2 ="[H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1";		
		// Combination of 1 and 2
		String constraintsInvalid_3 ="[$([H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),$([H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)]";
				
		String athocyanidin = "[$([#8]-[#6;R1]-1=[#6;R1](-[#8;R1]-[#6]-2=[#6]-[#6](=O)-[#6]=[#6]-[#6]-2=[#6;R1]-1)-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8;A;H1X2][#6]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)[C;R1]1([#8;A;H1X2])[#8;R1]-[#6]-2=[#6]-[#6]([#8;A;H1X2])=[#6]-[#6]=[#6]-2-[#6;R1]=[#6;R1]1-[#8]),"
				+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1](-"
				+ "[R0;#1,$([O][H]),$([O]-[C]([H])[H])])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])]";
		
		SmartsPatternCDK otherFlavonoidsPattern = new SmartsPatternCDK("["
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]1=,:[#6][#6][#6]2=,:[#6]([#8]1)[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#8]-[#6]=,:2[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6](-[#8])[#6]=,:2-[#6]-[#6]1[H]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]=,:1[#6][#6]2=,:[#6]([#8][#6]=,:1[H])[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#6]-[#6]=,:2[#6](-[#8])=,:[#6]([H])[#6](-[#8])=,:[#6]([H])[#6]=,:2-[#8]C1([H])[H])"
		                         + "]");
		                      
		SmartsPatternCDK anthocyanidinPattern = new SmartsPatternCDK(athocyanidin);
		SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](-[#8])=O");

		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[*,#1;OX2H1,$(NC(=O)C)])-[#6](!@-[#8])-[#6]-1!@-[#8]),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[OX2H1,$(NC(=O)C)])-[#6]-[#6]-1!@-[#8])"
				+ "]"
				);
		
		SmartsPatternCDK oMethylPattern 		=  new SmartsPatternCDK("[#6;A;H3X4][#8;X2R0]-[#6;R1]");		
		SmartsPatternCDK flavonoidPattern 		=  new SmartsPatternCDK("[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6]-1-,=[#6;R1]-[#6;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#8;R1]-1),"
				+ "$(O=[#6]1[#6]=,:[#6]([#8][#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]12)-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)"
				+ "]");		
		
		SmartsPatternCDK isoflavonoidPattern 	=  new SmartsPatternCDK("[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6;R1]-1-,=[#6]-[#8;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#6;R1]-1");
		
		glycosylMoietyPattern.match(molecule);
		oMethylPattern.match(molecule);
		sulfatedRadicalPattern.match(molecule);
		
//		System.out.println((flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0));
//		System.out.println(glycosylMoietyPattern.getUniqueMatchingAtoms().size());
//		System.out.println(glycosylMoietyPattern.getUniqueMatchingAtoms().size());
//		System.err.println(smartsPatternValid.hasSMARTSPattern(molecule) > 0);
		
		
		if( (flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0)  
				&& (glycosylMoietyPattern.getUniqueMatchingAtoms().size() +
				sulfatedRadicalPattern.getUniqueMatchingAtoms().size() >= 2)){
			
			// I removed oMethylPattern.getUniqueMatchingAtoms().size() from the constraint  because some compounds have 2 or more methyl groups when
			// undergoing reduction.
			
			polyphenol = false;
//			System.err.println("(flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0)   && (glycosylMoietyPattern.getUniqueMatchingAtoms().size() + oMethylPattern.getUniqueMatchingAtoms().size() + sulfatedRadicalPattern.getUniqueMatchingAtoms().size() >= 2)");
		}
		
		
		else if ( (ChemStructureExplorer.findAllOccurences(constraintsInvalid_3,molecule).size() <= 3) && (smartsPatternValid.hasSMARTSPattern(molecule) > 0)) {
			polyphenol = true;
//			System.err.println("BLA");
		}
		
		return polyphenol;
	}

//	public static boolean isInvalidPhaseIICandidate(IAtomContainer molecule) throws SMARTSException{
//		boolean isInvalidPhaseIICandidate = false;
//		
//		if(isInvalidPhaseIIMetabolite(molecule)){
//			isInvalidPhaseIICandidate = false;
//		} else {
//			
//		}
//		
//		
//		
//		return isInvalidPhaseIICandidate;
//		
//	}
	
	
	public static boolean isPotentialPhaseIISubstrateByReactionPatternMatching(IAtomContainer substrate) throws Exception{
		boolean phaseII = false;
		
		if(containsCarbon(substrate)){
			IAtomContainer subs = substrate.clone();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(subs);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(subs);
			IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
			LinkedHashMap<String, String> sfeatures = ChemStructureFingerprinter.getPhaseIIfingerpint();
			
			
//			if(isInvalidPhaseIISubstrateByExlcusion(subs)){
//				
//			}else{
				for(String s : sfeatures.values()){
					SmartsPattern pattern = SmartsPattern.create(s, builder);
					if(pattern.matches(subs)){
						phaseII = true;
						break;
					}			
				}
//			}
			

		}
		
		return phaseII;
	}
	
	public static boolean isInvalidPhaseIISubstrateByExlcusion(IAtomContainer substrate) throws Exception{
		return (ChemicalClassFinder.isEtherLipid(substrate) || ChemicalClassFinder.isGlyceroLipid(substrate) || 
				ChemicalClassFinder.isGlycerophosphoLipid(substrate) || ChemicalClassFinder.isSphingoLipid(substrate)
				||ChemicalClassFinder.isAcylCoAConjugate(substrate) || ChemicalClassFinder.isGlutathioneConjugate(substrate) || 
				ChemicalClassFinder.isOligoOrPolysaccharide(substrate) || ChemicalClassFinder.isTetrapyrrole(substrate)
				);
	}

	public static boolean isInvalidPhase2Metabolite(IAtomContainer mol) throws SMARTSException, CDKException, CloneNotSupportedException{
		
		// This applies to humans
		
		//		R1) Mol Pharm. 2012 Apr 2; 9(4): 862–873. 10.1021/mp200400s
		//		R2) Curr Drug Metab. 2011 November ; 12(9): 900–916 (Regioselective Sulfation and Glucuronidation of Phenolics: Insights into the Structural Basis of Conjugation)
		//		R3) J. Agric. Food Chem. 2009, 57, 10134–10142 (DOI:10.1021/jf901450z)
		//		R4) Free Radicals Research, Vol. 35, pp. 941-952 (http://dx.doi.org/10.1080/10715760100301441)
		//		R5) Biol. Chem., Vol. 386, pp. 279–283, March 2005 (DOI 10.1515/BC.2005.033)
		//		R6) Regioselective Monosulfation and Disulfation of the Phytoestrogens Daidzein and Genistein by Human Liver Sulfotransferases
		//		R7) JOURNAL OF FUNCTIONAL FOODS 7 (2014) 43–53 (Absorption and metabolism of proanthocyanidins)
		//		R8) J. Agric. Food Chem., 2012, 60 (14), pp 3592–3598 (Characterization of Sulfated Quercetin and Epicatechin Metabolites)
		//		R9) Br J Nutr. 2015 Feb 14;113(3):454-63. doi: 10.1017/S0007114514003511. Epub 2015 Jan 9. (Phenolic sulfates as new and highly abundant metabolites in human plasma after ingestion of a mixed berry fruit purée)
		//		Personal communication & expert panel validation
		
		IAtomContainer molecule = mol.clone();
//		molecule = ChemStructureManipulator.preprocessContainer(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		SmilesGenerator s =  new SmilesGenerator();
		
//		System.out.println("STRUCTURE TO VALIDATE: " + s.create(molecule));
		
		boolean bad =  false;
		
		String constraints = "["
				// Glycylglucuronide and glucuronidylglycine conjugates
				+ "$([H]N(C)CC(=O)OC1OC(C(O)C(O)C1O)C(O)=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#7]-[#6]-[#6]([#8,#7;A])=O),"
				
				// glutathionylglycine conjugates
				+ "$(CNC(CCC(=O)NCC(O)=O)C(O)=O),"
								
				// glutamatylglucuronide and glucuronidylglutamate conjugates
				+ "$(CNC(CCC(=O)OC1OC(C(O)C(O)C1O)C(O)=O)C(O)=O),"
				+ "$([H]N(C(CCC(O)=O)C(O)=O)C(=O)C1OC(O)C(O)C(O)C1O),"
				
				// bad flavan-3-ol metabolites
				
				// 3-O-sulfate
				+ "$([H][#8]S(=O)(=O)[#8]-[#6]-,:1-,:C([H])([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				
				//6-O-Sulfate
				+ "$([H][#8]-[#6]-,:1-,:C([H])([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#8]S(=O)(=O)[#8][H]),"
				// 8-O-Sulfate
				+ "$([H][#8]-[#6]-,:1-,:C([H])([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8]S(=O)(=O)[#8][H])[#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
//				// 4'-O-Glucuronide
//				+ "$([H]C-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8]-,:C([H])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:2)-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-,:[#6]-,:1[#8;A;H1X2]),"
//				// 5-O-Glucuronide
//				+ "$([H]C-,:1([H])-,:[#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)[#8]-,:C([H])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[#6]-,:1[#8;A;H1X2]),"
				// 6-O-Glucuronide
				+ "$([H]C-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R2]=,:2[#8]-,:C([H])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[#6]-,:1[#8;A;H1X2])-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				 // 8-O-Glucuronide
				+ "$([H]C-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8]-[#6]-3-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-3-[#8])-[#6](-[#8])=O)[#6;R2]=,:2[#8]-,:C([H])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[#6]-,:1[#8;A;H1X2]),"
							
			    // bad flavqnone metabolites
				// 5-O-Sulfonate
				+ "$([H]C-,:1([H])-,:[#6](=O)[#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]S([#8])(=O)=O)[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				// 6-O-Sulfonate
				+ "$([H]C-,:1([H])-,:[#6](=O)[#6;R2]=,:2[#6;R1]=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#8]S([#8])(=O)=O),"
				// 8-O-Sulfonate
				+ "$([H]C-,:1([H])-,:[#6](=O)[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8]S([#8])(=O)=O)[#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				// 8-O-Glucuronide
				+ "$([H]C-,:1([H])-,:[#6](=O)[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8]-[#6]-3-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-3-[#8])-[#6](-[#8])=O)[#6;R2]=,:2[#8]-,:C-,:1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"

				// bad flavone metabolites
				// 4'-O-Sulfate
				+ "$([H][#6]1=,:[#6]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8]S([#8])(=O)=O),"
				// 5-O-Sulfate
				+ "$([H][#6]1=,:[#6]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8]S([#8])(=O)=O)[#6;R2]=,:2[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				// 6-O-Sulfate
				+ "$([H][#6]1=,:[#6]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R2]=,:2[#6]1=O)-[#8]S([#8])(=O)=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				// 8-O-Sulfate
				+ "$([H][#6]1=,:[#6]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]S([#8])(=O)=O)[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				// 8-O-Glucuronide
				+ "$([H][#6]1=,:[#6]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"			
				
				 // bad flavonol metabolites				
				 // 4'-O-Sulfate
				+ "$([#8;X2]-[#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]),"
				 // 6-O-Sulfate
				+ "$([H][#8]-[#6]1=,:[#6]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R2]=,:2[#6]1=O)-[#8]S([#8;A;X2H1,X1-])(=O)=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				 // 8-O-Sulfate
				+ "$([H][#8]-[#6]1=,:[#6]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]S([#8;A;X2H1,X1-])(=O)=O)[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"
				 // 8-O-Glucuronide
				+ "$([H][#8]-[#6]1=,:[#6]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)[#6]1=O)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),"

				// bad isoflavone metabolites
				// 2-O-Sulfonate
				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]=,:1[#8;R1][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R1](=[O;X1])[#6;R1]=,:1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," 
				// 5-O-Sulfonate
				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:3)[#6;R1](=[O;X1])[#6;R2]1=,:2)," 
				// 6-O-Sulfonate
				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:3)[#6;R1](=[O;X1])[#6;R2]=,:2[#6;R1]=,:1)," 
				// 8-O-Sulfonate
				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]1[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)[#6;R1]2=[O;X1])," 
				// 2-O-Glucuronide
				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6]=,:2[#8;R1][#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#6;R1](=[O;X1])[#6;R1]=,:2-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," 
				// 3'-O-Glucuronide
				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:2)-[#6;R1]2=,:[#6][#8;R1][#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#6;R1]2=[O;X1])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," 
				// 6-O-Glucuronide
				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6]=,:[#6;R1](-[#6;R1]=,:4[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:4)[#6;R1](=[O;X1])[#6;R2]=,:3[#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H]),"
				// 8-O-Glucuronide
				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]3=,:[#6;R2]2[#8;R1][#6]=,:[#6;R1](-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)[#6;R1]3=[O;X1])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H]),"  
				
				
				
				// phloroglucinol glucuronide, sulfate, methyl conjugates
				+ "$([H][#6;R1]=,:1[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]([H])[#6;R1](-[#8;R0]S([#8;A;X2H1,X1-])(=O)=O)=,:[#6;R1]([H])[#6;R1]=,:1-[#8;R0]-[#1,#6,#16])," // sulfate
				+ "$([H][#6;R1]=,:1[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]([H])[#6;R1](-[#8;R0]C([H])([H])[H])=,:[#6;R1]([H])[#6;R1]=,:1-[#8;R0]-[#1,#6,#16])," // methyl
				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8;R0]-[#6;R1]=,:2[#6;R1]([H])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]([H])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:2[H])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // glucuronide				
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]=,:[#6]-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)," // 1,2-diglucuronide
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8][#6;A]~[#6;A][#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)," // 1,2-diglucuronide
				
				// Other bad patterns
				// glycosylsulfate
				+"$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8,#7,#16;A])-[#6;R1]([#8,#7,#16;A])-[#6;R1]-1-[#8]S([#8])(=O)=O),"
				// glycosylsulfate
				+ "$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8,#7,#16;A])-[#6;R1](!@-[#8]S([#8])(=O)=O)-[#6;R1]-1!@-[#8,#7,#16;A]),"
				// glycosylsulfate
				+ "$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](-[#8]S([#8])(=O)=O)-[#6;R1]([#8,#7,#16;A])-[#6;R1]-1!@-[#8,#7,#16;A])" 				
				
				+ "]";
				
		
//		String constraints = "["
//				// Glycylglucuronide and glucuronidylglycine conjugates
//				+ "$([H]N(C)CC(=O)OC1OC(C(O)C(O)C1O)C(O)=O),"
//				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#7]-[#6]-[#6]([#8,#7;A])=O),"
//				
//				// glutathionylglycine conjugates
//				+ "$(CNC(CCC(=O)NCC(O)=O)C(O)=O),"
//								
//				// glutamatylglucuronide and glucuronidylglutamate conjugates
//				+ "$(CNC(CCC(=O)OC1OC(C(O)C(O)C1O)C(O)=O)C(O)=O),"
//				+ "$([H]N(C(CCC(O)=O)C(O)=O)C(=O)C1OC(O)C(O)C(O)C1O),"
//				
//				// bad flavan-3-ol metabolites
//				+ "$([H][C;R1]-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1]-,:[C;R1]([H;R0])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[C;R1]-,:1([H])[#8;X2R0][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 3-O-sulfate
//				+ "$([H][C;R1]-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1]-,:[C;R1]([H;R0])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[C;R1]-,:1([H;R0])[#8;X2])-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 6-O-sulfate
//				+ "$([H][C;R1]-,:1([H])-,:[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6;R2]=,:2[#8;R1]-,:[C;R1]([H;R0])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-,:[C;R1]-,:1([H;R0])[#8;X2]),"  // 8-O-sulfate
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:2)[C;R1]-,:2([H;R0])-,:[#8;R1][#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3-,:[C;R1]([H])([H])-,:[C;R1]-,:2([H;R0])[#8;X2])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 4'-O-Glucuronide
//				+ "$([H][#8]-[#6]1C([H])([#8])C([H])([#8]C([H])([#8]-[#6;R1]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1]-,:[C;R1]([H;R0])([#6;R1]=,:4[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:4)-,:[C;R1]([H;R0])([#8;X2][H])-,:[C;R1]([H])([H])-,:[#6;R2]2=,:3)C1([H])[#8][H])[#6](-[#8])=O)," // 5-O-Glucuronide
//				+ "$([H;R0][C;R1]1([#8;R1]-[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8;R0]-[#6]-3-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-3-[#8])-[#6](-[#8])=O)[#6;R2]=,:2[#6;A;H2X4R1][#6;R1]1[#8;A;H1X2])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 5-O-Glucuronide
//
////				+ "$([H;R0][C;R1]1([#8;R1]-[#6;R2]-2=[#6;R2]([#6;A;H2X4R1][#6;R1]1[#8;A;H1X2])-[#6;R1](-[#8;R0]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O)=[#6;R1]-[#6;R1]=[#6;R1]-2)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // 5-O-Glucuronide
//				
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1]-,:[C;R1]([H;R0])([#6;R1]=,:4[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:4)-,:[C;R1]([H;R0])([#8;X2])-,:[C;R1]([H])([H])-,:[#6;R2]=,:3[#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 6-O-Glucuronide
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R2]=,:3-[#8;R1][C;R1]([H;R0])([#6;R1]=,:4[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:4)[C;R1]([H;R0])([#8;X2])[C;R1]([H])([H])[#6;R2]=,:3[#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 6-O-Glucuronide
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]3=,:[#6;R2]2-[#8;R1][C;R1]([H;R0])([#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)[C;R1]([H;R0])([#8;X2][H])[C;R1]3([H])[H])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H]),"  // 8-O-Glucuronide
//	
//				// bad flavqnone metabolites				
//				+ "$([H][C;R1]1([H])[#6;R1](=[O;X1])-[#6]=,:2[#6](=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 5-O-Sulfonate
//				+ "$([H][C;R1]1([H])[#6;R1](=[O;X1])-[#6]=,:2[#6]=,:[#6]([#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 6-O-Sulfonate
//				+ "$([H][C;R1]1([H])[#6;R1](=[O;X1])-[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 8-O-Sulfonate
//
//				// bad flavonol metabolites				
//				+ "$([#8;X2]-[#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 4'-O-Sulfate
//				+ "$([#8;X2]-[#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R2]=,:2[#6;R1]1=[O;X1])-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 6-O-Sulfate
//				+ "$([#8;X2]-[#6;R1]1=,:[#6;R1]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 8-O-Sulfate
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6]2=,:[#6][#6]=,:[#6][#6]3=,:[#6]2[#8][#6;R1](-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)=,:[#6;R1](-[#8;X2])[#6;R1]3=[O;X1])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 8-O-Glucuronide
//				
//				// bad flavone metabolites
//				// 
//				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])," // 4'-O-Sulfate
//				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6;R2]=,:2[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 5-O-Sulfate
//				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R2]=,:2[#6;R1]1=[O;X1])-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 6-O-Sulfate
//				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6;R2]2=,:[#6;R2]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]2-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6;R1]1=[O;X1])-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 8-O-Sulfate
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]3=,:[#6;R2]2[#8][#6;R1](=,:[#6;R1]([H])[#6;R1]3=[O;X1])-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 8-O-Glucuronide
//							
//				// bad isoflavone metabolites
//				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]=,:1[#8;R1][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R1](=[O;X1])[#6;R1]=,:1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // 2-O-Sulfonate
//				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:3)[#6;R1](=[O;X1])[#6;R2]1=,:2)," // 5-O-Sulfonate
//				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:3)[#6;R1](=[O;X1])[#6;R2]=,:2[#6;R1]=,:1)," // 6-O-Sulfonate
//				+ "$([#8;A;X2H1,X1-][S;X4](=[O;X1])(=[O;X1])[#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]1[#8;R1][#6;R1]=,:[#6;R1](-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)[#6;R1]2=[O;X1])," // 8-O-Sulfonate
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6]=,:2[#8;R1][#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#6;R1](=[O;X1])[#6;R1]=,:2-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 2-O-Glucuronide
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:2)-[#6;R1]2=,:[#6][#8;R1][#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#6;R1]2=[O;X1])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 3'-O-Glucuronide
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6]=,:[#6;R1](-[#6;R1]=,:4[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:4)[#6;R1](=[O;X1])[#6;R2]=,:3[#6;R1]=,:2)C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // 6-O-Glucuronide
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8]-[#6;R1]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]3=,:[#6;R2]2[#8;R1][#6]=,:[#6;R1](-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)[#6;R1]3=[O;X1])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H]),"  // 8-O-Glucuronide
//				
//				// phloroglucinol glucuronide, sulfate, methyl conjugates
//				+ "$([H][#6;R1]=,:1[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]([H])[#6;R1](-[#8;R0]S([#8;A;X2H1,X1-])(=O)=O)=,:[#6;R1]([H])[#6;R1]=,:1-[#8;R0]-[#1,#6,#16])," // sulfate
//				+ "$([H][#6;R1]=,:1[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]([H])[#6;R1](-[#8;R0]C([H])([H])[H])=,:[#6;R1]([H])[#6;R1]=,:1-[#8;R0]-[#1,#6,#16])," // methyl
//				+ "$([H][#8]-[#6](=O)C1([H])[#8]C([H])([#8;R0]-[#6;R1]=,:2[#6;R1]([H])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]([H])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:2[H])C([H])([#8][H])C([H])([#8][H])C1([H])[#8][H])," // glucuronide				
//				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]=,:[#6]-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)," // 1,2-diglucuronide
//				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8][#6;A]~[#6;A][#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)," // 1,2-diglucuronide
//				
//				// Other bad patterns
//
//				+"$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8,#7,#16;A])-[#6;R1]([#8,#7,#16;A])-[#6;R1]-1-[#8]S([#8])(=O)=O)," // glycosylsulfate
//				+ "$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8,#7,#16;A])-[#6;R1](!@-[#8]S([#8])(=O)=O)-[#6;R1]-1!@-[#8,#7,#16;A])," // glycosylsulfate
//				+ "$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](-[#8]S([#8])(=O)=O)-[#6;R1]([#8,#7,#16;A])-[#6;R1]-1!@-[#8,#7,#16;A])" // glycosylsulfate
//				+ "]";
//				
				
				// benzo[c]chromen‐6‐one (Urolithins must be hydroxylated at one or more positions)
				// Episin, J.C. (2013); Biological Significance of Urolithins, the Gut Microbial Ellagic Acid-Derived Metabolites: The Evidence So Far; 
				// Evidence-Based Complementary and Alternative Medicine; Volume 2013, Article ID 270418; 15 pages; http://dx.doi.org/10.1155/2013/270418
//				+ "$([H][#6;R1]1=,:[#6;R2]2[#6;R2]3=,:[#6;R2]([#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1[H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]3[H]),"
				
				
		SmartsPatternCDK _6H_benzo_c_chromen_6_one =  new SmartsPatternCDK("[H][#6;R1]1=,:[#6;R2]2[#6;R2]3=,:[#6;R2]([#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1[H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]3[H]");
				// + "[H][#6;R1]-1=[#6;R1]([H])-[#6;R2]-2=[#6;R2](-[#6;R1]([H])=[#6;R1]-1[H])-[#6;R2]-1=[#6;R1]([H])-[#6;R1]([H])=[#6;R1]([H])-[#6;R1]([H])=[#6;R2]-1-[#8;R1]-[#6;R1]-2=[O;R0]");
		
		String athocyanidin = "[$([#8]-[#6;R1]-1=[#6;R1](-[#8;R1]-[#6]-2=[#6]-[#6](=O)-[#6]=[#6]-[#6]-2=[#6;R1]-1)-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8;A;H1X2][#6]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)[C;R1]1([#8;A;H1X2])[#8;R1]-[#6]-2=[#6]-[#6]([#8;A;H1X2])=[#6]-[#6]=[#6]-2-[#6;R1]=[#6;R1]1-[#8]),"
				+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1](-"
				+ "[R0;#1,$([O][H]),$([O]-[C]([H])[H])])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])]";
		
		SmartsPatternCDK otherFlavonoidsPattern = new SmartsPatternCDK("["
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]1=,:[#6][#6][#6]2=,:[#6]([#8]1)[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#8]-[#6]=,:2[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6](-[#8])[#6]=,:2-[#6]-[#6]1[H]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]=,:1[#6][#6]2=,:[#6]([#8][#6]=,:1[H])[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#6]-[#6]=,:2[#6](-[#8])=,:[#6]([H])[#6](-[#8])=,:[#6]([H])[#6]=,:2-[#8]C1([H])[H])"
		                         + "]");
		                      
		
		
		/**
		 * Personal communication (Claudine Manach)
		 * For polyphenols with a catechol group, the observed phase II biotransformations are:
		 * 1) Sulfonation or glucuronidation of the 3’ (major site of biotransformation) or 4’-OH, which is frequent. 
		 * In most cases the metabolite produced circulates and is excreted as such. A second glucuronidation can occur, 
		 * but less frequently, on another non-adjacent OH of the molecule leading to sulfo-glucuronides or diglucuronides. 
		 * Disulfonation is not reported for flavonoids. Combination of 3 glucuronidation / sulfonation does not occur.
		 * 
		 * 2) Another frequent biotransformation is the methylation of 3’ or 4’-OH by COMT. Methylation occurs on one or 
		 * the other of the OH groups of a catechol, but not on both at the same time.  Methylation is frequently combined 
		 * with a glucuronidation or a sulfonation that occurs on another OH of the molecule.

		 */
		
		
		SmartsPatternCDK adjacentSulfateGlucuronideSubstitutionsOnCatechol = new SmartsPatternCDK("["
				+ "$([H][#8]S(=O)(=O)[#8]-[#6;R1]-1=[#6;R1](-[#8;R0]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6;R1]=[#6;R1]([H])-[#6;R1](-[#6])=[#6;R1]-1[H]),"
				+ "$([H][#6;R1]1=,:[#6;R1][#6;R1](-[#8;R0]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)=,:[#6;R1](-[#8]S([#8;A;X2H1,X1-])(=O)=O)[#6;R1]([H])=,:[#6;R1]1-[#6]),"
				+ "$([H][#8]S(=O)(=O)[#8]-[#6;R1]-1=[#6;R1](-[#8;R0]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6;R1]([H])=[#6;R1](-[#6])-[#6;R1]([H])=[#6;R1]-1),"
				+ "$([H][#6;R1]1=,:[#6;R1][#6;R1](-[#8]S([#8;A;X2H1,X1-])(=O)=O)=,:[#6;R1](-[#8;R0]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)[#6;R1]([H])=,:[#6;R1]1-[#6]),"
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]~[#6]-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)"
				+ "]");

		
		SmartsPatternCDK _2p_glucuronidatedCatechol = new SmartsPatternCDK(
				"[$([#6]-[#6;R1]-1=[#6;R1](-[*,#1;!$([OX2H1])])-[#6;R1](-[*,#1;!$([OX2H1])])=[#6;R1](-[*,#1;!$([OX2H1])])-[#6;R1](-[#8;R0])=[#6;R1]-1-[#8]-[#6](=O)-[#6;R1]-1-[#8;R1]-[#6;R1]([#8;A;H1X2])-[#6;R1]([#8;A;H1X2])-[#6;R1]([#8;A;H1X2])-[#6;R1]-1[#8;A;H1X2].[#6]-[#6;R1]-1=[#6;R1](-[#1,!#6])-[#6;R1](-[#1,!#6])=[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])=[#6;R1]-1-[#8]-[#6](=O)-[#6]-1-[#8]-[#6]([#8;A;H1X2])-[#6]([#8;A;H1X2])-[#6]([#8;A;H1X2])-[#6]-1[#8;A;H1X2])"
				+ "]");

		SmartsPatternCDK _2p_sulfatedCatechol = new SmartsPatternCDK(
				 "[$([H][#8]S(=O)(=O)[#8]-[#6;R1]-1=[#6;R1](-[#8;R0])-[#6;R1](-[#1,OX2H1,$([O]-[CX4H3])])=[#6;R1](-[#1,OX2H1,$([O]-[CX4H3])])-[#6;R1](-[#1,OX2H1,$([O]-[CX4H3])])=[#6;R1]-1-[#6])"
				+ "]");
			
		SmartsPatternCDK smartsPattern = new SmartsPatternCDK(constraints);
		SmartsPatternCDK anthocyanidinPattern = new SmartsPatternCDK(athocyanidin);
		SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](-[#8])=O");
		SmartsPatternCDK _12_ONS_diglucuronidatedPattern = new SmartsPatternCDK("[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](=O)[#8,#7,#16;A][#6]~[#6][#8,#7,#16;A][#6](=O)-[#6;R1]-1-[#8;R1]-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0]");

		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#8]-[#6;R0]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8])-[#6;R1]-1!@-[#8]),"
				+ "$([#6]!@-[#6;R1]-1-[#8;R1]-[#6;R1](!@-[#8])-[#6;R1](!@-[#8])-[#6;R1](!@-[#8])-[#6;R1]-1!@-[#8])"
				+ "]"
				);
		
		SmartsPatternCDK oMethylPattern 		=  new SmartsPatternCDK("[#6;A;H3X4][#8;X2R0]-[#6;R1]");		
		SmartsPatternCDK flavonoidPattern 		=  new SmartsPatternCDK("[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6]-1-,=[#6;R1]-[#6;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#8;R1]-1),"
				+ "$(O=[#6]1[#6]=,:[#6]([#8][#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]12)-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)"
				+ "]");
		SmartsPatternCDK isoflavonoidPattern 	=  new SmartsPatternCDK("[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6;R1]-1-,=[#6]-[#8;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#6;R1]-1");
		
		
		// sulfated glucuronide moieties (where sulfate is attached to glucuronide) and glucuronidated glucuronide (two or more glucuronides linked together through an ether bond)
		SmartsPatternCDK linkedPolyGlucuronidatedPattern = new SmartsPatternCDK("[$([#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0]-[#6;R1]-2-[#8]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-2-[#8;R0])-[#6](-[#8;R0])=O)-[#6;R1]-1-[#8;R0])-[#6;R0](-[#8;R0])=[O;R0]),"
				+ "$([#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0]-[#6;R1]-2-[#8]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-2-[#8;R0])-[#6](-[#8;R0])=O)-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6;R0](-[#8;R0])=[O;R0]),"
				+ "$([#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0]-[#6;R1]-1-[#8]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](-[#8;R0])=O)-[#6;R0](-[#8;R0])=[O;R0])"
				+ "$(OC1OC(C(OC2OC(C(O)C(O)C2O)C(O)=O)C(O)C1O)C(O)=O),$(OC1OC(C(O)C(OC2OC(C(O)C(O)C2O)C(O)=O)C1O)C(O)=O),"
				+ "$(OC1OC(C(O)C(O)C1OC1OC(C(O)C(O)C1O)C(O)=O)C(O)=O),"
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]-2-[#6]-[#8]-[#6](-[#6](-[#8])-[#6]-2-[#8])-[#6]([#8;A;X2H1,X1-])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6]([#8;A;X2H1,X1-])=O),"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8;R0]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O)-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6]([#8;A;X2H1,X1-])=O),"
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]-2-[#6]-[#8]-[#6](-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6;R0](=[O;R0])-[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6;R0](-[#8;R0])=[O;R0]),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O)"
				+ "]");
		
		SmartsPatternCDK linkedSulfoGlucuronidatedPattern = new SmartsPatternCDK("["
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#6]-1-[#8])-[#6](-[#8])=O),$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](=O)-[#8][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]),"
				+ "$(OC1C(O)C(OC2COC(C(O)C2O)C(O)=O)OC(C1O)C(O)=O),"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8]S([#8;A;X2H1,X1-])(=O)=O)-[#6]-1-[#8])-[#6](-[#8])=O)"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8]S([#8;A;X2H1,X1-])(=O)=O)-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#6]-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]S([#8;A;X2H1,X1-])(=O)=O)-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]-[#6]-[#6]-1-[#8]S([#8;A;X2H1,X1-])(=O)=O)-[#6](-[#8])=O),"
				+ "$(OC1OC(C(O)C(OS(O)(=O)=O)C1O)C(O)=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]S([#8])(=O)=O)-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8]S([#8])(=O)=O)-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6](-[#8]S([#8])(=O)=O)-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#6](-[#8])-[#6](-[#8]S([#8])(=O)=O)-[#8]-[#6](-[#6]-1-[#8])-[#6](-[#8])=O)"				
				
				// Sulfated hexosyl
				+ "$(OC1OC(COS(O)(=O)=O)C(O)C(O)C1O)"
				+ "]");
		
		SmartsPatternCDK glutathione = new SmartsPatternCDK("[#7;R0]-[#6;R0](-[#6;R0]-[#6;R0]-[#6;R0](=[O;R0])-[#7;R0]-[#6;R0](-[#6;R0]-[#16;R0])-[#6;R0](=[O;R0])-[#7;R0]-[#6;R0]-[#6;R0](-[#8;R0])=[O;R0])-[#6;R0](-[#8;R0])=[O;R0]");
		SmartsPatternCDK glutamate = new SmartsPatternCDK("NC(CCC(O)=O)C(O)=O");
		SmartsPatternCDK nAcetylatedProduct = new SmartsPatternCDK("[$([H]C([H])([H])[#6](=O)-[#8][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;a]),"
				+ "$([H][#7;A;X3]([#6;a])[#7;A;X3+0,X4+;!$([N]~[!#6])][#8]-[#6](=O)C([H])([H])[H]),"
				+ "$([#6;a]-[#7;X3R0]-[#8]-[#6](-[#6;H3X4])=O),"
				+ "$([H][#7;A;X3]([#6])[#6](=O)C([H])([H])[H]),"
				+ "$([H][#7;X3](-[#6])-[#7;X3]([H])-[#6;X3](=O)[#6;A;X4]([H])([H])[H]),"
				+ "$([#6;H3X4]-[#6](=O)-[#7;X3;H0,H1][S;X4]([#6])(=[O;X1])=[O;X1]),"
				+ "$([#6;a][#6;A;H2X4][#6;A;H2X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6](-[#6])=O),"
				+ "$([H]C([H])([H])[#6](=O)-[#8;X2]-[#7]-[*;a])"
				+ "]");
		
		
		SmartsPatternCDK oligosaccharides = new SmartsPatternCDK("["
				+ "$([#8]-[#6;R1]-1-[#6;R1]-[#8;R1]-[#6;R1](-[#6;R0]-[#8;R0]-[#6;R1]-2-[#8;R1]-[#6;R1]-[#6;R1](-[#8])-[#6;R1](-[#8])-[#6;R1]-2-[#8])-[#6;R1](-[#8])-[#6;R1]-1-[#8]),"
				+ "$([#8]-[#6]!@-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](!@-[#8;X2]!@-[#6]-2-[#8]-[#6]-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6]-1-[#8]),"
				+ "$([#8]-[#6]!@-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6]-[#6]-1-[#8]-[#8;X2]!@-[#6]-1-[#8]-[#6]-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]),"
				+ "$([#8]-[#6]!@-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#6]-[#6]-1-[#8])-[#8]-[#8;X2]!@-[#6]-1-[#8]-[#6]-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]),"
				+ "$([#8]-[#6]-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6](-[#8])-[#6]-1-[#8]),"
				+ "$([#8]-[#6]-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6]-1-[#8]-[#6](-[#6]-[#8]-[#6]-2-[#8]-[#6](-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])=O)-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]),"
				// Glucuronidated hexosyl
				+ "$(OCC1OC(O)C(O)C(OC2OC(C(O)C(O)C2O)C(O)=O)C1O),"
				+ "$(OC1OC(COC2OC(C(O)C(O)C2O)C(O)=O)C(O)C(O)C1O)"
				+ "]");
		
		SmartsPatternCDK badSterolconjugates = new SmartsPatternCDK("["
				+ "$([#8]-[#6@@H]-1-[#6@@H](-[#8])-[#6@H](-[#8]-[#6]~2-[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~4-,=[#6,#8,#7,#16;A;R2]-,=5-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=5-,=[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]~4~[#6,#8,#7,#16]~2~3)-[#8]-[#6@@H](-[#6@H]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6@@H]-1-[#6@@H](-[#8])-[#6@H](-[#8]-[#6]~2~[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=3-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=3-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~4~[#6,#8,#7,#16]~2~3)-[#8]-[#6@@H](-[#6@H]-1-[#8])-[#6](-[#8])=O),"
				+ "$([#8]-[#6@@H]-1-[#6@@H](-[#8])-[#6@H](-[#8]-[#6,#8,#7,#16]-,=2~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]~3~[#6,#8,#7,#16]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16;A;R2]-,=3-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=2-,=3)-[#8]-[#6@@H](-[#6@H]-1-[#8])-[#6](-[#8])=O),"
					// bad sulated steroid
				+ "$([#8]S(=O)(=O)[#8]-[#6,#8,#7,#16]~1~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~2~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16;A;R2]-,=4-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=4-,=[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]~3~[#6,#8,#7,#16]~1~2),"
				+ "$([#8]S(=O)(=O)[#8]C~1~2~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~1~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~1-,=[#6,#8,#7,#16;A;R2]-,=3-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=3-,=[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]~2~1),"
				+ "$([#8]S(=O)(=O)[#8][#6,#8,#7,#16;A]~1~[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16;A]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16;A;R2]~2~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~3~[#6,#8,#7,#16]~1~2)"
//				+ ""
//				+ ""
				+ "]");
		
		
		
		// 5-(3,4-dihydroxyphenyl)valeric acid must be dehydroxylated at the 4-position of the alkyl chain
		SmartsPatternCDK _34_dhpva = new SmartsPatternCDK("[#8;A;H1X2]!@-[#6]1=,:[#6][#6]=,:[#6]([#6]=,:[#6]1!@-[#8;A;H1X2])[#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#6;A;H2X4][#6;X3]([#8;A;X1-,X2H1])=O");
		
//		System.out.println("IS POLYPHENOL: " + (flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0));

		glycosylMoietyPattern.match(molecule);
		oMethylPattern.match(molecule);
		sulfatedRadicalPattern.match(molecule);
		
//		System.err.println(smartsPattern.hasSMARTSPattern(molecule) > 0);
//		System.err.println("SIZE: " + smartsPattern.match(molecule));
//		System.err.println("SIZE: " + smartsPattern.getUniqueMatchingAtoms().size());

		if (smartsPattern.hasSMARTSPattern(molecule) > 0 || 
				linkedPolyGlucuronidatedPattern.hasSMARTSPattern(molecule) > 0 ||
				linkedSulfoGlucuronidatedPattern.hasSMARTSPattern(molecule) > 0
				) {
			bad = true;
		}	
		else if( (flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0)  
				&& ( (glycosylMoietyPattern.getUniqueMatchingAtoms().size() + oMethylPattern.getUniqueMatchingAtoms().size() + 
				sulfatedRadicalPattern.getUniqueMatchingAtoms().size() ) > 2 || oligosaccharides.hasSMARTSPattern(molecule)>0) ){
			bad = true;
		}
		
		else if (_12_ONS_diglucuronidatedPattern.hasSMARTSPattern(molecule) > 0){
			// O-,N-,S-glucuronidation at two adjacent carbon atoms.
			bad = true;
			
		}
		
		else if(_6H_benzo_c_chromen_6_one.hasSMARTSPattern(molecule) > 0 
				// || _34_dhpva.hasSMARTSPattern(molecule) > 0
				
				){
			bad = true;
		}
		

		else if(anthocyanidinPattern.hasSMARTSPattern(molecule) > 0 && sulfatedRadicalPattern.hasSMARTSPattern(molecule) > 0){
			// No Anthocyanidin sulfate
			bad = true;
		} 
		
		else if( isoflavonePattern.hasSMARTSPattern(molecule) <= 0 && otherFlavonoidsPattern.hasSMARTSPattern(molecule) > 0 && (sulfatedRadicalPattern.hasSMARTSPattern(molecule)>0 && sulfatedRadicalPattern.getUniqueMatchingAtoms().size() > 1)){
			bad = true;
		}


//		else if( (isPolyphenolOrDerivative(molecule) && (isoflavonePattern.hasSMARTSPattern(molecule) == 0)) && sulfatedRadicalPattern.hasSMARTSPattern(molecule) >= 2){

		
		else if (isMetabolizablePolyphenolOrDerivative(molecule)){
//			System.out.println("====> Is metabolizable polyphenol or phenolic derivative");
//			System.out.println("glucuronidePattern.hasSMARTSPattern(molecule): " + glucuronidePattern.hasSMARTSPattern(molecule));
//			System.out.println("sulfatedRadicalPattern.hasSMARTSPattern(molecule): " + sulfatedRadicalPattern.hasSMARTSPattern(molecule));
//			System.out.println("glucuronidePattern.hasSMARTSPattern(molecule) + sulfatedRadicalPattern.hasSMARTSPattern(molecule): " + glucuronidePattern.hasSMARTSPattern(molecule) + sulfatedRadicalPattern.hasSMARTSPattern(molecule));
//			System.out.println(sulfatedRadicalPattern.getUniqueMatchingAtoms(molecule));
//			System.out.println(sulfatedRadicalPattern.getUniqueMatchingAtoms().size());
//			System.out.println(sulfatedRadicalPattern.getUniqueMatchingAtoms());
			// 2 or more linked glucuronide molecules
			
			
			
			
			if(linkedPolyGlucuronidatedPattern.hasSMARTSPattern(molecule)>0){
				bad = true;
			}
			// no more than 2 glucuronide moieties
			// No disulfates, except for isoflavones
					
			
			else if ( glucuronidePattern.hasSMARTSPattern(molecule)>0 && 			
					(glucuronidePattern.getUniqueMatchingAtoms().size() > 2) || (isoflavonePattern.hasSMARTSPattern(molecule) == 0 && glucuronidePattern.getUniqueMatchingAtoms().size() >= 2 )){
				bad = true;
			}
			else if ( glucuronidePattern.hasSMARTSPattern(molecule)>0 &&  sulfatedRadicalPattern.hasSMARTSPattern(molecule)>0 &&
					( glucuronidePattern.getUniqueMatchingAtoms().size() + sulfatedRadicalPattern.getUniqueMatchingAtoms().size() > 2)){
				// System.out.println("glucuronidePattern + sulfatedRadicalPattern");
				bad = true;
			}
////			else if (_2p_glucuronidatedCatechol.hasSMARTSPattern(molecule) > 0 || _2p_sulfatedCatechol.hasSMARTSPattern(molecule) > 0 || adjacentSulfateGlucuronideSubstitutionsOnCatechol.hasSMARTSPattern(molecule)>0){
			else if (adjacentSulfateGlucuronideSubstitutionsOnCatechol.hasSMARTSPattern(molecule)>0){
//				System.out.println("adjacentSulfateGlucuronideSubstitutionsOnCatechol");
				bad = true;
			}			
//			// if it contains glutamate, n-acetyl, n-glutathione transferease()
//			else if ( glutathione.hasSMARTSPattern(molecule) > 0 || glutamate.hasSMARTSPattern(molecule) > 0 || nAcetylatedProduct.hasSMARTSPattern(molecule) > 0){
//				bad =  true;
//			}

		}
		
				
		// ADD OTHER ISOFLAVONES DICONJUGATES (MUST BE COMBINATIONS OF SULFATE AND GLUCURONIDE)
		// See FLAVONOIDS Chemistry, Biochemistry and Applications. P. 328 (6.4.5)
		
//		if(bad == true){
//			System.err.println("This metabolite is invalid and will be removed.");
//		}
		
		return bad;
	}

	public static boolean isPhaseIPolyphenolCandidateOrDerivative(IAtomContainer molecule) throws SMARTSException, CDKException, CloneNotSupportedException{
		boolean ppc = false;
		
		if(!(ChemicalClassFinder.isGlycosylatedCompound(molecule) || ChemicalClassFinder.isSulfatedCompound(molecule) || ChemicalClassFinder.isGlycinatedCompound(molecule))
				&& ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(molecule)) {
			ppc = true;
		}
		
		return ppc;
		
	}
	
	public static boolean isMixture(IAtomContainer molecule) throws CDKException{
		// compound is not a mixture (checkConnectivity returns 2 or more atomContainers)
		boolean mixture = ConnectivityChecker.partitionIntoMolecules(molecule).getAtomContainerCount()>1;
		return mixture;	
	}

	public static boolean isPpsValid(IAtomContainer molecule) throws CDKException{
		// http://eawag-bbd.ethz.ch/predict/notbepredicted.html
		
		String inchikey = molecule.getProperty("InChIKey");
		if(inchikey ==null){
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();			
		}
		
		boolean valid = (!(isPpsCofactor(inchikey) || isPpsDeadEndCompound(inchikey)))  &&  AtomContainerManipulator.getNaturalExactMass(molecule)<1000.0 && 
				containsCarbon(molecule) && !isMixture(molecule) && 
				( (numberOfAtomWithAtomicNumber(molecule,1) + numberOfAtomWithAtomicNumber(molecule, 6) + 
						numberOfAtomWithAtomicNumber(molecule, 7) + numberOfAtomWithAtomicNumber(molecule, 8) + 
						numberOfAtomWithAtomicNumber(molecule, 15) + numberOfAtomWithAtomicNumber(molecule, 16)
						+ numberOfAtomWithAtomicNumber(molecule, 9) + numberOfAtomWithAtomicNumber(molecule, 17) 
						+ numberOfAtomWithAtomicNumber(molecule, 35) + numberOfAtomWithAtomicNumber(molecule, 53) 
						== molecule.getAtomCount() ) );
				
				//  && !isCompoundHighlyfluorinated(molecule)

		return valid;
	}
	
	public static boolean isBtValid(IAtomContainer molecule) throws CDKException{
		// http://eawag-bbd.ethz.ch/predict/notbepredicted.html
		
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen1 = factory.getInChIGenerator(molecule);
		String inchikey = gen1.getInchiKey();
		
		boolean valid = AtomContainerManipulator.getNaturalExactMass(molecule)<1000.0 && containsCarbon(molecule) && !isMixture(molecule) ;

		return valid;
	}	


	
	
	
//	public static boolean isPerfluoroDerivative(IAtomContainer molecule){
//		boolean pfd = false;		
//		// Perfluoroacetone transformation
////		Pattern smp = SmartsPattern.create("");
////		if()								
//		return pfd;
//	}

	/**
	 * 
	 * @param molecule
	 * @param atomicNumber
	 * @return the number of atoms with the given atomic number 
	 */
	public static int numberOfAtomWithAtomicNumber(IAtomContainer molecule,int atomicNumber){
		int atCount = 0;
		
		for(int l = 0; l <molecule.getAtomCount(); l++){
			if(molecule.getAtom(l).getAtomicNumber() == atomicNumber)
				atCount++;
		}
		return atCount;
	}

	public static boolean isPpsCofactor(String inchikey){
		return (ppsCofactors.containsKey(inchikey));	
	}
	
	public static boolean isPpsCofactor(IAtomContainer molecule) throws CDKException{
		String inchikey = molecule.getProperty("InChIKey");
		if( inchikey != null){
			return (ppsCofactors.containsKey(inchikey));
		} else{
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
			return (ppsCofactors.containsKey(inchikey));
		}
			
	}
	
	public static boolean isPpsDeadEndCompound(String inchikey){
		return (ppsDeadEndCompounds.containsKey(inchikey));	
	}
	
	
	public static boolean isPpsDeadEndCompound(IAtomContainer molecule) throws CDKException{
		String inchikey = molecule.getProperty("InChIKey");
		if( inchikey != null){
			return (ppsDeadEndCompounds.containsKey(inchikey));	
		} else{
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
			return (ppsDeadEndCompounds.containsKey(inchikey));	
		}
		
		
	}
	
	public static boolean isBioTransformerValidStrict(IAtomContainer molecule) throws CDKException{
		
		String inchikey = molecule.getProperty("InChIKey");
		if(inchikey ==null){
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
//			System.out.println("IS BIOTRANSFORMER VALID INCHIKEY: " + inchikey);
		}
 
//		System.out.println(!( isPpsCofactor(inchikey) || isPpsDeadEndCompound(inchikey) || 
//				isCompoundInorganic(molecule) || isStandardAminoAcid(molecule)));
		return !(isMixture(molecule) || isPpsCofactor(inchikey) || isPpsDeadEndCompound(inchikey) || 
				isCompoundInorganic(molecule) || isStandardAminoAcid(molecule)); 
	}
	

	
	public static boolean isBioTransformerValid(IAtomContainer molecule) throws CDKException{
		
		String inchikey = molecule.getProperty("InChIKey");
		if(inchikey ==null){
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
//			System.out.println("IS BIOTRANSFORMER VALID INCHIKEY: " + inchikey);
		}
//		boolean a = isPpsCofactor(inchikey);
//		boolean b = isPpsDeadEndCompound(inchikey);
//		boolean c = isStandardAminoAcid(molecule);
		return !(isPpsCofactor(inchikey) || isPpsDeadEndCompound(inchikey) || 
				isStandardAminoAcid(molecule)); 
	}
	
	public static void addInChIandKey(IAtomContainer molecule) throws CDKException{
		if(molecule.getProperty("InChIKey") == null){
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			String inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());			
		}

	}
	
	public static double getMajorIsotopeMass(IAtomContainer molecule){	
        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);		
		return MolecularFormulaManipulator.getMajorIsotopeMass(formula);
	}

	public static String getMolecularFormula(IAtomContainer molecule){	
        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);		
		return MolecularFormulaManipulator.getString(formula);
	}
	
	public static LinkedHashMap<String, String> computePhysicoChemicalProperties(IAtomContainer molecule) throws CDKException {
		LinkedHashMap<String, String> properties = new LinkedHashMap<String, String> ();
		
		IDescriptorResult xlogp 		= null;
		IDescriptorResult alogp 		= null;
		IDescriptorResult weight 	= null;
		Double majorIsotopeMass 	= null;
		
		
		xlogp = xLogpDescriptor.calculate(molecule).getValue();
		ALOGPDescriptor aLogpDescriptor = new ALOGPDescriptor();
		alogp = aLogpDescriptor.calculate(molecule).getValue();

		
		// Calculate the weight of specified element type in the supplied
//		WeightDescriptor weightD = new WeightDescriptor();
//		weight = weightD.calculate(molecule).getValue();
		
		// Get the summed major isotopic mass of all elements from an MolecularFormula.
		majorIsotopeMass = getMajorIsotopeMass(molecule);
		
//		DescriptorValue pka = PKASmartsDescriptor
		
		IDescriptorResult hbaCount 	= hbaDCountDescriptor.calculate(molecule).getValue();
		IDescriptorResult hbdCount 	= hbdCountDescriptor.calculate(molecule).getValue();
		IDescriptorResult rbCount	= rbCountDescriptor.calculate(molecule).getValue();
		
		

		properties.put("Major Isotope Mass" , majorIsotopeMass.toString());
		properties.put("ALogP", alogp.toString().split(",")[0]);
		properties.put("XLogP", xlogp.toString());		
//		properties.put("hbaCount", hbaCount.toString());
//		properties.put("hbdCount", hbdCount.toString());
//		properties.put("rbCount", rbCount.toString());
//		System.out.println("properties: " + properties);
//		properties.put("Molecular weight" , weight.toString());
		return properties;
		
		
	}
	
	public static boolean isCompoundInorganic(IAtomContainer molecule){
		//		boolean kingdom = false;		
		return !(containsCarbon(molecule));
	}
	
	public static boolean isUnneccessaryMetabolite(IAtomContainer molecule) throws SMARTSException{

		String patterns = "["
				+ "$([H][H]),"	// Dihydrogen
				+ "$(O=C=O)," 	// Carbon dioxide
				+ "$([#6;A;H2X3]=[O;X1]),"   	// Carbon monoxide
				+ "$([OX2H2])," 	// Water
				+ "$([H])," // Hydrogen
				+ "$([OX1]=[SX1]),"
				+ "$([NX3H3]),"
				+ "$([CX4H3])" // Methane
				+ "$([F,Cl,Br,I;-]),"
				+ "$([F,Cl,Br,I;H1X1]),"
				+ "$([#8;X2H1,X1-][S;X4]([#8;X2H1,X1-])(=[O;X1])=[O;X1])," // Sulfate
				+ "$([#8;X2H1,X1-][P;X4]([#8;X2H1,X1-])([#8;X2H1,X1-])=[O;X1])" // Phosphate
				+ "]";
		
		SmartsPatternCDK cdkPattern =  new SmartsPatternCDK(patterns);
		return (cdkPattern.hasSMARTSPattern(molecule)>0);
	}

	public static boolean isStandardAminoAcid(IAtomContainer molecule) throws CDKException{
		boolean aa = false;
		
		String inchikey = molecule.getProperty("InChIKey");
		if(inchikey ==null){
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();			
		}
		
		String strippedInChiKey = inchikey.split("-")[0];

		return standardAminoAcids.containsKey(strippedInChiKey);
	}
	
	
	private static LinkedHashMap<String, String[] > ppsCofactors;
	
	static{
		ppsCofactors = new LinkedHashMap<String, String[]>();
		ppsCofactors.put("ZSLZBFCDCINBPY-NGQYVWNKSA-N", new String[]{"acetyl-CoASH","CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)N1C=NC2=C(N)N=CN=C12"});
		ppsCofactors.put("QGZKDVFQNNGYKY-UHFFFAOYSA-N", new String[]{"ammonia","N"});
		ppsCofactors.put("BVKZGUZCCUSVTD-UHFFFAOYSA-M", new String[]{"bicarbonate","OC([O-])=O"});
		ppsCofactors.put("CURLTUGMZLYLDI-UHFFFAOYSA-N", new String[]{"carbon dioxide","O=C=O"});
		ppsCofactors.put("UGFAIRIUMAVXCW-UHFFFAOYSA-N", new String[]{"carbon monoxide","[C-]#[O+]"});
		ppsCofactors.put("BVKZGUZCCUSVTD-UHFFFAOYSA-N", new String[]{"carbonate","OC(O)=O"});
		ppsCofactors.put("RGJOEKWQDUBAIZ-UHFFFAOYSA-N", new String[]{"CoenzymeASH","CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS"});
		ppsCofactors.put("OTMSDBZUPAUEDD-UHFFFAOYSA-N", new String[]{"ethane","CC"});
		ppsCofactors.put("RWSXRVCMGQZWBV-UHFFFAOYSA-N", new String[]{"glutathione","NC(CCC(=O)NC(CS)C(=O)NCC(O)=O)C(O)=O"});
		ppsCofactors.put("RWSOTUBLDIXVET-UHFFFAOYSA-N", new String[]{"hydrogen sulfide","S"});
		ppsCofactors.put("IOVCWXUNBOPUCH-UHFFFAOYSA-M", new String[]{"nitrite","[O-]N=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-N", new String[]{"phosphate","OP(O)(O)=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-L", new String[]{"phosphate dianion","OP([O-])([O-])=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-K", new String[]{"phosphate trianion","[O-]P([O-])([O-])=O"});
		ppsCofactors.put("QAOWNCQODCNURD-UHFFFAOYSA-L", new String[]{"sulfate","[O-]S([O-])(=O)=O"});
		ppsCofactors.put("LSNNMFCWUKXFEE-UHFFFAOYSA-L", new String[]{"sulfite","[O-]S([O-])=O"});
		ppsCofactors.put("LSNNMFCWUKXFEE-UHFFFAOYSA-N", new String[]{"sulfurous acid","OS(O)=O"});
		ppsCofactors.put("XPRXJFAFJFZWTC-UHFFFAOYSA-M", new String[]{"trioxidosulfate","[O]S([O-])=O"});
	
	}
	
	
	/**
	 * A dictionary with dead-end compounds. Remember that AMBIT does not return stereo-specific configurations
	 */
	private static LinkedHashMap<String, String[] > ppsDeadEndCompounds;

	static{
		ppsDeadEndCompounds = new LinkedHashMap<String, String[]>();
		ppsDeadEndCompounds.put("GTZCVFVGUGFEME-IWQZZHSRSA-M", new String[]{"cis-aconitate","OC(=O)\\C=C(\\CC([O-])=O)C(O)=O"});
		ppsDeadEndCompounds.put("IKHGUXGNUITLKF-UHFFFAOYSA-N", new String[]{"acetaldehyde","CC=O"});
		ppsDeadEndCompounds.put("QTBSBXVTEAMEQO-UHFFFAOYSA-N", new String[]{"acetate","CC(O)=O"});
		ppsDeadEndCompounds.put("WDJHALXBUFZDSR-UHFFFAOYSA-N", new String[]{"acetoacetate","CC(=O)CC(O)=O"});
		ppsDeadEndCompounds.put("CSCPPACGZOOCGX-UHFFFAOYSA-N", new String[]{"acetone","CC(C)=O"});
		ppsDeadEndCompounds.put("HSFWRNGVRCDJHI-UHFFFAOYSA-N", new String[]{"acetylene","C#C"});
		ppsDeadEndCompounds.put("GFFGJBXGBJISGV-UHFFFAOYSA-N", new String[]{"adenine","NC1=C2N=CNC2=NC=N1"});
		ppsDeadEndCompounds.put("WNLRTRBMVRJNCN-UHFFFAOYSA-N", new String[]{"adipate","OC(=O)CCCCC(O)=O"});
		ppsDeadEndCompounds.put("QNAYBMKLOCPYGJ-UHFFFAOYSA-N", new String[]{"alanine","CC(N)C(O)=O"});
		ppsDeadEndCompounds.put("WQZGKKKJIJFFOK-DVKNGEFBSA-N", new String[]{"alpha-D-glucose","OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"});
		ppsDeadEndCompounds.put("VUTBELPREDJDDH-UHFFFAOYSA-N", new String[]{"4-amino-5-hydroxymethyl-2-methylpyrimidine","CC1=NC(N)=C(CO)C=N1"});
		ppsDeadEndCompounds.put("TYQCGQRIZGCHNB-JLAZNSOCSA-N", new String[]{"ascorbate","OC[C@H](O)[C@H]1OC(O)=C(O)C1=O"});
		ppsDeadEndCompounds.put("OHJMTUPIZMNBFR-UHFFFAOYSA-N", new String[]{"biuret","NC(=O)NC(N)=O"});
		ppsDeadEndCompounds.put("LRHPLDYGYMQRHN-UHFFFAOYSA-N", new String[]{"1-butanol","CCCCO"});
		ppsDeadEndCompounds.put("FERIUCNNQQJTOY-UHFFFAOYSA-N", new String[]{"butyrate","CCCC(O)=O"});
		ppsDeadEndCompounds.put("KXDHJXZQYSOELW-UHFFFAOYSA-N", new String[]{"carbamate","NC(O)=O"});
		ppsDeadEndCompounds.put("XLJMAIOERFSOGZ-UHFFFAOYSA-N", new String[]{"cyanate","OC#N"});
		ppsDeadEndCompounds.put("XFXPMWWXUTWYJX-UHFFFAOYSA-N", new String[]{"cyanide","[C-]#N"});
		ppsDeadEndCompounds.put("UFULAYFCSOUIOV-UHFFFAOYSA-N", new String[]{"cysteamine","NCCS"});
		ppsDeadEndCompounds.put("XUJNEKJLAYXESH-UHFFFAOYSA-N", new String[]{"cysteine","NC(CS)C(O)=O"});
		ppsDeadEndCompounds.put("OPTASPLRGRRNAP-UHFFFAOYSA-N", new String[]{"cytosine","NC1=NC(=O)NC=C1"});
		ppsDeadEndCompounds.put("WQZGKKKJIJFFOK-GASJEMHNSA-N", new String[]{"D-glucose","OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"});
		ppsDeadEndCompounds.put("GHVNFZFCNZKVNT-UHFFFAOYSA-N", new String[]{"decanoate","CCCCCCCCCC(O)=O"});
		ppsDeadEndCompounds.put("QSJXEFYPDANLFS-UHFFFAOYSA-N", new String[]{"diacetyl","CC(=O)C(C)=O"});
		ppsDeadEndCompounds.put("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", new String[]{"ethanol","CCO"});
		ppsDeadEndCompounds.put("HZAXFHJVJLSVMW-UHFFFAOYSA-N", new String[]{"ethanolamine","NCCO"});
		ppsDeadEndCompounds.put("LYCAIKOWRPUZTN-UHFFFAOYSA-N", new String[]{"ethylene glycol","OCCO"});
		ppsDeadEndCompounds.put("WSFSSNUMVMOOMR-UHFFFAOYSA-N", new String[]{"formaldehyde","C=O"});
		ppsDeadEndCompounds.put("ZHNUHDYFZUAESO-UHFFFAOYSA-N", new String[]{"formamide","NC=O"});
		ppsDeadEndCompounds.put("BDAGIHXWWSANSR-UHFFFAOYSA-N", new String[]{"formate","OC=O"});
		ppsDeadEndCompounds.put("VZCYOOQTPOCHFL-OWOJBTEDSA-N", new String[]{"fumarate","OC(=O)\\C=C\\C(O)=O"});
		ppsDeadEndCompounds.put("JFCQEDHGNNZCLN-UHFFFAOYSA-N", new String[]{"glutarate","OC(=O)CCCC(O)=O"});
		ppsDeadEndCompounds.put("PEDCQBHIVMGVHV-UHFFFAOYSA-N", new String[]{"glycerol","OCC(O)CO"});
		ppsDeadEndCompounds.put("DHMQDGOQFOQNFH-UHFFFAOYSA-N", new String[]{"glycine","NCC(O)=O"});
		ppsDeadEndCompounds.put("WGCNASOHLSPBMP-UHFFFAOYSA-N", new String[]{"glycoladehyde","O=CCO"});
		ppsDeadEndCompounds.put("AEMRFAOFKBGASW-UHFFFAOYSA-N", new String[]{"glycolate","OCC(O)=O"});
		ppsDeadEndCompounds.put("HHLFWLYXYJOTON-UHFFFAOYSA-N", new String[]{"glyoxylate","OC(=O)C=O"});
		ppsDeadEndCompounds.put("UYTPUPDQBNUYGX-UHFFFAOYSA-N", new String[]{"guanine","NC1=NC2=C(N=CN2)C(=O)N1"});
		ppsDeadEndCompounds.put("AFENDNXGAFYKQO-UHFFFAOYSA-N", new String[]{"2-hydroxybutyrate","CCC(O)C(O)=O"});
		ppsDeadEndCompounds.put("WHBMMWSBFZVSSR-UHFFFAOYSA-N", new String[]{"3-hydroxybutyrate","CC(O)CC(O)=O"});
		ppsDeadEndCompounds.put("SJZRECIVHVDYJC-UHFFFAOYSA-N", new String[]{"4-hydroxybutyrate","OCCCC(O)=O"});
		ppsDeadEndCompounds.put("ALRHLSYJTWAHJZ-UHFFFAOYSA-N", new String[]{"3-hydroxypropanoate","OCCC(O)=O"});
		ppsDeadEndCompounds.put("HHDDCCUIIUWNGJ-UHFFFAOYSA-N", new String[]{"hydroxypyruvate","O=C(O)C(=O)CO"});
		ppsDeadEndCompounds.put("KIAHVTZFUFPMSB-UHFFFAOYSA-N", new String[]{"imidazole pyruvate","CC(=O)C(=O)ON1C=CN=C1"});
		ppsDeadEndCompounds.put("KQNPFQTWMSNSAP-UHFFFAOYSA-N", new String[]{"isobutyrate","CC(C)C(O)=O"});
		ppsDeadEndCompounds.put("GWYFCOCPABKNJV-UHFFFAOYSA-N", new String[]{"isovalerate","CC(C)CC(O)=O"});
		ppsDeadEndCompounds.put("TYEYBOSBBBHJIV-UHFFFAOYSA-N", new String[]{"2-ketobutyrate","CCC(=O)C(=O)O"});
		ppsDeadEndCompounds.put("TYEYBOSBBBHJIV-UHFFFAOYSA-N", new String[]{"L-2-oxobutyrate","CCC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("JVTAAEKCZFNVCJ-UHFFFAOYSA-N", new String[]{"lactate","CC(O)C(O)=O"});
		ppsDeadEndCompounds.put("POULHZVOKOAJMA-UHFFFAOYSA-N", new String[]{"lauric acid","CCCCCCCCCCCC(O)=O"});
		ppsDeadEndCompounds.put("BJEPYKJPYRNKOW-UHFFFAOYSA-N", new String[]{"malate","OC(CC(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("FSQQTNAZHBEJLS-UPHRSURJSA-N", new String[]{"maleamate","NC(=O)\\C=C/C(O)=O"});
		ppsDeadEndCompounds.put("VZCYOOQTPOCHFL-UPHRSURJSA-N", new String[]{"maleate","OC(=O)\\C=C/C(O)=O"});
		ppsDeadEndCompounds.put("OFOBLEOULBTSOW-UHFFFAOYSA-N", new String[]{"malonate","C(C(=O)O)C(=O)O"});
		ppsDeadEndCompounds.put("OFOBLEOULBTSOW-UHFFFAOYSA-N", new String[]{"malonate semialdehyde","OC(=O)CC(O)=O"});
		ppsDeadEndCompounds.put("VNWKTOKETHGBQD-UHFFFAOYSA-N", new String[]{"methane","C"});
		ppsDeadEndCompounds.put("OKKJLVBELUTLKV-UHFFFAOYSA-N", new String[]{"methanol","CO"});
		ppsDeadEndCompounds.put("BAVYZALUXZFZLV-UHFFFAOYSA-N", new String[]{"methylamine","CN"});
		ppsDeadEndCompounds.put("YYPNJNDODFVZLE-UHFFFAOYSA-N", new String[]{"3-methylcrotonate","CC(C)=CC(O)=O"});
		ppsDeadEndCompounds.put("HNEGQIOMVPPMNR-IHWYPQMZSA-N", new String[]{"2-methylmaleate","C\\C(=C\\C(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("ZIYVHBGGAOATLY-UHFFFAOYSA-N", new String[]{"methylmalonate","CC(C(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("PVNIIMVLHYAWGP-UHFFFAOYSA-N", new String[]{"niacin","OC(=O)C1=CC=CN=C1"});
		ppsDeadEndCompounds.put("WWZKQHOCKIZLMA-UHFFFAOYSA-N", new String[]{"ocatanoate","CCCCCCCC(=O)O"});
		ppsDeadEndCompounds.put("MUBZPKHOEPUJKR-UHFFFAOYSA-N", new String[]{"oxalate","OC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("KHPXUQMNIQBQEV-UHFFFAOYSA-N", new String[]{"oxaloacetate","OC(=O)CC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("KPGXRSRHYNQIFN-UHFFFAOYSA-N", new String[]{"2-oxoglutarate","OC(=O)CCC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("NOXRYJAWRSNUJD-UHFFFAOYSA-N", new String[]{"2-oxopent-4-enoate","OC(=O)C(=O)CC=C"});
		ppsDeadEndCompounds.put("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", new String[]{"propane","CCC"});
		ppsDeadEndCompounds.put("KFZMGEQAYNKOFK-UHFFFAOYSA-N", new String[]{"2-propanol","CC(C)O"});
		ppsDeadEndCompounds.put("XBDQKXXYIPTUBI-UHFFFAOYSA-N", new String[]{"propionate","CCC(O)=O"});
		ppsDeadEndCompounds.put("DNIAPMSPPWPWGF-UHFFFAOYSA-N", new String[]{"propylene glycol","CC(O)CO"});
		ppsDeadEndCompounds.put("LCTONWCANYUPML-UHFFFAOYSA-N", new String[]{"pyruvate","CC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("FSYKKLYZXJSNPZ-UHFFFAOYSA-N", new String[]{"sarcosine","CNCC(O)=O"});
		ppsDeadEndCompounds.put("KDYFGRWQOYBRFD-UHFFFAOYSA-N", new String[]{"succinate","OC(=O)CCC(O)=O"});
		ppsDeadEndCompounds.put("AYFVYJQAPQTCCC-GBXIJSLDSA-N", new String[]{"threonine","C[C@@H](O)[C@H](N)C(O)=O"});
		ppsDeadEndCompounds.put("RWQNBRDOKXIBIV-UHFFFAOYSA-N", new String[]{"thymine","CC1=CNC(=O)NC1=O"});
		ppsDeadEndCompounds.put("ISAKRJDGNUQOIC-UHFFFAOYSA-N", new String[]{"uracil","O=C1NC=CC(=O)N1"});
		ppsDeadEndCompounds.put("XSQUKJJJFZCRTK-UHFFFAOYSA-N", new String[]{"urea","NC(N)=O"});
		ppsDeadEndCompounds.put("LRFVTYWOQMYALW-UHFFFAOYSA-N", new String[]{"xanthine","O=C1NC2=C(NC=N2)C(=O)N1"});
//		ppsDeadEndCompounds.put("", new String[]{"",""});
	}
	
		
	/**
	 * A dictionary of amino acids. considered only the first part of the inchikeys
	 */	
	private static LinkedHashMap<String, String[] > standardAminoAcids;

	static{
		standardAminoAcids = new LinkedHashMap<String, String[]>();
		standardAminoAcids.put("QNAYBMKLOCPYGJ", new String[]{"alanine","CC(N)C(O)=O"});
		standardAminoAcids.put("ODKSFYDXXFIFQN", new String[]{"arginine","NC(CCCNC(N)=N)C(O)=O"});
		standardAminoAcids.put("DCXYFEDJOCDNAF", new String[]{"asparagine","NC(CC(N)=O)C(O)=O"});
		standardAminoAcids.put("CKLJMWTZIZZHCS", new String[]{"aspartic acid","NC(CC(O)=O)C(O)=O"});
		standardAminoAcids.put("XUJNEKJLAYXESH", new String[]{"cysteine","NC(CS)C(O)=O"});
		standardAminoAcids.put("WHUUTDBJXJRKMK", new String[]{"glutamic acid","NC(CCC(O)=O)C(O)=O"});
		standardAminoAcids.put("ZDXPYRJPNDTMRX", new String[]{"glutamine","NC(CCC(N)=O)C(O)=O"});
		standardAminoAcids.put("DHMQDGOQFOQNFH", new String[]{"glycine","NCC(O)=O"});
		standardAminoAcids.put("HNDVDQJCIGZPNO", new String[]{"histidine","NC(CC1=CN=CN1)C(O)=O"});
		standardAminoAcids.put("AGPKZVBTJJNPAG", new String[]{"isoleucine","CC[C@H](C)[C@H](N)C(O)=O"});
		standardAminoAcids.put("ROHFNLRQFUQHCH", new String[]{"leucine","CC(C)CC(N)C(O)=O"});
		standardAminoAcids.put("KDXKERNSBIXSRK", new String[]{"lysine","NCCCCC(N)C(O)=O"});
		standardAminoAcids.put("FFEARJCKVFRZRR", new String[]{"methionine","CSCCC(N)C(O)=O"});
		standardAminoAcids.put("COLNVLDHVKWLRT", new String[]{"phenylalanine","N[C@@H](CC1=CC=CC=C1)C(O)=O"});
		standardAminoAcids.put("ONIBWKKTOPOVIA", new String[]{"proline","OC(=O)C1CCCN1"});
		standardAminoAcids.put("MTCFGRXMJLQNBG", new String[]{"serine","N[C@@H](CO)C(O)=O"});
		standardAminoAcids.put("AYFVYJQAPQTCCC", new String[]{"threonine","C[C@@H](O)[C@H](N)C(O)=O"});
		standardAminoAcids.put("QIVBCDIJIAJPQS", new String[]{"tryptophan",""});
		standardAminoAcids.put("QIVBCDIJIAJPQS", new String[]{"tyrosine","N[C@@H](CC1=CNC2=C1C=CC=C2)C(O)=O"});
		standardAminoAcids.put("KZSNJWFQEVHDMF", new String[]{"valine","CC(C)[C@H](N)C(O)=O"});
		standardAminoAcids.put("ZKZBPNGNEQAJSX", new String[]{"selenocysteine","N[C@@H](C[SeH])C(O)=O"});
		standardAminoAcids.put("ZFOMKMMPBOQKMC", new String[]{"pyrrolysine","C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N)C(O)=O"});
	}

	
	/**
	 * Cofactors of mammalian metabolism
	 * https://en.wikipedia.org/wiki/Cofactor_(biochemistry)
	 */		
	
	// 
	private static LinkedHashMap<String, String[] > mammalianCofactors;
	static{
		mammalianCofactors = new LinkedHashMap<String, String[]>();
		mammalianCofactors.put("BAWFJGJZGIEFAR", new String[]{"NAD+","O=C(N)c1ccc[n+](c1)[C@@H]2O[C@@H]([C@@H](O)[C@H]2O)COP([O-])(=O)OP(=O)(O)OC[C@H]5O[C@@H](n4cnc3c(ncnc34)N)[C@H](O)[C@@H]5O"});
		mammalianCofactors.put("XJLXINKUBYWONI", new String[]{"NADPH","NC(=O)c1ccc[n+](c1)C1OC(COP([O-])(=O)OP(O)(=O)OCC2OC(C(OP(O)(O)=O)C2O)n2cnc3c(N)ncnc23)C(O)C1O"});
		mammalianCofactors.put("YTNIXZGTHTVJBW", new String[]{"Flavin mononucleotide","Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)N(CC(O)C(O)C(O)COP(O)(O)=O)c2cc1C"});
		mammalianCofactors.put("VWWQXMAJTJZDQX", new String[]{"Flavin adenine dinucleotide","Cc1cc2N=C3C(=O)NC(=O)N=C3N(CC(O)C(O)C(O)COP(O)(=O)OP(O)(=O)OCC3OC(C(O)C3O)n3cnc4c(N)ncnc34)c2cc1C"});
		mammalianCofactors.put("XXFFZHQKSZLLIT", new String[]{"Coenzyme F420","CC(C(=O)NC(CCC(=O)NC(CCC(=O)O)C(=O)O)C(=O)O)OP(=O)(O)OCC(C(C(CN1C2=CC(=O)C=CC2=CC3=C1NC(=O)NC3=O)O)O)O"});
//		mammalianCofactors.put("", new String[]{"Adenosine 5'-triphosphate",""});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
	
	}
	private static LinkedHashMap<String, String[] > PolyhenolDeadEndAglyconeMetabolites;
	static{
		mammalianCofactors = new LinkedHashMap<String, String[]>();
		mammalianCofactors.put("YCIMNLLNPGFGHC", new String[]{"catechol","OC1=CC=CC=C1O"});
		mammalianCofactors.put("WQGWDDDVZFFDIG", new String[]{"benzene‐1,2,3‐triol","OC1=CC=CC(O)=C1O"});
		mammalianCofactors.put("WPYMKLBDIGXBTP", new String[]{"benzoic acid","OC(=O)C1=CC=CC=C1"});
		mammalianCofactors.put("QCDYQQDYXPDABM", new String[]{"phloroglucinol","OC1=CC(O)=CC(O)=C1"});
		mammalianCofactors.put("", new String[]{"",""});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
//		mammalianCofactors.put("", new String[]{});
	
	}	
//	static{
//		standardAminoAcids.put("QNAYBMKLOCPYGJ-UHFFFAOYNA-N", new String[]{"alanine","CC(N)C(O)=O"});
//		standardAminoAcids.put("ODKSFYDXXFIFQN-UHFFFAOYNA-N", new String[]{"arginine","NC(CCCNC(N)=N)C(O)=O"});
//		standardAminoAcids.put("DCXYFEDJOCDNAF-UHFFFAOYSA-N", new String[]{"asparagine","NC(CC(N)=O)C(O)=O"});
//		standardAminoAcids.put("CKLJMWTZIZZHCS-UHFFFAOYNA-N", new String[]{"aspartic acid","NC(CC(O)=O)C(O)=O"});
//		standardAminoAcids.put("XUJNEKJLAYXESH-UHFFFAOYNA-N", new String[]{"cysteine","NC(CS)C(O)=O"});
//		standardAminoAcids.put("WHUUTDBJXJRKMK-UHFFFAOYSA-N", new String[]{"glutamic acid","NC(CCC(O)=O)C(O)=O"});
//		standardAminoAcids.put("ZDXPYRJPNDTMRX-UHFFFAOYSA-N", new String[]{"glutamine","NC(CCC(N)=O)C(O)=O"});
//		standardAminoAcids.put("DHMQDGOQFOQNFH-UHFFFAOYSA-N", new String[]{"glycine","NCC(O)=O"});
//		standardAminoAcids.put("HNDVDQJCIGZPNO-UHFFFAOYSA-N", new String[]{"histidine","NC(CC1=CN=CN1)C(O)=O"});
//		standardAminoAcids.put("AGPKZVBTJJNPAG-WHFBIAKZSA-N", new String[]{"isoleucine","CC[C@H](C)[C@H](N)C(O)=O"});
//		standardAminoAcids.put("ROHFNLRQFUQHCH-UHFFFAOYNA-N", new String[]{"leucine","CC(C)CC(N)C(O)=O"});
//		standardAminoAcids.put("KDXKERNSBIXSRK-UHFFFAOYNA-N", new String[]{"lysine","NCCCCC(N)C(O)=O"});
//		standardAminoAcids.put("FFEARJCKVFRZRR-UHFFFAOYNA-N", new String[]{"methionine","CSCCC(N)C(O)=O"});
//		standardAminoAcids.put("COLNVLDHVKWLRT-QMMMGPOBSA-N", new String[]{"phenylalanine","N[C@@H](CC1=CC=CC=C1)C(O)=O"});
//		standardAminoAcids.put("ONIBWKKTOPOVIA-UHFFFAOYNA-N", new String[]{"proline","OC(=O)C1CCCN1"});
//		standardAminoAcids.put("MTCFGRXMJLQNBG-REOHCLBHSA-N", new String[]{"serine","N[C@@H](CO)C(O)=O"});
//		standardAminoAcids.put("AYFVYJQAPQTCCC-GBXIJSLDSA-N", new String[]{"threonine","C[C@@H](O)[C@H](N)C(O)=O"});
//		standardAminoAcids.put("QIVBCDIJIAJPQS-VIFPVBQESA-N", new String[]{"tryptophan",""});
//		standardAminoAcids.put("QIVBCDIJIAJPQS-VIFPVBQESA-N", new String[]{"tyrosine","N[C@@H](CC1=CNC2=C1C=CC=C2)C(O)=O"});
//		standardAminoAcids.put("KZSNJWFQEVHDMF-BYPYZUCNSA-N", new String[]{"valine","CC(C)[C@H](N)C(O)=O"});
//		standardAminoAcids.put("ZKZBPNGNEQAJSX-REOHCLBHSA-N", new String[]{"selenocysteine","N[C@@H](C[SeH])C(O)=O"});
//		standardAminoAcids.put("ZFOMKMMPBOQKMC-KXUCPTDWSA-N", new String[]{"pyrrolysine","C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N)C(O)=O"});
//	}
	
	
	public static int lipinskiViolations(IAtomContainer molecule, double mass, double xlogp, int hbaCount, int hbdCount) {
		// Tice Rules; https://doi.org/10.1002/1526-4998(200101)57:1<3::AID-PS269>3.0.CO;2-6
		
		int lViolations = 0;
		if(mass>500.0) {
			lViolations++;
		}
		// CDK does not offer Moriguchi's mlog calculation capabilities. 
		// R Guha recommends using XLogP instead of ALogP (http://blog.rguha.net/?p=896).
		if(xlogp>5.0) {
			lViolations++;
		}
		if(hbaCount>10) {
			lViolations++;
		}
		if(hbdCount>5) {
			lViolations++;
		}
		return lViolations;
	}

	public static int postEmmergencyHerbicideLikenessViolations(IAtomContainer molecule, double mass, double alogp, int hbaCount, int hbdCount, int rbCount) {
		// Tice Rules; https://doi.org/10.1002/1526-4998(200101)57:1<3::AID-PS269>3.0.CO;2-6
		int iViolations = 0;
		if(mass<150.0 && mass>500.0) {
			iViolations++;
		}
		if(alogp>5.0) {
			iViolations++;
		}
		if(hbaCount<2 && hbaCount>12) {
			iViolations++;
		}
		if(hbdCount>3) {
			iViolations++;
		}		
		if(rbCount>12) {
			iViolations++;
		}
		
		return iViolations;
	}
	
	public static int insecticideLikenessViolations(IAtomContainer molecule, double mass, double alogp, int hbaCount, int hbdCount, int rbCount) {
		// Tice Rules; https://doi.org/10.1002/1526-4998(200101)57:1<3::AID-PS269>3.0.CO;2-6
		int pehViolations = 0;
		
		if(mass<150.0 && mass>500.0) {
			pehViolations++;
		}
		if(alogp<0.0 && alogp>6.5) {
			pehViolations++;
		}
		if(hbaCount<1 && hbaCount>8) {
			pehViolations++;
		}
		if(hbdCount>2) {
			pehViolations++;
		}		
		if(rbCount>12.0) {
			pehViolations++;
		}		
		return pehViolations;
	}	
	
	public static int mddrDrugLikenessViolations(IAtomContainer molecule, int rigidBondCount, int ringCount, int rbCount) {
		//  Oprea,T.I. J. Comput. Aid. Mol. Des. 2000, 14, 251.
		int mddrdliolations = 0;
		
		if(rigidBondCount<3) {
			mddrdliolations++;
		}
		if(rigidBondCount<18) {
			mddrdliolations++;
		}
		if(rbCount<6) {
			mddrdliolations++;
		}
	
		return mddrdliolations;
	}	
	public static LinkedHashMap<String,Integer> calculateLikenessViolations(IAtomContainer molecule) throws CDKException {
		LinkedHashMap<String, Integer> violations = new LinkedHashMap<String, Integer> ();
		String mass = molecule.getProperty("Major Isotope Mass");
		String hbaCount = molecule.getProperty("hbaCount");
		String hbdCount = molecule.getProperty("hbdCount");
		String rbCount  = molecule.getProperty("rbCount");
		String alogp  	= molecule.getProperty("ALogP") ;
		String xlogp	= molecule.getProperty("XLogP");

		double mass_;
		if(Objects.equal(mass, null)) {
			mass_ =  Double.valueOf(getMajorIsotopeMass(molecule));
		}else {
			mass_ = Double.valueOf(mass);
		}

		
		if(Objects.equal(alogp, null)) {
			ALOGPDescriptor aLogpDescriptor = new ALOGPDescriptor();
			alogp = aLogpDescriptor.calculate(molecule).getValue().toString().split(",")[0];
		}		
		if(Objects.equal(xlogp, null)) {
			xlogp = xLogpDescriptor.calculate(molecule).getValue().toString();
		}	
		if(Objects.equal(hbaCount, null)) {
			hbaCount = hbaDCountDescriptor.calculate(molecule).getValue().toString();
//			System.err.println("hbaCount: new");
		}		
		if(Objects.equal(hbdCount, null)) {
			hbdCount = hbdCountDescriptor.calculate(molecule).getValue().toString();
//			System.err.println("hbdCount: new");
		}	
		if(Objects.equal(rbCount, null)) {
			rbCount = rbCountDescriptor.calculate(molecule).getValue().toString();
//			System.err.println("rbCount: new");
		}		
//		System.err.println("hbaCount: "+ hbaCount);
////		System.err.println("int hbaCount: "+ Integer.valueOf(hbaCount));
//		System.err.println("hbdCount: "+ hbdCount);
//		System.err.println("rbCount: "+ rbCount);
		violations.put("Lipinski_Violations", lipinskiViolations(molecule, mass_, Double.valueOf(xlogp), Integer.valueOf(hbaCount), Integer.valueOf(hbdCount)));
		violations.put("Insecticide_Likeness_Violations", insecticideLikenessViolations(molecule, mass_, Double.valueOf(alogp), Integer.valueOf(hbaCount), Integer.valueOf(hbdCount),  Integer.valueOf(rbCount)));
		violations.put("Post_Em_Herbicide_Likeness_Violations", insecticideLikenessViolations(molecule, mass_, Double.valueOf(alogp), Integer.valueOf(hbaCount), Integer.valueOf(hbdCount),  Integer.valueOf(rbCount)));
		
		return violations;
	}
	
	
//	public int agLikenessViolations(IAtomContainer molecule) {
//		int aViolations = 0;
//		
//		return aViolations;
//	}		
}
