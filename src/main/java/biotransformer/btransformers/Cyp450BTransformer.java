/**
 * This class implements the class of CYP450Biotransformers, which simulate the transformation of molecules by CYP450 enzymes.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.btransformers;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import bioTransformerAPI.BioTransformerAPI;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.Utilities;
import cyProduct.cyProductMain;
import exception.BioTransformerException;
import reactantpredictor.BioTransformerAPIs;


public class Cyp450BTransformer extends Biotransformer {
	//public final String[] cyp450Enzymes =  {"CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4"};
	public final String[] cyp450Enzymes =  {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"};
	public Cyp450BTransformer(BioSystemName bioSName) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {
		super(bioSName);
		setCyp450EnzymesAndReactionList();
		
	}
	
	private void setCyp450EnzymesAndReactionList(){
		ArrayList<MetabolicReaction> react = new ArrayList<MetabolicReaction>();
		
		/*
		 * Adding enzymes
		 */
		for(Enzyme enz : this.bSystem.getEnzymesList()){
			if(enz.getName().contains("CYP")){
				this.enzymesList.add(enz);
				react.addAll(enz.getReactionSet());
			}	
		}
		
		/*
		 * adding a list of unique Cyp450 mediated reactions
		 */
		Set<MetabolicReaction> hs = new HashSet<MetabolicReaction>();
		hs.addAll(react);
		react.clear();
		react.addAll(hs);
//		System.out.println("react: " + react.size());
		this.reactionsByGroups.put("cypReactions",react);	
		
		for(MetabolicReaction reaction : react){
			this.reactionsHash.put(reaction.getReactionName(), reaction);
		}
	}

	/**
	 * Added by Yannick Djoumbou Feunang
	 * This function allows the use to select 1 out of the 3 available modes
	 * If the input molecule is a cyp450 substrate, then continue
	 * Otherwise it returns a empty set
	 * @param substrates
	 * @param preprocess
	 * @param filter
	 * @param threshold
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> predictCyp450BiotransformationsByMode(IAtomContainerSet substrates, int mode, boolean preprocess, boolean filter, double threshold)  throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer substrate : substrates.atomContainers()){
			biotransformations.addAll(predictCyp450BiotransformationsByMode(substrate, mode, preprocess, filter, threshold) );					
		}		
		return biotransformations;
	}
	
	
	/**
	 * Added by Siyang
	 * This function allows the use to select 1 out of the 3 available modes
	 * If the input molecule is a cyp450 substrate, then continue
	 * Otherwise it returns a empty set
	 * @param substrate
	 * @param preprocess
	 * @param filter
	 * @param threshold
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> predictCyp450BiotransformationsByMode(IAtomContainer substrate, int mode, boolean preprocess, boolean filter, double threshold)  throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		BioTransformerAPIs cypReact = new BioTransformerAPIs();
		ArrayList<Integer> substrateForEnzymes = cypReact.predictReactantForEnzymes(substrate, this.cyp450Enzymes);
		if(substrateForEnzymes.isEmpty()){
			return biotransformations;
		}		
		if(mode ==  1){
			biotransformations = this.predictCyp450Biotransformations(substrate, false, preprocess, filter, threshold);
		}
		else if(mode == 2){
			biotransformations = this.predictCypProductBiotransformation(substrate, substrateForEnzymes, false);
		}
		else if(mode == 3){			
			biotransformations = this.predictByCyProductAndCyp450(substrate, preprocess, filter,  threshold, substrateForEnzymes);
		}
		else throw new Exception("You have selected the predictCyp450BiotransformationsByMode module, please input the mode you  want to run");
		return biotransformations;
	}
	/**
	 * Added by Siyang
	 * This function will use CyProduct to predict biotransformations
	 * @param substrate
	 * @param useCypReact
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> predictCypProductBiotransformation(IAtomContainer substrate, ArrayList<Integer> enzymeIdx, boolean useCypReact) throws  Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		AtomContainerManipulator.suppressHydrogens(substrate);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(substrate);
		//System.out.println(this.smiGen.create(substrate + " starts");
		IAtomContainer oneMole = cyProductMain.readSMILES_oneMole(this.smiGen.create(substrate));
		oneMole.setProperties(substrate.getProperties());
		
		//System.out.println(this.smiGen.create(substrate + " is done");
		IAtomContainerSet removedNonsense = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<String> enzymeList = new ArrayList<>();
		for(int i = 0; i < enzymeIdx.size(); i++){
			enzymeList.add(this.cyp450Enzymes[enzymeIdx.get(i)]);
		}
		IAtomContainerSet cyProduct_results = BioTransformerAPI.runOnePrediction(oneMole, enzymeList, true);
		for(int i = 0; i < cyProduct_results.getAtomContainerCount(); i++){
			if(cyProduct_results.getAtomContainer(i).getAtomCount() > 1){
				removedNonsense.addAtomContainer(cyProduct_results.getAtomContainer(i));
				String inChiKey = this.inchiGenFactory.getInChIGenerator(cyProduct_results.getAtomContainer(i)).getInchiKey();
				cyProduct_results.getAtomContainer(i).setProperty("InChIKey", inChiKey);
				cyProduct_results.getAtomContainer(i).setProperty("SMILES", this.smiGen.create(cyProduct_results.getAtomContainer(i)));
				cyProduct_results.getAtomContainer(i).setProperty("InChI", this.inchiGenFactory.getInChIGenerator(cyProduct_results.getAtomContainer(i)).getInchi());
			}
			
		}
		cyProduct_results = removedNonsense;
		IAtomContainerSet substrates_temp = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);			
		substrate.setProperty("SMILES", this.smiGen.create(substrate));
		substrate.setProperty("InChI", this.inchiGenFactory.getInChIGenerator(substrate).getInchi());
		substrate.setProperty("InChIKey", this.inchiGenFactory.getInChIGenerator(substrate).getInchiKey());
		substrates_temp.addAtomContainer(substrate);
		ArrayList<Biotransformation> cyProduct_extra_Biotransformations = convertCyProductToBioTransformation(substrates_temp, cyProduct_results);
		biotransformations.addAll(cyProduct_extra_Biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	/**
	 * Added by Siyang
	 * Convert the cyProduct predicted Results to BioTransformations
	 * @param containers
	 * @param nrOfSteps
	 * @param scoreThreshold
	 * @param outputFolder
	 * @param annotate
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> convertCyProductToBioTransformation(IAtomContainerSet substrates, IAtomContainerSet cyProductMolecules) throws Exception{
		//Biotransformation(IAtomContainerSet substrates, String reactionType, ArrayList<String> enzymeNames, 
		//IAtomContainerSet products, BioSystemName bsysName)
		ArrayList<Biotransformation> convertedBioTrans = new ArrayList<>();
		for(int i = 0; i < cyProductMolecules.getAtomContainerCount(); i++){
			IAtomContainerSet notIncludedMolecuels = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			IAtomContainer oneMolecule = cyProductMolecules.getAtomContainer(i);
			notIncludedMolecuels.addAtomContainer(oneMolecule);
			String reactionType = oneMolecule.getProperty("ReactionType");
			String enzymeNameString = oneMolecule.getProperty("Enzyme");
			String[] enzymeArray = enzymeNameString.split(" ");
			ArrayList<String> enzymeNames = new ArrayList<>();
			for(int j = 0; j < enzymeArray.length; j++){
				enzymeNames.add(enzymeArray[j]);
			}
			BioSystemName bioSys = BioSystemName.HUMAN;
			Biotransformation bioTrans = new Biotransformation(substrates, reactionType, enzymeNames, notIncludedMolecuels, bioSys);
			convertedBioTrans.add(bioTrans);
		}
		return convertedBioTrans;
	}
	
	/**
	 * Modified by Siyang. 
	 * Apply CyProduct here. This function will predict the union of BioTransformer cyp450 prediction and CyProduct
	 * @param substrates
	 * @param preprocess
	 * @param filter
	 * @param nrOfSteps
	 * @param threshold
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> predictByCyProductAndCyp450(IAtomContainer substrate, boolean preprocess, boolean filter, double threshold, ArrayList<Integer> enzymeIdx) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		containers.addAtomContainer(substrate);
		ArrayList<Biotransformation> currentBiotransformations = this.predictCyp450Biotransformations(containers, false, preprocess, filter, threshold);
		//ArrayList<Biotransformation> currentBiotransformations = this.predictCyp450Biotransformations(containers, preprocess, filter, threshold);
		if(!currentBiotransformations.isEmpty()){
			biotransformations.addAll(currentBiotransformations);
			containers.removeAllAtomContainers();
			containers = extractProductsFromBiotransformations(currentBiotransformations);				
		}
		ArrayList<String> enzymeList = new ArrayList<>();
		for(int i = 0; i < enzymeIdx.size(); i++){
			enzymeList.add(this.cyp450Enzymes[enzymeIdx.get(i)]);
		}
		AtomContainerManipulator.suppressHydrogens(substrate);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(substrate);
		IAtomContainer oneMole = cyProductMain.readSMILES_oneMole(this.smiGen.create(substrate));
		oneMole.setProperties(substrate.getProperties());
		IAtomContainerSet cyProduct_results = BioTransformerAPI.runOnePrediction(oneMole, enzymeList, false);
		//IAtomContainerSet cyProduct_results = BioTransformerAPI.runOnePrediction(substrates.getAtomContainer(t), enzymeList);
		IAtomContainerSet removedNonsense = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < cyProduct_results.getAtomContainerCount(); i++){
			if(cyProduct_results.getAtomContainer(i).getAtomCount() > 1){
				removedNonsense.addAtomContainer(cyProduct_results.getAtomContainer(i));
				String inChiKey = this.inchiGenFactory.getInChIGenerator(cyProduct_results.getAtomContainer(i)).getInchiKey();
				cyProduct_results.getAtomContainer(i).setProperty("InChIKey", inChiKey);
				cyProduct_results.getAtomContainer(i).setProperty("SMILES", this.smiGen.create(cyProduct_results.getAtomContainer(i)));
				cyProduct_results.getAtomContainer(i).setProperty("InChI", this.inchiGenFactory.getInChIGenerator(cyProduct_results.getAtomContainer(i)).getInchi());
			}
			
		}
		cyProduct_results = removedNonsense;
		IAtomContainerSet bioTransformerResults = extractProductsFromBiotransformations(biotransformations);	
		IAtomContainerSet substrates_temp = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);			
		substrate.setProperty("SMILES", this.smiGen.create(substrate));
		substrate.setProperty("InChI", this.inchiGenFactory.getInChIGenerator(substrate).getInchi());
		substrate.setProperty("InChIKey", this.inchiGenFactory.getInChIGenerator(substrate).getInchiKey());
		substrates_temp.addAtomContainer(substrate);
		IAtomContainerSet notIncluded_results = notIncludedMolecules(substrates_temp, cyProduct_results, bioTransformerResults);
		bioTransformerResults.add(notIncluded_results);
		ArrayList<Biotransformation> cyProduct_extra_Biotransformations = convertCyProductToBioTransformation(substrates_temp, notIncluded_results);
		biotransformations.addAll(cyProduct_extra_Biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	/**
	 * Added by Siyang
	 * Find the metabolites in set A and not in set B
	 * The input are two IAtomcontainerSets, one is set A and the other is set B
	 * @param containers
	 * @param nrOfSteps
	 * @param scoreThreshold
	 * @param outputFolder
	 * @param annotate
	 * @throws Exception
	 */
	public IAtomContainerSet notIncludedMolecules(IAtomContainerSet substrates, IAtomContainerSet cyProductResults, IAtomContainerSet cyp450BioReuslts) throws Exception{
		//ArrayList<String> inChiKeysOfCyp450 = new ArrayList<>();
		IAtomContainerSet notIncludedMolecuels = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i <cyProductResults.getAtomContainerCount(); i++){
			if(!ChemStructureExplorer.atomContainerInclusionHolds(cyp450BioReuslts, cyProductResults.getAtomContainer(i))){
				notIncludedMolecuels.addAtomContainer(cyProductResults.getAtomContainer(i));
			}
		}
//		for(int i = 0; i < cyp450BioReuslts.getAtomContainerCount(); i++){
//			//atomContainerInclusionHolds(cyp450BioReuslts, )
//			String identifier = this.inchiGenFactory.getInChIGenerator(cyp450BioReuslts.getAtomContainer(i)).getInchiKey().split("-")[0];
//			if(!inChiKeysOfCyp450.contains(identifier)) inChiKeysOfCyp450.add(identifier);
//		}
//		
//		for(int i = 0; i < cyProductResults.getAtomContainerCount(); i++){
//			String identifier = this.inchiGenFactory.getInChIGenerator(cyProductResults.getAtomContainer(i)).getInchiKey().split("-")[0];
//			if(!inChiKeysOfCyp450.contains(identifier)){
//				notIncludedMolecuels.addAtomContainer(cyProductResults.getAtomContainer(i));
//				inChiKeysOfCyp450.add(identifier);
//			}
//		}
		return notIncludedMolecuels;
		
	}
	
	public ArrayList<Biotransformation> predictCyp450Biotransformations(IAtomContainer substrate, boolean predictESSpecificity, boolean preprocess, boolean filter, double threshold) throws Exception{
		
		try{
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
			
			ChemStructureExplorer.addInChIandKey(substrate);
			if(ChemStructureExplorer.isCompoundInorganic(substrate) || ChemStructureExplorer.isMixture(substrate)){
				throw new IllegalArgumentException(substrate.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(substrate)){
//				ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(substrate);
//				System.err.println("will compute chemClasses.");
//				
//				if (ChemStructureExplorer.getMajorIsotopeMass(substrate) > 1500.0 || 
//						chemClasses.contains(ChemicalClassName.GLYCOSYLATED_COMPOUND) ||
//						chemClasses.contains(ChemicalClassName.GLUTATHIONE_CONJUGATE) ||
//						chemClasses.contains(ChemicalClassName.SULFATED_COMPOUND) ||
//						chemClasses.contains(ChemicalClassName.ACYL_CoA_CONJUGATE) ||
//						chemClasses.contains(ChemicalClassName.TETRAPYRROLE) ||
//						chemClasses.contains(ChemicalClassName.SACCHARIDE) ||
//						chemClasses.contains(ChemicalClassName.ETHER_LIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROLIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROPHOSPHOLIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL) ||
//						chemClasses.contains(ChemicalClassName.SPHINGOLIPID)){
//					
//				}
//				else{
					String[] cyp450s =  {"CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", 
							"CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4"};
					
					ArrayList<Enzyme> cyp450Enzymes = new ArrayList<Enzyme>();
					for(int i = 0; i < cyp450s.length; i++){
						cyp450Enzymes.add(this.bSystem.getEnzymeHash().get(cyp450s[i]));
					}
					
					biotransformations = this.metabolizeWithEnzymes(substrate, cyp450Enzymes, predictESSpecificity, preprocess, filter, threshold);
	//				for(String en : cyp450s) { 
	////					System.out.println("Predicting metabolism for " + en.name());
	//					biotransformations.addAll(this.metabolizeWithEnzyme(substrate, en, preprocess, filter, threshold));
	//				}
//				}			
				
			}
			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	
	public ArrayList<Biotransformation> predictCyp450BiotransformationsForSpecificEnzymes(IAtomContainer substrate, ArrayList<String> enzymeNames, boolean preprocess, boolean filter, double threshold) throws Exception{
		
		try{
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
			
			ChemStructureExplorer.addInChIandKey(substrate);
			if(ChemStructureExplorer.isCompoundInorganic(substrate) || ChemStructureExplorer.isMixture(substrate)){
				throw new IllegalArgumentException(substrate.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} else if(ChemStructureExplorer.isBioTransformerValid(substrate)) {
//				ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(substrate);
				
//				if (ChemStructureExplorer.getMajorIsotopeMass(substrate) > 1500.0 || 
//						chemClasses.contains(ChemicalClassName.GLYCOSYLATED_COMPOUND) ||
//						chemClasses.contains(ChemicalClassName.GLUTATHIONE_CONJUGATE) ||
//						chemClasses.contains(ChemicalClassName.SULFATED_COMPOUND) ||
//						chemClasses.contains(ChemicalClassName.ACYL_CoA_CONJUGATE) ||
//						chemClasses.contains(ChemicalClassName.TETRAPYRROLE) ||
//						chemClasses.contains(ChemicalClassName.SACCHARIDE) ||
//						chemClasses.contains(ChemicalClassName.ETHER_LIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROLIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROPHOSPHOLIPID) ||
//						chemClasses.contains(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL) ||
//						chemClasses.contains(ChemicalClassName.SPHINGOLIPID)){	
//				}
//				else{
					ArrayList<Enzyme> cyp450Enzymes = new ArrayList<Enzyme>();
					for(int i = 0; i < enzymeNames.size(); i++){
						cyp450Enzymes.add(this.bSystem.getEnzymeHash().get(enzymeNames.get(i)));
					}
					
					biotransformations = this.metabolizeWithEnzymes(substrate, cyp450Enzymes, preprocess, filter, threshold);
	
//				}			
				
			}
			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	public ArrayList<Biotransformation> predictCyp450Biotransformations(IAtomContainerSet substrates, boolean preprocess, boolean filter, double threshold) throws Exception{
		return predictCyp450Biotransformations(substrates, true, preprocess, filter, threshold);
	}	

	public ArrayList<Biotransformation> predictCyp450Biotransformations(IAtomContainerSet substrates, boolean predictESSpecificity, boolean preprocess, boolean filter, double threshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		for(IAtomContainer sub : substrates.atomContainers()){
			biotransformations.addAll(predictCyp450Biotransformations(sub, predictESSpecificity, preprocess, filter, threshold) );					
		}
		
		return biotransformations;
	}
	
	public ArrayList<Biotransformation> predictCyp450Biotransformations(String sdfFileName, boolean preprocess, boolean filter, double threshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		IAtomContainerSet containers = FileUtilities.parseSdf(sdfFileName);
		biotransformations = predictCyp450Biotransformations(containers, preprocess, filter, threshold);

		return biotransformations;
	}
	
	
	public ArrayList<Biotransformation> predictCyp450BiotransformationChain(IAtomContainer substrate, boolean preprocess, boolean filter, int nrOfSteps, double threshold) throws Exception{
		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		molecules.addAtomContainer(substrate);
		return this.predictCyp450BiotransformationChain(molecules, preprocess, filter, nrOfSteps, threshold);
	}
	
	public ArrayList<Biotransformation> predictCyp450BiotransformationChain(IAtomContainerSet substrates, boolean preprocess, boolean filter, int nrOfSteps, double threshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = substrates;
		int counter = 0;
		while(nrOfSteps>0){
			counter++;
			ArrayList<Biotransformation> currentBiotransformations = this.predictCyp450Biotransformations(containers, preprocess, filter, threshold);
			nrOfSteps--;
			if(!currentBiotransformations.isEmpty()){
				biotransformations.addAll(currentBiotransformations);
				containers.removeAllAtomContainers();
				containers = extractProductsFromBiotransformations(currentBiotransformations);				
			}
			else{
				break;
			}
		}

		return Utilities.selectUniqueBiotransformations(biotransformations);
	}	

	
	public ArrayList<Biotransformation> predictCyp450BiotransformationChainByMode(IAtomContainer substrate, boolean preprocess, boolean filter, int nrOfSteps, double threshold, int cyp450Mode) throws Exception{
		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		molecules.addAtomContainer(substrate);		
		return predictCyp450BiotransformationChainByMode(molecules, preprocess, filter, nrOfSteps, threshold, cyp450Mode);
	}
	
	public ArrayList<Biotransformation> predictCyp450BiotransformationChainByMode(IAtomContainerSet substrates, boolean preprocess, boolean filter, int nrOfSteps, double threshold, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = substrates;
		int counter = 0;
		while(nrOfSteps>0){
			counter++;
			ArrayList<Biotransformation> currentBiotransformations = this.predictCyp450BiotransformationsByMode(containers, cyp450Mode, preprocess, filter, threshold);
			nrOfSteps--;
			if(!currentBiotransformations.isEmpty()){
				biotransformations.addAll(currentBiotransformations);
				containers.removeAllAtomContainers();
				containers = extractProductsFromBiotransformations(currentBiotransformations);				
			}
			else{
				break;
			}
		}

		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	
	public void  simulateCyp450MetabolismAndSaveToSDF(IAtomContainerSet containers, int nrOfSteps, Double scoreThreshold, String outputFolder, boolean annotate) throws Exception {
		
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){

				String identifier = molecule.getProperty(CDKConstants.TITLE);
				if(identifier == null){
					identifier = molecule.getProperty("Name");
					if(identifier == null){
						identifier = molecule.getProperty("InChiKey");
						if(identifier == null){
							identifier = this.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
						}
					}
				}
				System.out.println(identifier);
				if(identifier.contains("/") || identifier.contains(":")){
					System.out.println("The identifier contains characters that cannot be used to create a file. / and : would be replaced wit - and _, respectively");
					identifier = identifier.replace(":", "-").replace("/", "_");
				}
				ArrayList<Biotransformation> biotransformations = this.predictCyp450BiotransformationChain(molecule, true, true, nrOfSteps, scoreThreshold);
//				System.out.println(biotransformations.size() + " biotransformations.");
				this.saveBioTransformationProductsToSdf(biotransformations, outputFolder + "/" + identifier + "_CYP450_based_metabolites.sdf", annotate);			
			}
		}		
	}
	
	
	public void  simulateCyp450MetabolismAndSaveToSingleSDF(IAtomContainerSet containers, int nrOfSteps, Double scoreThreshold, String outputFileName, boolean annotate) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){
				biotransformations.addAll(this.predictCyp450BiotransformationChain(molecule, true, true, nrOfSteps, scoreThreshold));
//				System.out.println(biotransformations.size() + " biotransformations.");
			}
			
			this.saveBioTransformationProductsToSdf(biotransformations, outputFileName, annotate);			

		}		
	}
	
	public void printStatistics(){
//		System.out.println(this.reactionsByGroups.get("cypReactions").size());
		int count = 0;
		for(Enzyme e: this.enzymesList){
//			System.out.println(e.getName() + " : " + e.getReactionsNames().size());
			count = count + e.getReactionsNames().size();
		}
		System.out.println("Humber of enzymes: " + this.enzymesList.size());
		System.out.println("Humber of biotransformation rules: " + this.reactionsByGroups.get("cypReactions").size());
		System.out.println("Humber of enzyme-biotransformation rules associations: " + count);
	}
	
}
