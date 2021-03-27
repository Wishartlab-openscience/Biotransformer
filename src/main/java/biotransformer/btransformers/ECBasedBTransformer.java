/**
 * This class implements the class of EC-Based biotransformers, which simulate the transformation of molecules
 * by enzymes, based the NC-IUBMB's EC-based classification.
 * 
 * NC-IUBMB : Nomenclature Committee of the International Union of Biochemistry and Molecular Biology
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */
package biotransformer.btransformers;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import ambit2.smarts.query.SMARTSException;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemicalClassFinder;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;
import biotransformer.utils.FileUtilities;
import exception.BioTransformerException;

/**
 * @author Djoumbou Feunang, Yannick
 *
 */
public class ECBasedBTransformer extends Biotransformer {

	/**
	 * @throws ParseException 
	 * @throws IOException 
	 * @throws CDKException 
	 * @throws URISyntaxException 
	 * @bioSName Name of the biosystem to simulate.
	 * 
	 */
	public ECBasedBTransformer(BioSystemName bioSName) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {		
		// TODO Auto-generated constructor stub
		super(bioSName);
		setEnzymesList();
		setReactionsList();

	}
	
	/**
	 * Collect the list of enzymes associated with the current biosystem from the database.
	 */
	
	private void setEnzymesList(){
		for(Enzyme enz : this.bSystem.getEnzymesList()){
			if(enz.getName().contains("EC_") || enz.getName().contains("UDP_GLUCURONOSYLTRANSFERASE") || 
				enz.getName().contains("SULFOTRANSFERASE") || enz.getName().contains("ACETYLTRANSFERASE") ||
				enz.getName().contentEquals("ABKAR1") || enz.getName().contentEquals("ABKAR2") || enz.getName().contentEquals("CYPB5_CYPBR5")){
				if(!enz.getReactionSet().isEmpty()){
					this.enzymesList.add(enz);
				}
				
			}
		}
	}
	
	/**
	 * Collect the list of metabolic reactions associated with the current biosystem, inferred from the list of enzymes.
	 */
	
	private void setReactionsList(){
//		this.reactionsList.put("standardizationReactions", MReactionSets.standardizationReactions);
		this.enzymesByreactionGroups.put("ecBasedDeconjugations", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("ecBasedReactionsNonDeconjugative", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("ecBasedConjugations", new ArrayList<Enzyme>());
		this.reactionsByGroups.put("ecBasedDeconjugations", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("ecBasedReactionsNonDeconjugative", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("ecBasedConjugations", new ArrayList<MetabolicReaction>());
		LinkedHashMap<String, MetabolicReaction> tempR = new LinkedHashMap<String, MetabolicReaction>();
		ArrayList<String> rNames =  new ArrayList<String>();
		for(Enzyme enz : this.enzymesList){
//			System.err.println(enz.getName());
			if (
					// carboxylesterase
					enz.getName().contentEquals("EC_3_1_1_1") || enz.getName().contentEquals("EC_3_1_1_20") || 
					
					// sulfatase
					
					enz.getName().contentEquals("EC_3_1_6_1") || 
					
					// glycosidase
					enz.getName().contentEquals("EC_3_2_1_X") || enz.getName().contentEquals("EC_3_2_1_20") || 
					enz.getName().contentEquals("EC_3_2_1_21") || enz.getName().contentEquals("EC_3_2_1_23") || 
					enz.getName().contentEquals("EC_3_2_1_31") || enz.getName().contentEquals("EC_3_2_1_40") ||
					enz.getName().contentEquals("EC_3_2_1_147") )
					
				
			{
				this.enzymesByreactionGroups.get("ecBasedDeconjugations").add(enz);
				for(MetabolicReaction m : enz.getReactionSet()){
					if(!rNames.contains(m.name)){
						rNames.add(m.name);
						this.reactionsByGroups.get("ecBasedDeconjugations").add(m);
						this.reactionsHash.put(m.name, m);
					}
				}
			} 
			else if(
					enz.getName().contains("EC_2_1_1") || 
					enz.getName().contains("EC_2_3_1_5") ||
					enz.getName().contains("EC_2_3_1_13") || 
					enz.getName().contains("EC_2_3_1_65") ||
					enz.getName().contains("EC_2_5_1_18") || 
					enz.getName().contains("EC_2_3_1_80") || 
					enz.getName().contains("EC_2_3_1_85") || 
					enz.getName().contains("EC_2_3_1_87") || 
					enz.getName().contains("EC_2_3_1_88") || 
					enz.getName().contains("EC_2_3_1_118") || 
					enz.getName().contains("EC_2_4_1_17") || 
//					enz.getName().contains("") || 
//					enz.getName().contains("") || 
//					enz.getName().contains("") || 
					
					enz.getName().contains("EC_2_8_2") ||
					enz.getName().contentEquals("SULFOTRANSFERASE") ||
					enz.getName().contentEquals("UDP_GLUCURONOSYLTRANSFERASE")
					
					){
				
				this.enzymesByreactionGroups.get("ecBasedConjugations").add(enz);
				for(MetabolicReaction m : enz.getReactionSet()){
					if(!rNames.contains(m.name)){
						rNames.add(m.name);
						this.reactionsByGroups.get("ecBasedConjugations").add(m);
						this.reactionsHash.put(m.name, m);
					}
				}
			}
			else {
//				System.out.println(enz.getName());
				this.enzymesByreactionGroups.get("ecBasedReactionsNonDeconjugative").add(enz);
				for(MetabolicReaction m : enz.getReactionSet()){
					if(!rNames.contains(m.name)){
//						System.err.println(m.name);
						rNames.add(m.name);
						this.reactionsByGroups.get("ecBasedReactionsNonDeconjugative").add(m);
						this.reactionsHash.put(m.name, m);
					}
				}
			}
			
		}
	}
	
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nrOfSteps -  number of steps
	 * @return an arraylist of biotransformations, which are instances of the EC-based metabolic reactions applied to the target, with a threshold of 0.0
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyEcBasedDeconjuations(IAtomContainer target, boolean preprocess, boolean filter, int nrOfSteps)
			throws Exception {
		return  applyEcBasedDeconjugations(target, preprocess, filter, nrOfSteps, 0.0);
	}

	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param scoreThreshold - minimum score threshold for selected metbaolic reactions
	 * @return an arraylist of biotransformations, which are instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> applyEcBasedDeconjugations(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
//		ArrayList<Biotransformation> biotransformations = applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.
//				get("ecBasedDeconjugations"), preprocess, filter, scoreThreshold);
		
		
		try {
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)){
				
				biotransformations = metabolizeWithEnzymesBreadthFirst(target,
							this.enzymesByreactionGroups.get("ecBasedDeconjugations"), preprocess, filter, nr_of_steps, scoreThreshold);			
			}

			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception - Throw an exception
	 */
	public ArrayList<Biotransformation> applyEcBasedDeconjuationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps)
			throws Exception {

		return  applyEcBasedDeconjuationsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */

	public ArrayList<Biotransformation> applyEcBasedDeconjuationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
//		ArrayList<Biotransformation> biotransformations = applyReactionsChainAndReturnBiotransformations(target, this.reactionsByGroups.get("ecBasedDeconjugations"), preprocess, filter, nr_of_steps, scoreThreshold);		

		
		try{
			ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();	
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)){
				ArrayList<Biotransformation> biotransformations = metabolizeWithEnzymesBreadthFirst(target,
						this.enzymesByreactionGroups.get("ecBasedDeconjugations"), preprocess, filter, nr_of_steps, scoreThreshold);
	
				int btCounter = 0;
				
				for(Biotransformation bt : biotransformations) {
					boolean goodBiontransfo =  true;
					btCounter++;
	//				System.out.println("Biotransformation nr. " +  btCounter);
	//				System.out.println("Biotransformation type. " +  bt.getReactionType());
					for( IAtomContainer at : bt.getProducts().atomContainers() ){
	//					System.out.println(at.getProperty("InChI"));
						if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
	//						System.err.println("is invalid phase 2 metabolite");
							goodBiontransfo = false;					
							break;
						}
										
					}
					
	//				System.out.println("goodBiontransfo: " + goodBiontransfo);
					if(goodBiontransfo){
						selectedBiotransformations.add(bt);
	//					System.err.println("Biotransformation " + btCounter + " has been added.");
					}
	//				System.err.println("selected biotransformations: " + selectedBiotransformations.size());
				}
				
			}

		
		return selectedBiotransformations;
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the EC-based metabolic reactions applied to the target, with a threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyEcBasedTransformations(IAtomContainer target, boolean preprocess, boolean filter)
			throws Exception {
		return  applyEcBasedTransformations(target, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param scoreThreshold - minimum threshold for reaction scores
	 * @return an arraylist of biotransformations, which are instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> applyEcBasedTransformations(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		
		
		try{
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)){
				if(!this.isDeconjugationCandidate(target)){
	//				System.out.println("Applying EC for " + this.smiGen.create(target));
	//				biotransformations = applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.
	//						get("ecBasedReactionsNonDeconjugative"), preprocess, filter, scoreThreshold);
					
					biotransformations = metabolizeWithEnzymes(target,
							this.enzymesByreactionGroups.get("ecBasedReactionsNonDeconjugative"), preprocess, filter, scoreThreshold);
					int btCounter = 0;
					ArrayList<Biotransformation> filteredBiotransformations = (ArrayList<Biotransformation>) biotransformations.clone();
					for(Biotransformation bt : biotransformations){
						btCounter++;
	//					System.out.println("Biotransformation nr. " +  btCounter);
						for(IAtomContainer at : bt.getProducts().atomContainers() ){
	//						System.out.println(at.getProperty("InChI"));
							if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
	//							System.err.println("is invalid phase 2 metabolite");
	//							biotransformations.remove(bt);
								filteredBiotransformations.remove(bt);
								break;
							}
						}	
					}
				}		
			}
			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyEcBasedTransformationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps)
			throws Exception {

		return  applyEcBasedTransformationsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyEcBasedTransformationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
			
		try {
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target) && !this.isDeconjugationCandidate(target)){
				int btCounter = 0;
	//			biotransformations = applyReactionsChainAndReturnBiotransformations(target, this.reactionsByGroups.get("ecBasedReactionsNonDeconjugative"), preprocess, filter, nr_of_steps, scoreThreshold);		
				biotransformations = metabolizeWithEnzymesBreadthFirst(target,
						this.enzymesByreactionGroups.get("ecBasedReactionsNonDeconjugative"), preprocess, filter, nr_of_steps, scoreThreshold);
	//			for(Biotransformation bt : biotransformations) {
	//				boolean goodBiontransfo =  true;
	//				btCounter++;
	//				System.out.println("Biotransformation nr. " +  btCounter);
	//				System.out.println("Biotransformation type. " +  bt.getReactionType());
	//				for( IAtomContainer at : bt.getProducts().atomContainers() ){
	//					System.out.println(at.getProperty("InChI"));
	//					if(ChemStructureExplorer.isInvalidPhaseIIMetabolite(at)){
	//						System.err.println("is invalid phase 2 metabolite");
	//						goodBiontransfo = false;					
	//						break;
	//					}	
	//				}
	//				
	//				System.out.println("goodBiontransfo: " + goodBiontransfo);
	//				if(goodBiontransfo){
	//					selectedBiotransformations.add(bt);
	//					System.err.println("Biotransformation " + btCounter + " has been added.");
	//				}
	//				System.err.println("selected biotransformations: " + selectedBiotransformations.size());
	//			}
			}
	
			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @param scoreThreshold - minimum threshold for reaction scores
	 * @return an arraylist of biotransformations, which are instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> applyEcBasedConjugations(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
//		ArrayList<Biotransformation> biotransformations = applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.
//				get("ecBasedConjugations"), preprocess, filter, scoreThreshold);
		
		try{
			ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();			
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)){
				ArrayList<Biotransformation> biotransformations = metabolizeWithEnzymes(target,
						this.enzymesByreactionGroups.get("ecBasedConjugations"), preprocess, filter, scoreThreshold);
				
				
				int btCounter = 0;
				for(Biotransformation bt : biotransformations) {
					boolean goodBiontransfo =  true;
					btCounter++;
	//				System.out.println("Biotransformation nr. " +  btCounter);
	//				System.out.println("Biotransformation type. " +  bt.getReactionType());
					for( IAtomContainer at : bt.getProducts().atomContainers() ){
	//					System.out.println(at.getProperty("InChI"));
						if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
	//						System.err.println("is invalid phase 2 metabolite");
							goodBiontransfo = false;					
							break;
						}	
					}
					
	//				System.out.println("goodBiontransfo: " + goodBiontransfo);
					if(goodBiontransfo){
						selectedBiotransformations.add(bt);
	//					System.err.println("Biotransformation " + btCounter + " has been added.");
					}
					System.err.println("selected biotransformations: " + selectedBiotransformations.size());
				}
				
			}
			return selectedBiotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyEcBasedConjuationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps)
			throws Exception {

		return  applyEcBasedDeconjuationsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the EC-based metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyEcBasedConjuationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
//		ArrayList<Biotransformation> biotransformations = applyReactionsChainAndReturnBiotransformations(target, this.reactionsByGroups.get("ecBasedConjugations"), preprocess, filter, nr_of_steps, scoreThreshold);				
		
		ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();		
		int btCounter = 0;
		
		if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
			throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
		} 
		else if(ChemStructureExplorer.isBioTransformerValid(target)){
			ArrayList<Biotransformation> biotransformations = metabolizeWithEnzymes(target,
					this.enzymesByreactionGroups.get("ecBasedConjugations"), preprocess, filter, scoreThreshold);

			for(Biotransformation bt : biotransformations) {
				boolean goodBiontransfo =  true;
				btCounter++;
//				System.out.println("Biotransformation nr. " +  btCounter);
//				System.out.println("Biotransformation type. " +  bt.getReactionType());
				for( IAtomContainer at : bt.getProducts().atomContainers() ){
//					System.out.println(at.getProperty("InChI"));
					if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
//						System.err.println("is invalid phase 2 metabolite");
						goodBiontransfo = false;					
						break;
					}
					
					
				}
				
//				System.out.println("goodBiontransfo: " + goodBiontransfo);
				if(goodBiontransfo){
					selectedBiotransformations.add(bt);
//					System.err.println("Biotransformation " + btCounter + " has been added.");
				}
//				System.err.println("selected biotransformations: " + selectedBiotransformations.size());
			}		
		}

		
		return selectedBiotransformations;
	}



	public ArrayList<Biotransformation> simulateECBasedPhaseIMetabolismStep(IAtomContainer molecule, boolean preprocess, boolean filter, Double scoreThreshold) {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();


			try {
//				if(ChemStructureExplorer.isBioTransformerValid(molecule)){
					// Rewrite this when the sChemCategoryPathways has been added & parsed
//			System.out.println("LIPID??");
					ChemicalClassName chemClassName = ChemicalClassFinder.findChemicalClass(molecule);
					if(!(chemClassName == ChemicalClassName.ETHER_LIPID || chemClassName == ChemicalClassName.GLYCEROLIPID || chemClassName == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
							chemClassName == ChemicalClassName.SPHINGOLIPID ||chemClassName == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL || 
							chemClassName == ChemicalClassName.C24_BILE_ACID || chemClassName == ChemicalClassName.C23_BILE_ACID
							)) {
//						System.out.println("NO");
						 if(isDeconjugationCandidate(molecule)){
//					 			System.out.println("IS DECONJUGATION CANDIDATE");
							 biotransformations.addAll(this.applyEcBasedDeconjuationsChain(molecule, true, true, 1, scoreThreshold));
						 }
						 else{
//					 		System.out.println("IS NOT A DECONJUGATION CANDIDATE");
							 biotransformations.addAll(this.applyEcBasedTransformationsChain(molecule, true, true, 1, scoreThreshold));
						 }
						 
					} else if(ChemStructureExplorer.isCompoundInorganic(molecule) || ChemStructureExplorer.isMixture(molecule)){
						throw new IllegalArgumentException(molecule.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
					}
					else if(ChemStructureExplorer.isBioTransformerValid(molecule)) {
						
						
						
						if(chemClassName == ChemicalClassName.ETHER_LIPID){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "ETHER_LIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}
						if(chemClassName == ChemicalClassName.GLYCEROLIPID){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "GLYCEROLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}
						if(chemClassName == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "INOSITOL_PHOSPHATE_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}							
						if(chemClassName == ChemicalClassName.GLYCEROPHOSPHOLIPID){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "GLYCEROPHOSPHOLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}		
						if(chemClassName == ChemicalClassName.SPHINGOLIPID){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "SPHINGOLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}
						if(chemClassName == ChemicalClassName.C24_BILE_ACID || chemClassName == ChemicalClassName.C23_BILE_ACID){
							biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "BILE_ACID_METABOLISM", preprocess, filter, 1, scoreThreshold));
						}
					}		

//				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		return biotransformations;
	}
	
	public ArrayList<Biotransformation> simulateECBasedPhaseIMetabolismStep(IAtomContainerSet molecules, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();		
		for(IAtomContainer ac : molecules.atomContainers()){
			biotransformations.addAll(this.simulateECBasedPhaseIMetabolismStep(ac, preprocess, filter, scoreThreshold));
		}	
		return biotransformations;
	}
	
	public ArrayList<Biotransformation> simulateECBasedPhaseIMetabolismChain(IAtomContainer molecule, boolean preprocess, boolean filter,int nrOfSteps, Double scoreThreshold) throws Exception{
		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		molecules.addAtomContainer(molecule);
		return simulateECBasedPhaseIMetabolismChain(molecules, preprocess, filter, nrOfSteps, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> simulateECBasedPhaseIMetabolismChain(IAtomContainerSet molecules, boolean preprocess, boolean filter,int nrOfSteps, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = molecules;
		int counter = 0;
		while(nrOfSteps>0){
			counter++;
			ArrayList<Biotransformation> currentBiotransformations = this.simulateECBasedPhaseIMetabolismStep(containers, preprocess, filter, scoreThreshold);
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
//		System.out.println("Stopped after " + counter + " steps.");
		return biotransformations;
	}

	
	public ArrayList<Biotransformation> simulateECBasedMetabolismStep(IAtomContainer molecule, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		try{
			if(ChemStructureExplorer.isMixture(molecule) || ChemStructureExplorer.isCompoundInorganic(molecule)){
	//			throw new IllegalArgumentException("The substrate must not be a mixture.");
				throw new IllegalArgumentException(molecule.getProperty("InChI")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
	//			return biotransformations;
			}
			else if(ChemStructureExplorer.isBioTransformerValid(molecule)){
				
				// Rewrite this when the sChemCategoryPathways has been added & parsed
	//			System.out.println("LIPID??");
				
				ChemicalClassName chemClassName = ChemicalClassFinder.findChemicalClass(molecule);
//				System.out.print("Chemical class: ");
//				System.out.println(chemClassName);
				
				if(!(chemClassName == ChemicalClassName.ETHER_LIPID || chemClassName == ChemicalClassName.GLYCEROLIPID || chemClassName == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
						chemClassName == ChemicalClassName.SPHINGOLIPID ||chemClassName == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL ||
						chemClassName == ChemicalClassName.C24_BILE_ACID || chemClassName == ChemicalClassName.C23_BILE_ACID)) {
					 
					 if(isDeconjugationCandidate(molecule)){				 
	
						 biotransformations.addAll(this.applyEcBasedDeconjuationsChain(molecule, true, true, 1, scoreThreshold));
					 }
					 else{
						 biotransformations.addAll(this.applyEcBasedTransformationsChain(molecule, true, true, 1, scoreThreshold));
						 if(this.isConjugationCandidate(molecule)){
							 
							// Inconsistency here. For some enzymes, I must apply the phase 2 filter
							 biotransformations.addAll(this.applyEcBasedConjuationsChain(molecule, true, true, 1, scoreThreshold));
						 }
					 }
					 
				} else {
					
	
	//				if(ChemicalClassFinder.isEtherLipid(molecule)){
					if(chemClassName == ChemicalClassName.ETHER_LIPID){
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "ETHER_LIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}
					if(chemClassName == ChemicalClassName.GLYCEROLIPID){
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "GLYCEROLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}
					if(chemClassName == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL){
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "INOSITOL_PHOSPHATE_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}
					if(chemClassName == ChemicalClassName.GLYCEROPHOSPHOLIPID){
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "GLYCEROPHOSPHOLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}		
					if(chemClassName == ChemicalClassName.SPHINGOLIPID){
	//					System.out.println(chemClassName + " see?");
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "SPHINGOLIPID_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}
					if(chemClassName == ChemicalClassName.C24_BILE_ACID || chemClassName == ChemicalClassName.C23_BILE_ACID){
						biotransformations.addAll(this.applyPathwaySpecificBiotransformationsChain(molecule, "BILE_ACID_METABOLISM", preprocess, filter, 1, scoreThreshold));
					}
				}				
			}
			
			return biotransformations;
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	public ArrayList<Biotransformation> simulateECBasedMetabolismStep(IAtomContainerSet molecules, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();		
		for(IAtomContainer ac : molecules.atomContainers()){
			biotransformations.addAll(this.simulateECBasedMetabolismStep(ac, preprocess, filter, scoreThreshold));
		}	
		return biotransformations;
	}
	
	public ArrayList<Biotransformation> simulateECBasedMetabolismChain(IAtomContainer molecule, boolean preprocess, boolean filter,int nrOfSteps, Double scoreThreshold) throws Exception{
		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		molecules.addAtomContainer(molecule);
		return this.simulateECBasedMetabolismChain(molecules, preprocess, filter, nrOfSteps, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> simulateECBasedMetabolismChain(IAtomContainerSet molecules, boolean preprocess, boolean filter,int nrOfSteps, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = molecules;
		int counter = 0;
		try {
			while(nrOfSteps>0){
				counter++;
				ArrayList<Biotransformation> currentBiotransformations = this.simulateECBasedMetabolismStep(containers, preprocess, filter, scoreThreshold);
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
	//		System.out.println("Stopped after " + counter + " steps.");
			return biotransformations;
		}
		catch(Exception e) {
			throw e;
		}
	}

	
	public void  simulateECBasedMetabolismAndSaveToSDF(IAtomContainerSet containers, int nrOfSteps, Double scoreThreshold, String outputFolder, boolean annotate) throws Exception {	
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){
				ChemStructureExplorer.addInChIandKey(molecule);
				String identifier = molecule.getProperty(CDKConstants.TITLE);
				if(identifier == null){
					identifier = molecule.getProperty("Name");
					if(identifier == null){
						identifier = molecule.getProperty("InChIKey");
					}
				}
				identifier = identifier.replace(":", "-").replace("/", "_");
//				System.out.println(identifier);
				ArrayList<Biotransformation> biotransformations = this.simulateECBasedMetabolismChain(molecule, true, true, nrOfSteps, scoreThreshold);
//				System.out.println(biotransformations.size() + " biotransformations.");
//				this.saveBioTransformationProductsToSdf(biotransformations, outputFolder + "/" + identifier + "_EC_based_metabolites.sdf");
				this.saveBioTransformationProductsToSdf(biotransformations, outputFolder + "/" + identifier + "_EC_based_metabolites.sdf", annotate);			
			}
		}		
	}
	
	
	public void  simulateECBasedMetabolismAndSaveToSDF(String sdfFileName, int nrOfSteps, Double scoreThreshold, String outputFolder, boolean annotate) throws Exception {
		IAtomContainerSet containers = FileUtilities.parseSdf(sdfFileName);
		simulateECBasedMetabolismAndSaveToSDF(containers, nrOfSteps, scoreThreshold, outputFolder, annotate);
	}
	
	public boolean isDeconjugationCandidate(IAtomContainer molecule) throws SMARTSException, CDKException, IOException{
		boolean dec = false;
		
		for(MetabolicReaction m : this.reactionsByGroups.get("ecBasedDeconjugations")){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molecule)){
				dec = true;
				break;
			}
		}
		
		return dec;		
	}
	
	
	public boolean isConjugationCandidate(IAtomContainer molecule) throws SMARTSException, CDKException, IOException{
		boolean dec = false;
		
		for(MetabolicReaction m : this.reactionsByGroups.get("ecBasedConjugations")){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molecule)){
				dec = true;
				break;
			}
		}
		
		return dec;		
	}
	
	
	public void printStatistics(){
		int count = 0;
		for(Enzyme e: this.enzymesList){
			count = count + e.getReactionsNames().size();
		}
		System.out.println("Number of enzymes: " + this.enzymesList.size());
		System.out.println("Number of biotransformation rules: " + (this.reactionsByGroups.get("ecBasedDeconjugations").size()
				+ this.reactionsByGroups.get("ecBasedReactionsNonDeconjugative").size() 
				+ this.reactionsByGroups.get("ecBasedConjugations").size()) );
		
		System.out.println("Number of enzyme-biotransformation rules associations: " + count);
	
	}
}
