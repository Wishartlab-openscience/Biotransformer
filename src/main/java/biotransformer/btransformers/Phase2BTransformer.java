/**
 * This class implements the class of phase II biotransformers, which simulate the transformation of molecules
 * by phase II enzymes.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.btransformers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SMARTSException;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.esaprediction.ESSpecificityPredictor;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.ChemicalClassFinder;
import biotransformer.utils.Utilities;
import exception.BioTransformerException;
import phase2filter.prediction.P2Filter;;

/**
 * @author Yannick Djoumbou Feunang
 *
 */
public class Phase2BTransformer extends Biotransformer{

	/**
	 * @throws ParseException 
	 * @throws IOException 
	 * @throws CDKException 
	 * 
	 */
	P2Filter p2Filter = null;
	
	public Phase2BTransformer(BioSystemName bioSName) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {
		super(bioSName);
		setPhase2EnzymesAndReactionList();
		if(this.getBioSystemName() == BioSystemName.HUMAN){
//			this.p2Filter = new P2Filter();
			this.p2Filter = new P2Filter((LinkedHashMap<String, Object>) this.bSystem.mlmodels.get("P2Filter"));
//			System.out.println("this.p2Filter : " + this.p2Filter.properties);
		}
		
//		printStatistics();
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * Create a list of Phase 2 enzymes for the corresponding biosystem, and group reactions into 
	 * different conjugation types (i.e. glucuronidation, sulfonation, methylation, acetylation, and glutathione
	 * transfer).
	 * @throws IOException 
	 * @throws JsonMappingException 
	 * @throws JsonGenerationException 
	 */
	private void setPhase2EnzymesAndReactionList() throws JsonGenerationException, JsonMappingException, IOException{
		
		this.enzymesByreactionGroups.put("Glucuronidation", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("Sulfonation", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("Acetylation", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("Methylation", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("Glutathione Transfer", new ArrayList<Enzyme>());
		this.enzymesByreactionGroups.put("Glycine Transfer", new ArrayList<Enzyme>());

			
		this.reactionsByGroups.put("Glucuronidation", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("Sulfonation", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("Acetylation", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("Methylation", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("Glutathione Transfer", new ArrayList<MetabolicReaction>());
		this.reactionsByGroups.put("Glycine Transfer", new ArrayList<MetabolicReaction>());
		
//		System.out.println(mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
//				this.reactionsByGroups));
//		System.out.println(this.bSystem.getEnzymesList());		
		
		if(this.bSystem.name == BioSystemName.HUMAN){
			for(Enzyme enz : this.bSystem.getEnzymesList()){
				if(enz.getName().contentEquals("EC_2_4_1_17")){				
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Glucuronidation").add(enz);
					this.reactionsByGroups.get("Glucuronidation").addAll(enz.getReactionSet());
//					methylations.addAll(enz.getReactionSet());
				} 
				else if(enz.getName().contentEquals("EC_2_8_2_1") ||
						enz.getName().contentEquals("EC_2_8_2_2")){
//					 || enz.getName().contentEquals("EC_2_8_2_15")
					this.enzymesList.add(enz);					
					this.enzymesByreactionGroups.get("Sulfonation").add(enz);
//					this.reactionsByGroups.put("Sulfonation", enz.getReactionSet());
					this.reactionsByGroups.get("Sulfonation").addAll(enz.getReactionSet());
					
//					System.err.println("\n\n-------------------------------------" + enz.getName() +"\n"+ enz.getReactionSet());
//					System.out.println( "Sulfonation\n\n" + 
//							this.reactionsByGroups.get("Sulfonation"));
				}
				else if(enz.getName().contentEquals("EC_2_1_1_6") || enz.getName().contentEquals("EC_2_1_1_9") || 
						enz.getName().contentEquals("EC_2_1_1_67") ||enz.getName().contentEquals("EC_2_1_1_96")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Methylation").add(enz);
					this.reactionsByGroups.get("Methylation").addAll(enz.getReactionSet());
//					methylations.addAll(enz.getReactionSet());
//					System.err.println("\n\n-------------------------------------" + enz.getName() +"\n"+ enz.getReactionSet());

				}
				else if(enz.getName().contentEquals("ACETYLTRANSFERASE")  || enz.getName().contentEquals("EC_2_3_1_5")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Acetylation").add(enz);
//					this.reactionsByGroups.put("Acetylation", enz.getReactionSet());
					this.reactionsByGroups.get("Acetylation").addAll(enz.getReactionSet());
//					System.err.println("\n\n-------------------------------------" + enz.getName() +"\n"+ enz.getReactionSet());

				}
				else if(enz.getName().contentEquals("EC_2_5_1_18")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Glutathione Transfer").add(enz);
//					this.reactionsByGroups.put("Glutathione Transfer", enz.getReactionSet());
					this.reactionsByGroups.get("Glutathione Transfer").addAll(enz.getReactionSet());
//					System.err.println("\n\n-------------------------------------" + enz.getName() +"\n"+ enz.getReactionSet());

				}
				
				else if(enz.getName().contentEquals("EC_2_3_1_13")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Glycine Transfer").add(enz);
//					this.reactionsByGroups.put("Glycine Transfer", enz.getReactionSet());
					this.reactionsByGroups.get("Glycine Transfer").addAll(enz.getReactionSet());
					
//					System.err.println(enz.getReactionSet().getClass());				
//					System.err.println("\n\n-------------------------------------" + enz.getName() +"\n"+ enz.getReactionSet());
//					System.out.println(
//							this.reactionsByGroups.get("Glycine Transfer"));
				}

			}
////			System.out.println("standardizationReactions : " + mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
////					this.reactionsByGroups.get("standardizationReactions")));
//			System.out.println("Glucuronidation : " + 
//					this.reactionsByGroups.get("Glucuronidation"));
//
//			System.out.println("\n--------------------------------------\nSulfonation : " + 
//					this.reactionsByGroups.get("Sulfonation"));
//
//			System.out.println("\n--------------------------------------\nMethylation : " + 
//					this.reactionsByGroups.get("Methylation"));
//
//			System.out.println("\n--------------------------------------\nAcetylation : " + 
//					this.reactionsByGroups.get("Acetylation"));
//
//			System.out.println("\n--------------------------------------\nGlutathione Transfer : " + 
//					this.reactionsByGroups.get("Glutathione Transfer"));
//
//			System.out.println("\n--------------------------------------\nGlycine Transfer : " + 
//					this.reactionsByGroups.get("Glycine Transfer"));
//			
			
		}else if (this.bSystem.name == BioSystemName.GUTMICRO){
			for(Enzyme enz : this.bSystem.getEnzymesList()){
				if(enz.getName().contains("UDP_GLUCURONOSYLTRANSFERASE")){				
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Glucuronidation").add(enz);
					this.reactionsByGroups.get("Glucuronidation").addAll(enz.getReactionSet());
				} 
				else if(enz.getName().contains("SULFOTRANSFERASE")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Sulfonation").add(enz);
					this.reactionsByGroups.get("Sulfonation").addAll(enz.getReactionSet());
				}
				else if(enz.getName().trim() =="EC_2_1_1_6" || enz.getName().trim() == "EC_2_1_1_9" || 
						enz.getName().contains("EC_2_1_1_67") ||enz.getName().contains("EC_2_1_1_96")){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("Methylation").add(enz);
					this.reactionsByGroups.get("Methylation").addAll(enz.getReactionSet());
					
				}
			}
		}
//		System.out.println(this.reactionsByGroups);
		
		for(String key : this.reactionsByGroups.keySet()){
//			System.out.println(key);
			if( this.reactionsByGroups.get(key) != null) {
				for(MetabolicReaction mr : this.reactionsByGroups.get(key)){
					if(! this.reactionsHash.containsKey(mr.getReactionName())){
						this.reactionsHash.put(mr.getReactionName(), mr);
					}
				}				
			}

		}
	}

	/**
	 * 
	 * @param target
	 * @param preprocess
	 * @param filter
	 * @return an arraylist of biotransformations, which are instances of the Phase II metabolic reactions applied to the 
	 * target, with the set minimum threshold of 0.0
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> applyPhase2Transformations(IAtomContainer target, boolean precheck, boolean preprocess, boolean filter) throws Exception {	
		return applyPhase2Transformations(target, precheck, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target
	 * @param preprocess
	 * @param filter
	 * @param scoreThreshold
	 * @return an arraylist of biotransformations, which are instances of the Phase II metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> applyPhase2Transformations(IAtomContainer target, boolean precheck, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//		System.out.println("molecule is aromatic: " + ChemStructureExplorer.isAromatic(target));
		ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();
		
		try{
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			}
			
			else{
				if(ChemStructureExplorer.isBioTransformerValid(target) && (precheck == true && isPotentialPhase2SubstrateByReactionPatternMatching(target)) || precheck == false){
//					System.out.println("test 1");
					if(this.getBioSystemName() == BioSystemName.HUMAN || this.getBioSystemName() == BioSystemName.GUTMICRO){
						if(ChemStructureExplorer.isInvalidPhase2Metabolite(target)){
//							System.out.println("test 2");
							return biotransformations;
						}
					}
					
					ArrayList<String> conjugation_types = new ArrayList<String>(this.reactionsByGroups.keySet());
					conjugation_types.remove("standardizationReactions");
					
					for(String conjugationType : conjugation_types){
						
						// System.out.println("Is this a polypenol or phenolic derivative? " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(target));
						if(!ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(target)){
//							System.err.println("Applying " + conjugationType + " reactions.");
	//						biotransformations.addAll(applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get(conjugationType), preprocess, filter, scoreThreshold));
							biotransformations.addAll(this.metabolizeWithEnzymes(target, this.enzymesByreactionGroups.get(conjugationType), preprocess, filter, scoreThreshold));
//							System.out.println(biotransformations.size());
						} else {
							
							if(conjugationType == "Glucuronidation" || conjugationType == "Methylation" || conjugationType == "Sulfonation" || 
									conjugationType == "Glycine Transfer"){
//								System.out.println("test 3");
//								System.err.println("Applying " + conjugationType + " reactions.");
	//							biotransformations.addAll(applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get(conjugationType), preprocess, filter, scoreThreshold));					
								biotransformations.addAll(this.metabolizeWithEnzymes(target, this.enzymesByreactionGroups.get(conjugationType), preprocess, filter, scoreThreshold));
							}
						}
					}
					
					if(this.getBioSystemName() == BioSystemName.HUMAN || this.getBioSystemName() == BioSystemName.GUTMICRO){
						for(Biotransformation bt : biotransformations){
							boolean g = true;
							for(IAtomContainer at : bt.getProducts().atomContainers() ){
//								System.out.println("test 4");
			//					at = ChemStructureManipulator.preprocessContainer(at);
								if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
	//								biotransformations.remove(bt);
//									System.out.println("test 5");
									g =  false;									
									break;
								}
//								else{
//									
//									System.err.println("Invalid Compound: " + this.smiGen.create(at));
//								}
							}
							
							if(g == true){
								selectedBiotransformations.add(bt);
							}
						}
					}
				
				}
			}
		} catch (Exception iae){
		iae.printStackTrace();
		return null;
		}
		

		
		return selectedBiotransformations;
	}
	
	public ArrayList<Biotransformation> applyPhase2Transformations(IAtomContainerSet targets, boolean precheck, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//		IAtomContainerSet filteredTargets = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		if(precheck == true){
//			filteredTargets = this.reduceSet(targets);
//		}else{
//			filteredTargets = targets;
//		}
		
//		for(IAtomContainer at : filteredTargets.atomContainers()){
		for(IAtomContainer at : targets.atomContainers()){
			ArrayList<Biotransformation> bt = applyPhase2Transformations(at, precheck, preprocess, filter, scoreThreshold);
			if(bt.size()>0){
				biotransformations.addAll(bt);
			}
		}
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	public ArrayList<Biotransformation> applyPhase2TransformationsChainAndReturnBiotransformations(IAtomContainer target,
			boolean precheck, boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return applyPhase2TransformationsChainAndReturnBiotransformations(target, precheck, preprocess, filter, 
				nr_of_steps, 0.0);
	}
	
	public ArrayList<Biotransformation> applyPhase2TransformationsChainAndReturnBiotransformations(IAtomContainer target,
			boolean precheck, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{

		AtomContainerSet startingSet = new AtomContainerSet();
		startingSet.addAtomContainer(target);
		
		return applyPhase2TransformationsChainAndReturnBiotransformations(startingSet,
				preprocess, precheck, filter, nr_of_steps, scoreThreshold);
	}

	public ArrayList<Biotransformation> applyPhase2TransformationsChainAndReturnBiotransformations(IAtomContainerSet targets,
			boolean precheck, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> products = new ArrayList<Biotransformation>();	
		
		
//		AtomContainerSet containers = (AtomContainerSet) targets;
		IAtomContainerSet containers = targets;
		int counter = 0;
		
		while(nr_of_steps>0){
			counter++;			
			ArrayList<Biotransformation> currentProducts = applyPhase2Transformations(containers, precheck, preprocess, filter, scoreThreshold);
			nr_of_steps--;
//			System.err.println(currentProducts.size() + " biotransformations at step " + counter);
			if(!currentProducts.isEmpty()){
				products.addAll(currentProducts);
				containers.removeAllAtomContainers();
				containers = extractProductsFromBiotransformations(currentProducts);
				for(IAtomContainer a : containers.atomContainers()){
					
//					if(ChemStructureExplorer.isCompoundInorganic(a)) {
//						
//					}else {
						ChemStructureManipulator.standardizeMoleculeWithCopy(a);
//					}
//					 AtomContainerManipulator.convertImplicitToExplicitHydrogens(a);				
				}
//				System.err.println("Number of compounds for upcoming setp " + (counter + 1) + ": " + containers.getAtomContainerCount());
			}
			else{
				break;
			}
		}
//		System.out.println("Stopped after " + counter + " steps.");
//		System.out.println("Number of phase 2 metabolites: " + products.size());
		return products;
	}

	public ArrayList<Biotransformation> applyReactionAndReturnBiotransformations(IAtomContainer target,
			MetabolicReaction reaction, boolean preprocess, Double scoreThreshold) throws Exception{
//		System.out.println(reaction.name);
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();		
		
		

//		IAtomContainer starget = this.standardizeMoleculeWithCopy(target);
		IAtomContainer starget = target.clone();
		
		if (preprocess) {
			try {
				starget = ChemStructureManipulator.preprocessContainer(starget);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		} else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(starget);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
		}
		
		if(isPotentialPhase2SubstrateByReactionPatternMatching(target)){
			ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
			ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			if(target.getProperty("InChI") == null || ((String) target.getProperty("InChI")).trim().length()==0){
				InChIGenerator gen0 = factory.getInChIGenerator(target);
				target.setProperty("InChI", gen0.getInchi());
				target.setProperty("InChIKey", gen0.getInchiKey());
				Utilities.addPhysicoChemicalProperties(target);
				target.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(target));
			}
			if(target.getProperty("SMILES") == null) {
				target.setProperty("SMILES", this.smiGen.create(AtomContainerManipulator.removeHydrogens(target)));
//				System.err.println(this.smiGen.create(AtomContainerManipulator.removeHydrogens(target)));			
			}
//			System.out.println("STARGET: " + this.smiGen.create(starget));
			IAtomContainerSet partialSet = generateAllMetabolitesFromAtomContainer(starget, reaction, true);
			Double score=0.0;
			AtomContainerSet subs = new AtomContainerSet();
			AtomContainerSet prod = new AtomContainerSet();
//			System.out.println("Partial set size: " +  partialSet.getAtomContainerCount());
					
			if(partialSet.getAtomContainerCount()>0){	
				
				

				if(target.getProperty("Score") !=null && !ChemStructureExplorer.isInvalidPhase2Metabolite(target)){					
					score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(reaction.name)  );
				}else if(!ChemStructureExplorer.isInvalidPhase2Metabolite(target)){
					score = this.bSystem.getReactionsORatios().get(reaction.name);
				}
//				System.out.println("score: " + score);
							
				if(score>=scoreThreshold){
					subs.addAtomContainer(target);
					
					for(IAtomContainer pc : partialSet.atomContainers()){
//						System.out.println(this.smiGen.create(pc));
//						pc = ChemStructureManipulator.preprocessContainer(pc);
						if(!ChemStructureExplorer.isInvalidPhase2Metabolite(pc)){
//							System.err.println("Adding this compound: " + this.smiGen.create(pc));
							
							if(pc.getProperty("InChI") == null){
								InChIGenerator gen = factory.getInChIGenerator(pc);
								pc.setProperty("InChI", gen.getInchi());
								pc.setProperty("InChIKey", gen.getInchiKey());	
								pc.setProperty("SMILES", this.smiGen.create(AtomContainerManipulator.removeHydrogens(pc)));
							}					
							Utilities.addPhysicoChemicalProperties(pc);
							pc.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(pc));
							prod.addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
							Biotransformation bioT = new Biotransformation(subs, reaction.name, null, prod, score, this.bSystem.name );
							results.add(bioT);
							
						}
//						else{
//							System.err.println("Invalid Compound: " + this.smiGen.create(pc));
//						}
					}
				}
			}
		}
		

		return results;
	}

	public ArrayList<Biotransformation> applyReactionsAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();		
//		IAtomContainer starget = this.standardizeMoleculeWithCopy(target);
		IAtomContainer starget = target.clone();
		
		if (preprocess) {
			try {
				starget = ChemStructureManipulator.preprocessContainer(starget);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
//				System.out.println("After preprocessing: " + this.smiGen.create(starget));
			}
			catch (Exception e) {
				System.out.println(e);
			}
		} else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(starget);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
		}
	if(isPotentialPhase2SubstrateByReactionPatternMatching(starget)){	
	//		System.err.println("Is this molecule aromatic? " + ChemStructureExplorer.isAromatic(starget));
			ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
			ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();		
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		
			for (MetabolicReaction i : reactions) {
				boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(i, starget);
				if (match_constraints) {
	//				System.out.println("Compound matches " + i.name);
					matchedReactions.add(i);
				}
			}		
			if(filter == false){
				filteredReactions = matchedReactions;		
			} else{
				filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
			}
	//		System.out.println("Number of reactions after filtering: " + filteredReactions.size());
	
			for(MetabolicReaction j : filteredReactions){
					
				ArrayList<Biotransformation> partial = applyReactionAndReturnBiotransformations(starget, j, false, scoreThreshold);
				
				if(partial.size()>0){
					results.addAll(partial);
				}
			}
		}
		return results;
	}
	
	public void applyReactionChainFromSDF(String inputFileName,
			boolean precheck, boolean preprocess, boolean filter, 
			int nr_of_steps, Double scoreThreshold) throws Exception{
			
			ArrayList<Biotransformation> allPhaseIImetabolites = new ArrayList<Biotransformation>();
			IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(inputFileName),
					builder);
			while (sdfr.hasNext()){
				IAtomContainer mol_ = sdfr.next();
				IAtomContainer mol = ChemStructureManipulator.preprocessContainer(mol_);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
								
				ArrayList<Biotransformation> molPhaseIImetabolites = this.applyPhase2TransformationsChainAndReturnBiotransformations(mol, precheck,
						preprocess, filter, nr_of_steps, scoreThreshold);
				allPhaseIImetabolites.addAll(molPhaseIImetabolites);	
				
				
			}
			
//			System.out.println("allPhaseIImetabolites.size() = " + allPhaseIImetabolites.size());
//			System.out.println(inputFileName);
//			System.out.println(inputFileName.replace(".sdf", ""));
			if(allPhaseIImetabolites.size()>0){
				this.saveBioTransformationProductsToSdf(allPhaseIImetabolites, inputFileName.replace(".sdf", "") + "_PhaseII_metabolites.sdf");			
			}
			sdfr.close();
		}
	
	public void applyReactionChainFromSdfToSingleSdf(String inputFileName,
			boolean precheck, boolean preprocess, boolean filter, 
			int nr_of_steps, Double scoreThreshold, String outputFileName, boolean annotate) throws Exception{

			ArrayList<Biotransformation> allPhaseIImetabolites = new ArrayList<Biotransformation>();
			IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(inputFileName),
					builder);
			int count = 0;
			while (sdfr.hasNext()){
				count++;
				
//				System.out.println("Count: " + count);
				IAtomContainer mol_ = sdfr.next();
//				System.out.println(this.smiGen.create(mol_));
				IAtomContainer mol = ChemStructureManipulator.preprocessContainer(mol_);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
								
				ArrayList<Biotransformation> molPhaseIImetabolites = this.applyPhase2TransformationsChainAndReturnBiotransformations(mol,
						precheck, preprocess, filter, nr_of_steps, scoreThreshold);

				allPhaseIImetabolites.addAll(this.selectBiotransformationWithValidProducts(molPhaseIImetabolites));		
			}
			if(allPhaseIImetabolites.size()>0){
				this.saveBioTransformationProductsToSdf(allPhaseIImetabolites, outputFileName, annotate);			
			}
			sdfr.close();
		}
	
	
	
	public boolean isPotentialHumanPhase2SubstrateByReactionPatternMatching(IAtomContainer substrate) throws Exception{
		boolean phaseII = false;
		
		
//		return (ChemStructureExplorer.getMajorIsotopeMass(substrate) < 900.0 && p2Filter.filter(substrate).get(0) == "T");
		return p2Filter.filter(substrate).get(0) == "T";

	}

	
	public boolean isPotentialNonHumanPhase2SubstrateByReactionPatternMatching(IAtomContainer substrate) throws Exception{
		boolean phaseII = false;


			
			ESSpecificityPredictor ess = new ESSpecificityPredictor(this.bSystem);
			
			// WATCH OUT: CERAMIDES CAN BE GLUCURONIDATED - SEE http://cancerres.aacrjournals.org/content/59/22/5768.long
			if(ChemicalClassFinder.isEtherLipid(substrate) || ChemicalClassFinder.isGlyceroLipid(substrate) || 
					ChemicalClassFinder.isGlycerophosphoLipid(substrate) || ChemicalClassFinder.isSphingoLipid(substrate)
					|| ChemicalClassFinder.isAcylCoAConjugate(substrate)){
				//do nothing - phaseII = false
			}
			else
			{
				for(Enzyme e : this.enzymesList){
					if(ess.isValidSubstrate(substrate, e.getName())){
						phaseII = true;
						break;
					}
				}

		}
			



		return phaseII;
	}

	
	// FIX THIS - ADD FUNCTIONS IN ESSpecificityPredictor THAT TAKES ENZYMES AS ARGUMENTS, TO AVOIR RECREATING ENZYMES VIA ENZYME NAMES.
	public boolean isPotentialPhase2SubstrateByReactionPatternMatching(IAtomContainer substrate) throws Exception{
		boolean phaseII = false;
		
		if(this.getBioSystemName() == BioSystemName.HUMAN){
			phaseII = p2Filter.filter(
					substrate 
					).get(0) == "T";
			
		}
		else{
			
			ESSpecificityPredictor ess = new ESSpecificityPredictor(this.bSystem);
			
			// WATCH OUT: CERAMIDES CAN BE GLUCURONIDATED - SEE http://cancerres.aacrjournals.org/content/59/22/5768.long
			if(ChemicalClassFinder.isEtherLipid(substrate) || ChemicalClassFinder.isGlyceroLipid(substrate) || 
					ChemicalClassFinder.isGlycerophosphoLipid(substrate) || ChemicalClassFinder.isSphingoLipid(substrate)){

				
				//do nothing - phaseII = false
			}
			else
			{
				for(Enzyme e : this.enzymesList){
//					System.out.println(e.getName());
					if(ess.isValidSubstrate(substrate, e.getName())){
						phaseII = true;
						break;
					}
				}
			}
		}
			



		return phaseII;
	}
	

	public IAtomContainerSet reduceSet(IAtomContainerSet molecules) throws CDKException, Exception{

		IAtomContainerSet filteredMolecules =  DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		if(this.getBioSystemName() == BioSystemName.HUMAN){
			filteredMolecules = this.p2Filter.returnFilteredPhaseIICandidates(molecules);
		}
		else{
			
			for(int i = 0; i < molecules.getAtomContainerCount(); i++){
				if(this.isPotentialNonHumanPhase2SubstrateByReactionPatternMatching(molecules.getAtomContainer(i))){
					filteredMolecules.addAtomContainer(molecules.getAtomContainer(i));
				}				
			}

		}
		
		return filteredMolecules;
	}
	
	public ArrayList<Biotransformation> selectBiotransformationWithValidProducts(ArrayList<Biotransformation> biotransformations) throws SMARTSException, CDKException, CloneNotSupportedException{
		
		ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();
		
		for(Biotransformation bt : biotransformations) {
			boolean goodBiontransfo =  true;

//			System.out.println("Biotransformation type. " +  bt.getReactionType());
			for( IAtomContainer at : bt.getProducts().atomContainers() ){
				
//				System.out.println(at.getProperty("InChI"));
//				at = ChemStructureManipulator.preprocessContainer(at);
				if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
//					System.err.println("is invalid phase 2 metabolite");
					goodBiontransfo = false;					
					break;
				}
//				else{
//					System.err.println("Invalid Compound: " + this.smiGen.create(at));
//				}
			}
			
//			System.out.println("goodBiontransfo: " + goodBiontransfo);
			if(goodBiontransfo){
				selectedBiotransformations.add(bt);
			}
		}
		
		return selectedBiotransformations;
	}
	
	
	
//	public boolean isPhaseIIMetabolite(IAtomContainer substrate){
//		boolean phaseII = false;
//		
//		
//		
//		
//		return phaseII;
//	}
	
	
	public void printStatistics(){
//		System.out.println(this.reactionsByGroups.get("cypReactions").size());
		int count = 0;
		for(Enzyme e: this.enzymesList){
//			System.out.println(e.getName() + " : " + e.getReactionsNames().size());
			count = count + e.getReactionsNames().size();
		}
		System.out.println("Humber of enzymes: " + this.enzymesList.size());
		System.out.println("Humber of biotransformation rules: " 
				+ (this.reactionsByGroups.get("Glucuronidation").size()
				+ this.reactionsByGroups.get("Sulfonation").size()
				+ this.reactionsByGroups.get("Acetylation").size()
				+ this.reactionsByGroups.get("Methylation").size()
				+ this.reactionsByGroups.get("Glutathione Transfer").size()
				+ this.reactionsByGroups.get("Glycine Transfer").size()
				
				));
		System.out.println("Number of enzyme-biotransformation rules associations: " + count);
	}
}
