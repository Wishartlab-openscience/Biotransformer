/**
 * This class implements the class of environmental microbial biotransformers, which simulate the transformation of molecules
 * by enzymes from microbial species found in the soil and aquatic environments. It implements rules and constraints developed
 * by the EAWAG/IUMBBD (http://eawag-bbd.ethz.ch/).
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.btransformers;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import exception.BioTransformerException;

public class EnvMicroBTransformer extends Biotransformer {

	/**
	 * @throws ParseException 
	 * @throws IOException 
	 * @throws CDKException 
	 * 
	 */
	
	public EnvMicroBTransformer() throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {
		super(BioSystemName.ENVMICRO);
		setEnzymesList();
		setReactionsList();

	}
	
	/**
	 * Collect the list of enzymes associated with the current biosystem from the database.
	 */
	
	private void setEnzymesList(){
		this.enzymesByreactionGroups.put("envMicroReactions", new ArrayList<Enzyme>());
		
		for(Enzyme enz : this.bSystem.getEnzymesList()){
//			if(enz.getName().contains("EC_") || enz.getName().contains("UDP_GLUCURONOSYLTRANSFERASE") || 
//				enz.getName().contains("SULFOTRANSFERASE") || enz.getName().contains("ACETYLTRANSFERASE") ||
//				enz.getName().contentEquals("ABKAR1") || enz.getName().contentEquals("ABKAR2") || enz.getName().contentEquals("CYPB5_CYPBR5")){
				

				if(!enz.getReactionSet().isEmpty()){
					this.enzymesList.add(enz);
					this.enzymesByreactionGroups.get("envMicroReactions").add(enz);
				}				
//			}
		}
	}
	
	/**
	 * Collect the list of metabolic reactions associated with the current biosystem, inferred from the list of enzymes.
	 */
	private void setReactionsList(){
//		this.reactionsList.put("standardizationReactions", MReactionSets.standardizationReactions);
		this.reactionsByGroups.put("envMicroReactions", new ArrayList<MetabolicReaction>(this.bSystem.getReactionsHash().values()));
		this.reactionsHash.putAll(this.bSystem.getReactionsHash());
	}
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the environmental microbial reactions applied to the target, with a threshold of 0.0
	 * @throws Exception
	 */	
	public ArrayList<Biotransformation>  applyEnvMicrobialTransformations(IAtomContainer target, boolean preprocess, boolean filter)
			throws Exception {
		return applyEnvMicrobialTransformations(target, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param scoreThreshold - minimum threshold for reaction scores
	 * @return an arraylist of biotransformations, which are instances of the environmental microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */	
	public ArrayList<Biotransformation> applyEnvMicrobialTransformations(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		
		try{
		
				if(ChemStructureExplorer.isPpsValid(target)){
				
				ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		//		biotransformations = applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get("envMicroReactions"), preprocess, filter, scoreThreshold);		
				
				biotransformations = metabolizeWithEnzymes(target,
						this.enzymesByreactionGroups.get("envMicroReactions"), preprocess, filter, scoreThreshold);
				
				return biotransformations;
			}
			else {
				throw new IllegalArgumentException("\n\n" + this.smiGen.create(target) + "\n\n"
						+ "For the prediction of environmental microbial metabolism, the compound must: "
						+ "1) be organic; 2) not be a mixture; 3) not be a cofactor or dead end compound, and 4) have a molecular mass of 1000 Da. or less.\nFor more"
						+ " information, consult the following link: http://eawag-bbd.ethz.ch/predict/notbepredicted.html#Unknowncomp\n\n\n");
			} 
		}
		catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * 
	 * @param targets
	 * @param preprocess
	 * @param filter
	 * @param scoreThreshold
	 * @return an arraylist of biotransformations, which are instances of the environmental microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */	
	public ArrayList<Biotransformation> applyEnvMicrobialTransformations(IAtomContainerSet targets, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		for(IAtomContainer atc : targets.atomContainers()){
			biotransformations.addAll(metabolizeWithEnzymes(atc,
					this.enzymesByreactionGroups.get("envMicroReactions"), preprocess, filter, scoreThreshold));
		}
//		biotransformations = applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get("envMicroReactions"), preprocess, filter, scoreThreshold);		
		
		return biotransformations;
	}	
	
	
	

//	biotransformations = metabolizeWithEnzymesBreadthFirst(target,
//			this.enzymesByreactionGroups.get("ecBasedDeconjugations"), preprocess, filter, scoreThreshold);

	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the environmental microbial metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyEnvMicrobialTransformationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps)
			throws Exception {
		return applyEnvMicrobialTransformationsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the environmental microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	

	public ArrayList<Biotransformation> applyEnvMicrobialTransformationsChain(IAtomContainer target, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
		
		try{
		
			if(ChemStructureExplorer.isPpsValid(target)){
				ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
				AtomContainerSet startingSet = new AtomContainerSet();
				startingSet.addAtomContainer(target);
				
//				System.out.println("this.enzymesByreactionGroups.get(\"envMicroReactions\")\n");
//				for ( Enzyme r : this.enzymesByreactionGroups.get("envMicroReactions")) {
//					System.out.println(r.toString());
//				}
		//		biotransformations = applyReactionsChainAndReturnBiotransformations(startingSet, this.reactionsByGroups.get("envMicroReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
				biotransformations = metabolizeWithEnzymesBreadthFirst(startingSet,
						this.enzymesByreactionGroups.get("envMicroReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
		
				
				return biotransformations;
			}
						
			else {
				throw new IllegalArgumentException(this.smiGen.create(target) + "\n"
						+ "For the prediction of environmental microbial metabolism, the compound must: "
						+ "1) be organic; 2) not be a mixture; 3) not be a cofactor or dead end compound, and 4) have a mmolecular mass of 1000 Da. or less.\nFor more"
						+ " information, consult the following link: http://eawag-bbd.ethz.ch/predict/notbepredicted.html#Unknowncomp");
			} 
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	public ArrayList<Biotransformation> applyEnvMicrobialTransformationsChain(IAtomContainerSet targets, boolean preprocess, boolean filter, int nr_of_steps, double scoreThreshold)
			throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		for(IAtomContainer target : targets.atomContainers()) {
			biotransformations.addAll(applyEnvMicrobialTransformationsChain(targets, preprocess, filter, nr_of_steps, scoreThreshold));
		}
		
		return biotransformations;
	}

	
	public void simulateEnvMicrobialDegradationAndSaveToSDF(IAtomContainer container, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.applyEnvMicrobialTransformationsChain(container, preprocess, filter, nrOfSteps, scoreThreshold);
		this.saveBioTransformationProductsToSdf(biotransformations, outputFileName, annotate);
	}
	
	public void simulateEnvMicrobialDegradationAndSaveToCSV(IAtomContainer container, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.applyEnvMicrobialTransformationsChain(container, preprocess, filter, nrOfSteps, scoreThreshold);
		this.saveBioTransformationProductsToCSV(biotransformations, outputFileName, annotate);
	}
	
	
	public ArrayList<Biotransformation> simulateEnvMicrobialDegradation(IAtomContainerSet containers, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations  =  new ArrayList<Biotransformation>();
		
		for(IAtomContainer atc: containers.atomContainers()){
			biotransformations.addAll(this.applyEnvMicrobialTransformationsChain(atc, preprocess, filter, nrOfSteps, scoreThreshold));
		}
		
		return biotransformations;		
	}
	
	public void simulateEnvMicrobialDegradationAndSaveToSDF(IAtomContainerSet containers, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold, String outputFolder, boolean annotate) throws Exception{
		
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
				identifier = identifier.replace(":", "-").replace("/", "_");
				System.out.println(identifier);
				this.simulateEnvMicrobialDegradationAndSaveToSDF(molecule, preprocess, filter, nrOfSteps, scoreThreshold, outputFolder + "/" + identifier + "_env_based_metabolites.sdf", annotate);			
			}
		}	
		
	}
	
	public ArrayList<Biotransformation> applyReactionsAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();	
		IAtomContainer starget = ChemStructureManipulator.standardizeMoleculeWithCopy(target, preprocess);
		
//		if (preprocess) {
//			try {
//				starget = ChemStructureManipulator.preprocessContainer(starget);
////				System.out.println("After preprocessing: " + this.smiGen.create(starget));
//			}
//			catch (Exception e) {
//				System.out.println(e);
//			}
//		}
		try{
			if(ChemStructureExplorer.isPpsValid(target)){
				System.out.println("The target is ppt valid");
				ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
				ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();		
				InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
				InChIGenerator gen0 = factory.getInChIGenerator(target);
				target.setProperty("InChI", gen0.getInchi());
				target.setProperty("InChIKey", gen0.getInchiKey());
			
				for (MetabolicReaction i : reactions) {
		//			System.out.println(i.name);
		//			System.out.println(i.getReactionSMIRKS());
					boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(i, starget);
					//			System.out.println(i.name);
					if (match_constraints) {
	//									System.out.println("Compound matches " + i.name + ": " + match_constraints);
						matchedReactions.add(i);
					}
				}		
				if(filter == false){
					filteredReactions = matchedReactions;		
				} else{
					filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
				}
				
				for(MetabolicReaction j : filteredReactions){
	//				System.out.println(j.name);
					IAtomContainerSet partialSet = generateAllMetabolitesFromAtomContainer(starget, j, false);
								
	//							System.out.println(partialSet.getAtomContainerCount());
					Double score=0.0;
					AtomContainerSet subs = new AtomContainerSet();
					AtomContainerSet prod = new AtomContainerSet();
					
					
					if(partialSet.getAtomContainerCount()>0){
						
						if(target.getProperty("Score") != null){	
							
							score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(j.name)  );
						}else{
							score = this.bSystem.getReactionsORatios().get(j.name);
						}
						//				System.out.println(score);
		
						if(score>=scoreThreshold){
							subs.addAtomContainer(target);
							for(IAtomContainer pc : partialSet.atomContainers()){
								InChIGenerator gen = factory.getInChIGenerator(pc);
								pc.setProperty("InChI", gen.getInchi());
								pc.setProperty("InChIKey", gen.getInchiKey());
								prod.addAtomContainer(pc);
							}
							
							Biotransformation bioT = new Biotransformation(subs, j.name, null, prod, score, this.bSystem.name );
							results.add(bioT);
						}
					}	
				}
			}else{
				throw new IllegalArgumentException(this.smiGen.create(target)+ "\n"
						+ "For the prediction of environmental microbial metabolism, the compound must:"
						+ "1) be organic; 2) not be a mixture; 3) not be a cofactor or dead end compound, and 4) have a mmolecular mass of 1000 Da. or less.\nFor more"
						+ " information, consult the following link: http://eawag-bbd.ethz.ch/predict/notbepredicted.html#Unknowncomp");				
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return results;
	}
	
	public void printStatistics(){
		int count = 0;
		for(Enzyme e: this.enzymesList){
			count = count + e.getReactionsNames().size();
		}
		System.out.println("Humber of enzymes: " + this.enzymesList.size());
		System.out.println("Humber of biotransformation rules: " + this.reactionsByGroups.get("envMicroReactions").size() );		
		System.out.println("Humber of enzyme-biotransformation rules associations: " + count);	
		System.out.println("No. of preference rules");
	}
}
