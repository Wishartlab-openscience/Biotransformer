package biotransformer.utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import exception.BioTransformerException;
import phase2filter.prediction.P2Filter;

public class UniversalBioTransformer {
	protected ECBasedBTransformer ecb;
	protected Cyp450BTransformer cyb;
	protected HGutBTransformer hgb;
	protected Phase2BTransformer p2b;
	protected EnvMicroBTransformer emb;
	protected LinkedHashMap<String, LinkedHashMap<String, String>> compoundDictionary		
													= new LinkedHashMap<String, LinkedHashMap<String, String>>();
	
	protected P2Filter p2filter 					= new P2Filter();
	protected LinkedHashMap<String, MetabolicReaction> combinedReactionsHash							
													= new LinkedHashMap<String, MetabolicReaction>();
	
	public SmilesParser smiParser;
	public SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
	
	
	public UniversalBioTransformer() throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {

		ecb			= new ECBasedBTransformer(BioSystemName.HUMAN);
		cyb 		= new Cyp450BTransformer(BioSystemName.HUMAN);
		hgb 		= new HGutBTransformer();
		p2b 		= new Phase2BTransformer(BioSystemName.HUMAN);
		emb 		= new EnvMicroBTransformer();
		smiParser	= ecb.getSmiParser();
		
		for(Map.Entry<String, MetabolicReaction> m : this.ecb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(m.getKey())){
				this.combinedReactionsHash.put(m.getKey(), m.getValue());
			}
		}
		
		for(Map.Entry<String, MetabolicReaction> n : this.cyb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(n.getKey())){
				this.combinedReactionsHash.put(n.getKey(), n.getValue());
			}
		}
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
//		System.err.println(this.combinedReactionsHash.containsKey("STEROL_16_HYDROXYLATION_PATTERN2"));
		
		for(Map.Entry<String, MetabolicReaction> p : this.hgb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(p.getKey())){
				this.combinedReactionsHash.put(p.getKey(), p.getValue());
			}
		}
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
		
		for(Map.Entry<String, MetabolicReaction> o : this.p2b.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(o.getKey())){
				this.combinedReactionsHash.put(o.getKey(), o.getValue());
			}
		}
		
		for(Map.Entry<String, MetabolicReaction> o : this.emb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(o.getKey())){
				this.combinedReactionsHash.put(o.getKey(), o.getValue());
			}
		}
	}


	public IAtomContainer standardizeMoleculeWithCopy(IAtomContainer target) throws Exception{
		return  ChemStructureManipulator.standardizeMoleculeWithCopy(target);
	}
	
	public IAtomContainerSet extractProductsFromBiotransformations(ArrayList<Biotransformation> biotransformations) throws Exception{
		return this.ecb.extractProductsFromBiotransformations(biotransformations);
	}
	
	public IAtomContainerSet extractProductsFromBiotransformationsWithTransformationData(ArrayList<Biotransformation> biotransformations, boolean annotate) throws Exception{
		return this.ecb.extractProductsFromBiotransformationsWithTransformationData(biotransformations, this.combinedReactionsHash, annotate);
	}

	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
	
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
}
