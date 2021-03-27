/**
 * This class implements the class of biosystems. They can represent either individual
 * species/organisms or a collection thereof
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */
package biotransformer.biosystems;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.apache.commons.io.FilenameUtils;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.JsonParser.Feature;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import ambit2.smarts.SMIRKSManager;
import biotransformer.biomolecule.Enzyme;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.JsonUtils;
import exception.BioTransformerException;

/**
 * @author Yannick Djoumbou Feunang
 *
 */
public class BioSystem {

	/**
	 * 
	 */
	// Contains a mapping of ezymes (names) to a set of metabolic reactions they catalyze.
	protected LinkedHashMap<String, MetabolicReaction> reactionsHash;
	
	protected LinkedHashMap<String, Double>	reactionsOcurrenceRatios;
	protected ArrayList<Enzyme>	enzymes;
	protected LinkedHashMap<String, Enzyme>	enzymesHash;
	
	// Update this to include references too
	protected LinkedHashMap<String, ArrayList<Enzyme>> metPathwaysHash;
	public LinkedHashMap<String, Object> reactionPrecedenceRules;
//	public MReactionsFilter mrFilter;
	protected SMIRKSManager smrkMan;		
	public BioSystemName name;
	
	// The following attribute contains location for different models.
	public LinkedHashMap<String, Object> mlmodels;

	
	
	public static void main (String[] args) throws JsonParseException, JsonMappingException, IOException, BioTransformerException{
		ObjectMapper mapper = new ObjectMapper();
		mapper.configure(Feature.ALLOW_COMMENTS, true);
		mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);
		
		BioSystem bt = new BioSystem(BioSystemName.GUTMICRO, mapper);
		 LinkedHashMap<String,Object> kbase = bt.loadBioSystemSpecificKnowledgebaseAndMetadata(mapper); 
//		System.out.println(mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
//		kbase));
	}
	
	public BioSystem(BioSystemName bsName, ObjectMapper mapper) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException {
		this.reactionsHash = new LinkedHashMap<String, MetabolicReaction> ();
		this.reactionsOcurrenceRatios = new LinkedHashMap<String, Double>();
		this.enzymes = new ArrayList<Enzyme>();
		this.enzymesHash = new  LinkedHashMap<String, Enzyme>();
		this.metPathwaysHash = new LinkedHashMap<String,ArrayList<Enzyme>>();
		this.smrkMan	= new SMIRKSManager(SilentChemObjectBuilder.getInstance());
		this.mlmodels = null;
		
		this.name 			= bsName;
		annotateBiosystem(mapper);


	}
	
	
	public BioSystem (BioSystemName bsName, ArrayList<Enzyme> enzymes) {
		this.name 			= bsName;
		this.enzymes		= enzymes;
	}
	
	
	public enum BioSystemName {
		HUMAN, ENVMICRO, GUTMICRO
	}

	
	private LinkedHashMap<Object, Object> ingestBioSystemSpecificTable(ObjectMapper mapper, String filePathname, 
			String kb_base_path, boolean use_base_files) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException {
		LinkedHashMap<Object,Object> table = new LinkedHashMap<Object,Object>();
		
		table = JsonUtils.ingestJson(filePathname, mapper, false);
		if(table.size() == 0 && use_base_files) {
			Path newpath = Paths.get(kb_base_path, FilenameUtils.getName(filePathname));
			table = JsonUtils.ingestJson(newpath.toString(), mapper);
		}
		return table;	
	}
	
	
	@SuppressWarnings("unchecked")
	private LinkedHashMap<String,Object> loadBioSystemSpecificKnowledgebaseAndMetadata(ObjectMapper mapper) 
			throws JsonParseException, JsonMappingException, 
			FileNotFoundException, IOException, BioTransformerException {
		String CONFIGPATH 						= "config.json";
		LinkedHashMap<String,Object> kbase 		= new LinkedHashMap<String,Object>();	
		LinkedHashMap<Object,Object> configData = JsonUtils.ingestJson(CONFIGPATH, mapper);	
		String kb_base_path 					= (String)((LinkedHashMap<Object, Object>) configData.
				get("database_locations")).get("kb_base_path");
		
		//		System.out.println(configData);
		
		/*
		 * Assigning the path to machine-learning models
		 */
		this.mlmodels = (LinkedHashMap<String, Object>)((LinkedHashMap<String, Object>) configData.
				get("abstracted_models")).get(this.name.toString());
		//		System.out.println(configData);
		
		/*
		 * Collecting data related to this biosystem.
		 */
		
		LinkedHashMap<Object, Object> kb_biosystem_config 	= 
				(LinkedHashMap<Object, Object>) ((LinkedHashMap<Object, Object>) 
						((LinkedHashMap<Object, Object>) configData.get("database_locations")).
				get("kb_biosystem_path")).get(this.name.toString());
				
		collectEnzymes(kbase, kb_biosystem_config, kb_base_path, mapper);
		collectMetabolicReactions(kbase, kb_biosystem_config, kb_base_path, mapper);
		collectEnzymeReactionMappings(kbase, kb_biosystem_config, kb_base_path, mapper);
		collectBiosystemsEnzymes(kbase, kb_biosystem_config, kb_base_path, mapper);
		collectReactionORatios(kbase, kb_biosystem_config, kb_base_path, mapper);	
		collectReactionPrecedenceRules(kbase, kb_biosystem_config, kb_base_path, mapper);
		
		// System.out.println(mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
		//							kbase));
		return kbase;
		
	}
	
	
	private void collectEnzymes(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {

		/*
		 * Ingesting the descriptions of enzymes.
		 */
		LinkedHashMap<Object, Object>  allEnzymes = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "enzymes.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));		
		//		System.out.println("allEnzymes : " + ((LinkedHashMap<Object, Object>) allEnzymes.get("enzymes")).size() + " enzymes." );
		if(allEnzymes.size() != 0) {
			kbase.put("allEnzymes", allEnzymes.get("enzymes"));
			if(((LinkedHashMap<Object, Object>) kbase.get("allEnzymes")).size() == 0) {
				BioTransformerException b = new BioTransformerException("No enzymes could be found for the " + this.name.toString() + " biosystem.\n"
						+ "\nAborting...\n\n");
				//				System.err.println(b);
				throw b;
			}		
		}else {
			FileNotFoundException me = new FileNotFoundException("Missing file \"enzymes.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");
			
			//			System.err.print(me);
			throw me;
		}
	}
	
	
	private void collectMetabolicReactions(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {
		/*
		 * Ingesting the descriptions of metabolic reactions.
		 */
		LinkedHashMap<Object, Object>  allReactions = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "metabolicReactions.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));		
		//		System.out.println("allEnzymes : " + ((LinkedHashMap<Object, Object>) allReactions.get("reactions")).size() + " reactions." );
		if(allReactions.size() != 0) {
			kbase.put("allReactions", allReactions.get("reactions"));
			
			if(((LinkedHashMap<Object, Object>) kbase.get("allReactions")).size() == 0) {
				BioTransformerException c = new BioTransformerException("No metabolic reactions could be found for the " + this.name.toString() + " biosystem.\n"
						+ "\nAborting...\n\n");
				//				System.err.println(c);
				throw c;
			}
		}else {
			FileNotFoundException mar = new FileNotFoundException("Missing file \"metabolicReactions.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");
			
			//			System.err.print(mar);
			throw mar;	
		}		

		if((boolean) kb_biosystem_config.get("has_pathways") == true) {
			LinkedHashMap<Object, Object>  pathways = ingestBioSystemSpecificTable(mapper, 
					Paths.get(kb_biosystem_config.get("path").toString(), "pathways.json").toString(), 
					kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));
			//			System.out.println(pathways.size());
			if(pathways.size() != 0) {
				kbase.put("pathways", pathways.get("metabolicPathways"));
				//				System.out.println("pathways : " + ((LinkedHashMap<Object, Object>) kbase.get("pathways")).size() + " pathways." );				
				if(((LinkedHashMap<Object, Object>) kbase.get("pathways")).size() == 0) {
					BioTransformerException d = new BioTransformerException("No metabolic pathways could be found for the " + this.name.toString() + " biosystem.\n"
							+ "\nAborting...\n\n");
					//					System.err.println(d);
					throw d;
				}				
			}
			else {
				
			}			
		}
	}
	
	
	private void collectEnzymeReactionMappings(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {

		/*
		 * Ingesting the enzyme-reactions mappings.
		 */
		LinkedHashMap<Object, Object>  enzymeReactionList = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "enzymeReactions.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));		
		//		System.out.println(enzymeReactionList.size() );
		if(enzymeReactionList.size() != 0) {
			
			kbase.put("enzymeReactionList", ((LinkedHashMap<String, ArrayList<String>>) enzymeReactionList.get("eReactionLists")) );
			//			System.out.println(((LinkedHashMap<String, ArrayList<String>>) kbase.get("enzymeReactionList")));
			//			System.out.println("enzymeReactionList : " + ((LinkedHashMap<String, Object>) kbase.get("enzymeReactionList")).size() + " enzyme-reaction mappings." );
			if(((HashMap<Object, Object>) kbase.get("enzymeReactionList")).size() == 0) {
				BioTransformerException e = new BioTransformerException("No metabolic enzyme-reactions mapping could be found for the " + this.name.toString() + " biosystem.\n"
						+ "\nAborting...\n\n");
				//				System.err.println(e);
				throw e;				
			}
		}else {
			FileNotFoundException mer = new FileNotFoundException("Missing file \"enzymeReactions.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");
			
			//			System.err.print(mer);
			throw mer;	
		}		
	}
	
	
	
	private void collectBiosystemsEnzymes(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {
		
		/*
		 * Ingesting the list of enzymes for the biosystem.
		 */
		LinkedHashMap<Object, Object>  bioSysEnzymeList = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "biosystemEnzymes.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));
		
		if(bioSysEnzymeList.size() != 0) {
			kbase.put("bioSysEnzymeList", ((LinkedHashMap<String, ArrayList<String>>) bioSysEnzymeList.get("enzymeLists")).get(this.name.toString()));
		//			System.out.println(mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
		//					(kbase.get("bioSysEnzymeList"))));
		//			System.out.println("bioSysEnzymeList : " + ((ArrayList<String>) kbase.get("bioSysEnzymeList")).size() + " enzyme mappings for this biosystem." );
			
			if(((ArrayList<String>) kbase.get("bioSysEnzymeList")).size() == 0) {				
				BioTransformerException g = new BioTransformerException("No metabolic enzyme mapping could be found for the " + this.name.toString() + " biosystem.\n"
						+ "\nAborting...\n\n");
				//				System.err.println(g);
				throw g;

			}
		}else {
			FileNotFoundException mbe = new FileNotFoundException("Missing file \"biosystemEnzymes.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");
			
			//			System.err.print(mbe);
			throw mbe;
		
		}			
	}
	
	
	private void collectReactionORatios(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {
				
		/*
		 * Ingesting the reaction occurrence rations (scores) for the biosystem.
		 */
		LinkedHashMap<Object, Object>  reactionsOcurrenceRatios = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "biosystemsReactionORatios.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));

		if(reactionsOcurrenceRatios.size() != 0) {		
			LinkedHashMap<String, Double> t = (LinkedHashMap<String, Double>) ((LinkedHashMap<Object, Object>)reactionsOcurrenceRatios.get("reactionsORatios")).get(this.name.toString());
			kbase.put("reactionsOcurrenceRatios", t);
			//			System.out.println(kbase.get("reactionsOcurrenceRatios"));
			//			System.out.println("reactionsOcurrenceRatios : " + ((LinkedHashMap<String, Double>) kbase.get("reactionsOcurrenceRatios")).size() + " reaction-score mappings for this biosystem." );
			
			if(((LinkedHashMap<String, Double>) kbase.get("reactionsOcurrenceRatios")).size() == 0) {
				BioTransformerException e = new BioTransformerException("No metabolic reaction-score mapping could be found for the " + this.name.toString() + " biosystem.\n"
						+ "\nAborting...\n\n");
				//				System.err.println(e);
				throw e;				
			}
		}else {
			FileNotFoundException mrr = new FileNotFoundException("Missing file \"biosystemsReactionORatios.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");
			
			//  System.err.print(mrr);
			throw mrr;	
		}		
	}
	
	private void collectReactionPrecedenceRules(
			LinkedHashMap<String,Object> kbase,
			LinkedHashMap<Object, Object> kb_biosystem_config,
			String kb_base_path,
			ObjectMapper mapper
			) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException {
		/*
		 * Ingesting the precedence rules for the biosystem.
		 */
		LinkedHashMap<Object, Object>  reactionPrecedenceRules = ingestBioSystemSpecificTable(mapper, 
				Paths.get(kb_biosystem_config.get("path").toString(), "reactionPrecedenceRules.json").toString(), 
				kb_base_path, (boolean) kb_biosystem_config.get("use_base_files"));

		if(reactionPrecedenceRules.size() != 0) {		
			LinkedHashMap<String, Object> t = (LinkedHashMap<String, Object>) ((LinkedHashMap<Object, Object>)reactionPrecedenceRules.get("precedenceRules")).get(this.name.toString());
			LinkedHashMap<String, Object> rels = (LinkedHashMap<String, Object>) t.get("relative");
//			String[] strict = (String[]) t.get("strict");
			ArrayList<String> strict = (ArrayList<String>) t.get("strict");
			//			System.out.println("rels : " + rels);
			//			System.out.println("strict : " + strict);
			LinkedHashMap<String, Object> n = new LinkedHashMap<String, Object>();
			n.put("relative", rels);
			n.put("strict", strict);
			kbase.put("reactionPrecedenceRules", n);
			//			System.out.println("reactionPrecedenceRules : " + rels.size() + " relative precedence rules, and "
			//					+ strict.size()
			//					+ " strict precedence rules for this biosystem." );
			
			if(((LinkedHashMap<String, Object>) kbase.get("reactionPrecedenceRules")).size() == 0) {
				System.err.println("No	precedence mapping could be found for the " + this.name.toString() + " biosystem.\n\n");
				
			}
		}else {
			FileNotFoundException mrr = new FileNotFoundException("Missing file \"reactionPrecedenceRules.json\". The file could not be found."
					+ " Make sure it is located in the path ("
					+ kb_biosystem_config.get("path") + ") "
					+ "specified in \"config.json\","			
					+ " or in the base folder (\"" + kb_base_path + "\")."
					+ "\nAborting...\n\n");

			throw mrr;	
		}		
	}
	
	
	
	@SuppressWarnings("unchecked")
	private void annotateBiosystem(ObjectMapper mapper) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException {
		
		/*
		 * Parsing the JSON files that constitute the database.
		 */				

		LinkedHashMap<String,Object> kbase = loadBioSystemSpecificKnowledgebaseAndMetadata(mapper);
		// System.out.println(kbase);

		/*
		 * Parsing the list of enzymes and their attributes.
		 */	
		
		/*
		 * create a unique list of reactions and create them.
		 * This helps because a reaction can be catalyzed by many 
		 * enzymes, and we do not one to create the same object many times.
		 */		
		/*
		 * Now for each enzyme associated with the biosystem, built the reactions arraylist and then create the enzyme.
		 */		
		
		LinkedHashMap<String,Object> allReactions 	= (LinkedHashMap<String,Object>) kbase.get("allReactions");
		LinkedHashMap<String,Object> allEnzymes 	= (LinkedHashMap<String, Object>) kbase.get("allEnzymes");
		
		for( String e : ((ArrayList<String>)kbase.get("bioSysEnzymeList")) ){			
			ArrayList<MetabolicReaction> enzymeSpecificReactionObjects = new ArrayList<MetabolicReaction>();
			// System.out.println(e);
			// System.out.println( ((LinkedHashMap<String, ArrayList<String>>)  kbase.get("enzymeReactionList")).get(e)   );
			for(String s : ((LinkedHashMap<String, ArrayList<String>>)  kbase.get("enzymeReactionList")).get(e) ){
				// System.out.println("s: " + s);
				if(this.reactionsHash.containsKey(s)){
					enzymeSpecificReactionObjects.add(this.reactionsHash.get(s));
				} 				
				else {						
					// System.out.println(mapper.writerWithDefaultPrettyPrinter().writeValueAsString(
					//	(LinkedHashMap<String,Object>) allReactions.get(s)));
					LinkedHashMap<String,Object> mrObj = (LinkedHashMap<String,Object>) allReactions.get(s);
					//	System.out.println(mrObj);
					String commonName = (String)mrObj.get("commonName");
					if(commonName == null || commonName.contentEquals("")){
						commonName = s;
					}
					String reactionBTMRID = (String)mrObj.get("btmrID");
					// System.err.println(reactionBTMRID);
					MetabolicReaction r = new MetabolicReaction(s, commonName, reactionBTMRID, (String)mrObj.get("smirks"), 
							(ArrayList<String>)mrObj.get("smarts"), (ArrayList<String>)mrObj.get("excludedReactantSmarts"), this.smrkMan);
					
					this.reactionsHash.put(s,r);
					enzymeSpecificReactionObjects.add(r);
					//	r.display();
				}
			}
			
			 LinkedHashMap<String, Object> enz = (LinkedHashMap<String, Object>) allEnzymes.get(e);
			 //	System.out.println(enz);
			 LinkedHashMap<String, Object> biosystems = (LinkedHashMap<String, Object>) enz.get("biosystems");
			 //	System.out.println(biosystems.get(this.name.toString()));
			 String description = (String) enz.get("description");

			 ArrayList<String> uniprot_ids = null;
			 ArrayList<String> cellularLocations = null;
			 if(biosystems.containsKey(this.name.toString())){
				 uniprot_ids = (ArrayList<String>)  ((LinkedHashMap<String, Object>) 
						 biosystems.get(this.name.toString())).get("uniprot_ids");
				 cellularLocations = (ArrayList<String>)  ((LinkedHashMap<String, Object>) 
						 biosystems.get(this.name.toString())).get("cellular_locations");
				 
			// System.out.println("uniprot_ids 		= " +  uniprot_ids);
			// System.out.println("cellularLocations 	= " + cellularLocations);				 
			 }

			 String acceptedName =  (String) enz.get("acceptedName");
			 if(acceptedName == null || acceptedName.contentEquals("")){
				 acceptedName =  (String)e;
			 }
			 
			Enzyme enzy = new Enzyme(e, description, uniprot_ids, cellularLocations, enzymeSpecificReactionObjects, acceptedName); 
			this.enzymes.add(enzy);
			this.enzymesHash.put(e, enzy);
			// enzy.display();
		}
		
		//		System.out.println(reactionsHash);
		//		System.err.println(reactionsHash.size());
		//		System.err.println(this.reactionsHash);
		
		/*
		 * Now generating a pathway map.
		 */

		LinkedHashMap<String, Object> allPathways = (LinkedHashMap<String, Object>) kbase.get("pathways");
		for(String pName : allPathways.keySet()) {
//						System.out.println(pName);
			LinkedHashMap<String,Object> lm = (LinkedHashMap<String, Object>) allPathways.get(pName);
//			System.out.println(lm.keySet().contains(this.name.toString()));
			
			if(lm.containsKey(this.name.toString())){
				LinkedHashMap<String,Object> bPaths = (LinkedHashMap<String,Object>) lm.get(this.name.toString());
				//				System.out.println(this.name.toString());
				//				System.out.println(lm);
				ArrayList<String> enz = (ArrayList<String>) bPaths.get("enzymes");
//				System.out.println(enz);
				
				this.metPathwaysHash.put(pName, new ArrayList<Enzyme>());
				
				for(String e_n : enz) {
					this.metPathwaysHash.get(pName).add(enzymesHash.get(e_n));
				}

			}

		}
		
		this.reactionsOcurrenceRatios = (LinkedHashMap<String, Double>) kbase.get("reactionsOcurrenceRatios");
		this.reactionPrecedenceRules = (LinkedHashMap<String, Object>) kbase.get("reactionPrecedenceRules");
//		this.mrFilter		= new MReactionsFilter(this.reactionPrecedenceRules);
		
	}

	public ArrayList<Enzyme> getEnzymesList(){
		return this.enzymes;
	}

	
	public SMIRKSManager getSmirksManager(){
		return this.smrkMan;
	}
	
	public  LinkedHashMap<String, MetabolicReaction> getReactionsHash(){
		return this.reactionsHash;
	}
	
	public LinkedHashMap<String, Double> getReactionsORatios(){
		return this.reactionsOcurrenceRatios;
	}


	public  LinkedHashMap<String, Enzyme> getEnzymeHash(){
		return this.enzymesHash;
	}
	
	public  LinkedHashMap<String, ArrayList<Enzyme>> getMetPathwaysHash(){
		return this.metPathwaysHash;
	}
	

}
