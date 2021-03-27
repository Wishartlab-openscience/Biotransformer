/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.transformation;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.JsonParser.Feature;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import ambit2.smarts.SMIRKSManager;
import biotransformer.utils.JsonUtils;

public class MReactionSets {

	public static ArrayList<MetabolicReaction> standardizationReactions;
	protected SMIRKSManager smrkMan	;

	protected ObjectMapper mapper;
	
	public MReactionSets() throws JsonParseException, JsonMappingException, FileNotFoundException, IOException {
		// TODO Auto-generated constructor stub
		standardizationReactions = new ArrayList<MetabolicReaction>();
		mapper = new ObjectMapper();
		mapper.configure(Feature.ALLOW_COMMENTS, true);
		mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);		
		
		smrkMan	= new SMIRKSManager(SilentChemObjectBuilder.getInstance());
		
		LinkedHashMap<Object,Object> standardizationR = JsonUtils.ingestJson("database/standardizationReactions.json", this.mapper);
		addStandardizationReactions(standardizationR);
//		System.out.println("Number of standardizaton reactions: " + this.standardizationReactions.size());
//		this.standardizationReactions.get(0).display();
	}
	
	
	protected void addStandardizationReactions(LinkedHashMap<Object,Object> standardizationR) {
		if(standardizationR.size()>1 &&  standardizationR.containsKey((Object)"reactions")) {
			LinkedHashMap<String, Object> reactions = (LinkedHashMap<String, Object>) standardizationR.get("reactions");
			
			for(Entry<String, Object> react : reactions.entrySet() ) {
				LinkedHashMap<String, Object> v = (LinkedHashMap<String, Object>) react.getValue();
//				ArrayList<String> smarts = (ArrayList<String>) v.get("smarts");
//				System.out.println("v : " + smarts);
				
				MetabolicReaction mr = new MetabolicReaction(
						react.getKey().toString(),
						v.get("commonName").toString(),
						v.get("smirks").toString(),
						(ArrayList<String>) v.get("smarts"),
						new ArrayList<String>(),
						smrkMan						
						);
//				System.out.println("standardizationReactions is null? " + (this.standardizationReactions == null));
				this.standardizationReactions.add(mr);
			}			
		}


		
		
	}


}
