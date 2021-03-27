/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */


package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.json.JSONException;
import org.json.JSONObject;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;

public class ChemdbRest {

	// https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
	public static String pubChemCompoundURL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/";

	public ChemdbRest() {
		// TODO Auto-generated constructor stub		
	}

	public static LinkedHashMap<String,ArrayList<String>> getSynonymsObjectViaInChIKey(String inchikey) throws UnirestException{
//		LinkedHashMap<String,ArrayList<String>> results = new LinkedHashMap<String,ArrayList<String>>();
//		results.put("CID", new ArrayList<String>());
//		results.put("Synonym", new ArrayList<String>());
		
		LinkedHashMap<String,ArrayList<String>> results = null;
		
		try{
			Unirest.setTimeouts(4000, 4000);
			String path = pubChemCompoundURL + "inchikey/" + inchikey +"/synonyms/json";
//			System.out.println(path);
			
			HttpResponse<JsonNode> jsonResponse = Unirest.post(path).header("accept", "application/json").asJson();
			JSONObject jObject = jsonResponse.getBody().getObject();
//			System.out.println("KEYS");
//			System.out.println(jObject.keySet().contains("InformationList"));
			
			if(jObject.keySet().contains("InformationList")){
				JSONObject informationList = new JSONObject(jObject.get("InformationList").toString());
				JSONObject information = new JSONObject(informationList.getJSONArray("Information").get(0).toString());
				
				if(information !=null){
					results = new LinkedHashMap<String,ArrayList<String>>();
					results.put("CID", new ArrayList<String>());
					results.get("CID").add(information.get("CID").toString());
					results.put("Synonyms", new ArrayList<String>());

					for(int i = 0; i < information.getJSONArray("Synonym").length(); i++){
						results.get("Synonyms").add(information.getJSONArray("Synonym").get(i).toString());
					}
				}				
			}
			return results;	
		}

		catch(JSONException e){
			System.err.println(e.getLocalizedMessage());
			return results;
		}
		catch(NullPointerException n){
			System.err.println(n.getLocalizedMessage());
			return results;
		}		
		catch(com.mashape.unirest.http.exceptions.UnirestException ue){
			System.err.println(ue.getLocalizedMessage());
			return results;
		}
	}
	
}
