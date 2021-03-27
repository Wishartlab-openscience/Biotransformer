package biotransformer.utils;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.JsonParser.Feature;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

public class JsonUtils {
	
	
	@SuppressWarnings("unchecked")
	public static LinkedHashMap<Object,Object> ingestJson(String filePathname, ObjectMapper mapper) throws JsonParseException, 
		JsonMappingException, FileNotFoundException, IOException{
	
		return ingestJson(filePathname, mapper, true);
	}
	
	@SuppressWarnings("unchecked")
	public static LinkedHashMap<Object,Object> ingestJson(String filePathname, ObjectMapper mapper, boolean reportMissingFile) throws JsonParseException, 
		JsonMappingException, FileNotFoundException, IOException{
		LinkedHashMap<Object,Object> jsonMap = new LinkedHashMap<Object,Object>();
		
//		System.out.println("Ingesting \""+ filePathname + "\"...");
		
		try {
			jsonMap = (LinkedHashMap<Object,Object>) mapper.readValue(new FileInputStream(filePathname), Map.class);
		}
		catch(FileNotFoundException f) {
			if(reportMissingFile == true) {
				System.err.println(f.getMessage());				
			}

			
		}
		catch(Exception e){
			System.err.println(e.getMessage());
			System.exit(1);		
		}
		
		return jsonMap;		
	}

	public static void main (String[] args) throws FileNotFoundException, IOException {
		ObjectMapper mapper = new ObjectMapper();
		mapper.configure(Feature.ALLOW_COMMENTS, true);
		mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);
		
		System.out.println(ingestJson("database/enzymes.json", mapper));
	}


}


