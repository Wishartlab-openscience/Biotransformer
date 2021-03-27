/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.biosystems;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

import exception.BioTransformerException;

public class EnvMicro extends BioSystem {

	public EnvMicro(ObjectMapper mapper) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException{
		super(BioSystemName.ENVMICRO, mapper);
	}

}
