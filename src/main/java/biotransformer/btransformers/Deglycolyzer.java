/**
 * This class implements the class of metabolic enzymes.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */
package biotransformer.btransformers;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;

import biotransformer.biosystems.BioSystem.BioSystemName;
import exception.BioTransformerException;

/**
 * @author Yannick Djoumbou Feunang
 *
 */
public class Deglycolyzer extends Biotransformer {

	/**
	 * @throws ParseException 
	 * @throws IOException 
	 * @throws CDKException 
	 * 
	 */
	public Deglycolyzer(BioSystemName bioSName) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException {
		super(bioSName);
		// TODO Auto-generated constructor stub
	}

}
