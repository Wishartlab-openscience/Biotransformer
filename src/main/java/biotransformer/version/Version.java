package biotransformer.version;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.DocumentBuilder;  
import org.w3c.dom.Document;  
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import java.io.File;
import java.io.IOException;

public class Version {
	public static String current; 
	
	public Version() throws ParserConfigurationException, SAXException, IOException {
		current = getVersion();
	}
	
	private static String getVersion() throws ParserConfigurationException, SAXException, IOException {
//		String version = null;		
		File pomfile = new File("pom.xml");
		
		//an instance of factory that gives a document builder  
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();  
		//an instance of builder to parse the specified xml file  
		DocumentBuilder db = dbf.newDocumentBuilder();  
		Document doc = db.parse(pomfile);  
		doc.getDocumentElement().normalize();	
		NodeList nodeList = doc.getElementsByTagName("version");
//		System.out.println("Version : " + nodeList.item(0).getTextContent());
//		System.out.println(doc.getElementsByTagName("artifactId").item(0).getTextContent());
//		System.out.println(doc.getElementsByTagName("repositories").item(0).getTextContent());
//	
		return nodeList.item(0).getTextContent();
	}
	
//	public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException {
//		getVersion();
//	}
}
