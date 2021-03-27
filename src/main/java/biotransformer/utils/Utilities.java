/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.transformation.Biotransformation;

public class Utilities {

	public Utilities() {
		// TODO Auto-generated constructor stub
	
	}
	
	public static ArrayList<String> removeDuplicateStrings(ArrayList<String> listOfStrings){
		ArrayList<String> unique =  new ArrayList<String>();		
		LinkedHashSet<String> set = new LinkedHashSet<String>(listOfStrings);
		unique = new ArrayList<String>(set);
		
		return unique;
	}
	
	public static void print(ArrayList<String> aList){
		for(String i : aList){
			System.out.println(i);
		}
	}

	public static IAtomContainerSet getCDKAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
	}
	
	public static void annotateAtomContainerWithProps(IAtomContainer molecule, LinkedHashMap<Object, Object> props){
		LinkedHashMap<Object, Object> p = new LinkedHashMap<Object, Object>();
		for(Entry<Object, Object> m: props.entrySet()){
			if(m.getKey() != "Synonym"){
				p.put(m.getKey(), m.getValue());
			}
			else {
				p.put(CDKConstants.TITLE,m.getValue());
			}
		}		
//		System.out.println(p);
		molecule.setProperties(p);	
	}


	public static ArrayList<Biotransformation> selectUniqueBiotransformations(ArrayList<Biotransformation> biotransformations){
		ArrayList<Biotransformation> unique_bts = new ArrayList<Biotransformation>();
		for(int i = 0; i < biotransformations.size(); i ++){
			if(!containsBiotransformation(unique_bts, biotransformations.get(i))){
				unique_bts.add(biotransformations.get(i));
			}
		}	
		return unique_bts;
	}
	
	public static boolean containsBiotransformation(ArrayList<Biotransformation> biotransformations, Biotransformation bt){
		boolean inc = false;
		for(int i = 0; i < biotransformations.size(); i ++){
			if(biotransformations.get(i).equals(bt)){
				inc = true;
				break;
			}
		}
		
		
		return inc;
	}
	public static IAtomContainerSet createEmptyAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	}
	
	
	public static String returnFirstCleanSynonym(String[] synonyms){
		String fcs = null;
		Pattern p = Pattern.compile("CHEBI:[0-9]+|UNII-|CHEMBL|ZINC|DB[0-9]+|[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-|[0-9]+-[0-9]+-[0-9]|^AC[0-9]+|%",  Pattern.CASE_INSENSITIVE);
		for(int i=0; i < synonyms.length ; i++){
			Matcher b = p.matcher(synonyms[i]);
//			System.out.println(synonyms[i] + ": " + b.find());
			if(!b.find()){
				fcs = synonyms[i];
//				System.out.println("WORKS: " + synonyms[i]);
				break;				
			}		
		}
		if(fcs == null){
			fcs = synonyms[0];
		}
		return fcs;		
	}

	public static String returnFirstCleanSynonym(ArrayList<String> synonyms){
		String fcs = null;
		// CHEMSPIDER|CID[0-9]+|SID[0-9]+|
		Pattern p = Pattern.compile("CHEBI:[0-9]+|UNII-|CHEMBL|ZINC|DB[0-9]+|[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-|[0-9]+-[0-9]+-[0-9]|^AC[0-9]+|%",  Pattern.CASE_INSENSITIVE);
		for(int i=0; i < synonyms.size() ; i++){
			Matcher b = p.matcher(synonyms.get(i).trim());
//			System.out.println(synonyms.get(i) + ": " + b.find());
			
			if(!b.find()){
				fcs = synonyms.get(i);
//				System.out.println("WORKS: " + synonyms.get(i));
				break;			
			}		
//			System.out.println(synonyms.get(i) + ": " + b.find());
		}
		if(fcs == null){
			fcs = synonyms.get(0);
		}
		return fcs;	
	}
//	public static IAtomContainerSet uniquefy(IAtomContainerSet molecules)
//			throws Exception {
//		if (molecules != null && (!molecules.isEmpty()) && molecules.getAtomContainerCount() > 1) {
//			
//			IAtomContainerSet uniqueContainer = DefaultChemObjectBuilder.getInstance()
//					.newInstance(IAtomContainerSet.class);
//			
//			uniqueContainer.addAtomContainer(molecules.getAtomContainer(0));
//			
//			for (int i = 1; i < molecules.getAtomContainerCount(); i++) {
//				if (! ( (molecules.getAtomContainer(i) == null) || atomContainerInclusionHolds(uniqueContainer,
//						molecules.getAtomContainer(i) ))) {
//					uniqueContainer.addAtomContainer(molecules.getAtomContainer(i));
//				}
//			}
//
//			return uniqueContainer;
//		}
//
//		else
//			return molecules;
//
//	}

	public static void addPhysicoChemicalProperties(IAtomContainer molecule) throws CDKException {
			
			LinkedHashMap<String, String> properties = ChemStructureExplorer.computePhysicoChemicalProperties(molecule);
			
			for(Map.Entry<String, String> prop : properties.entrySet()){
	//			System.out.println(prop.getValue());
				molecule.setProperty(prop.getKey(), String.format("%.8s", Double.valueOf(prop.getValue()))   );
			}
			
		}

}
