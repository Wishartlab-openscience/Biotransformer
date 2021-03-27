/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.similarity.Tanimoto;

import biotransformer.fingerprint.ChemStructureFingerprinter;
import biotransformer.fingerprint.SFingerprint;

public class ChemStructureSimilarity {


	public static double calculateTanimotoWithCustomSmarts(IAtomContainer molecule1, IAtomContainer molecule2, String[] smarts) throws Exception{
		SubstructureFingerprinter sf = new SubstructureFingerprinter(smarts);
		double tanimoto_score = 0;
		
		IBitFingerprint fingerprint1 = sf.getBitFingerprint(molecule1);
		IBitFingerprint fingerprint2 = sf.getBitFingerprint(molecule2);
		
		tanimoto_score  = Tanimoto.calculate(fingerprint1, fingerprint2);

		return tanimoto_score;
	}

	public static double calculateTanimotoWithCustomSmartsAndCounts(IAtomContainer molecule1, IAtomContainer molecule2, String[] smarts) throws Exception{
		SubstructureFingerprinter sf = new SubstructureFingerprinter(smarts);
		double tanimoto_score = 0;
		
		ICountFingerprint fingerprint1 = sf.getCountFingerprint(molecule1);
		ICountFingerprint fingerprint2 = sf.getCountFingerprint(molecule2);
		
		tanimoto_score  = Tanimoto.calculate(fingerprint1, fingerprint2);

		return tanimoto_score;
	}
	
	public static double calculateTanimoto(IAtomContainer molecule1, IAtomContainer molecule2, LinkedHashMap<String, String> queries) throws Exception{

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
		
		double tanimoto_score = 0;
		SFingerprint sfingerprint1 = cs.generateClassyfireFingerprintAsDouble(molecule1, queries);
		SFingerprint sfingerprint2 = cs.generateClassyfireFingerprintAsDouble(molecule2, queries);
		
		tanimoto_score  = Tanimoto.calculate(sfingerprint1.getBitValuesDouble(), sfingerprint2.getBitValuesDouble());
		
		return tanimoto_score;
	}
	
	public static double calculateTanimotoWithMACCS(IAtomContainer molecule1, IAtomContainer molecule2) throws Exception{

		double tanimoto_score = 0;
		MACCSFingerprinter mx = new MACCSFingerprinter();
		IBitFingerprint ifingerprint1 = mx.getBitFingerprint(molecule1);
		IBitFingerprint ifingerprint2 = mx.getBitFingerprint(molecule2);

		
		tanimoto_score  = Tanimoto.calculate(ifingerprint1,ifingerprint2);
		
		return tanimoto_score;
	}	
	
	public static String[] getSmarts() throws Exception{
		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
		LinkedHashMap<String, String> smartsMap = cs.getFingerprintPatterns();
		
		String[] smarts = new String[smartsMap.size()];

		int iterator = 0;
		
		System.out.println("Getting smarts");
		
		
		
		Iterator it = smartsMap.entrySet().iterator();
		while (it.hasNext()) {
		     Map.Entry pair = (Map.Entry)it.next();
		   
		     smarts[iterator] = (String) pair.getValue();
		     iterator++;
		}
		
		return smarts;
	}
	
	public static HashMap<String, Double> tanimotoListCreator(HashMap<String, IAtomContainer> compounds, LinkedHashMap<String, String> queries) throws Exception{
		ArrayList<IAtomContainer> compoundValues = new ArrayList<IAtomContainer>();
		ArrayList<String> compoundNames = new ArrayList<String>();
		HashMap<String, Double> tanimotoList = new HashMap<String, Double>();
		double tanimotoScore;
		
		System.out.println("Calculating Tanimoto");
		
		Iterator it = compounds.entrySet().iterator();
		while (it.hasNext()) {
		     Map.Entry pair = (Map.Entry)it.next();		     
		     compoundValues.add((IAtomContainer) pair.getValue());
		     compoundNames.add((String) pair.getKey());
		    
		}
		
		for (int i = 0; i < compoundValues.size(); i++){

			for (int j = i+1; j < compoundValues.size(); j++){
				tanimotoScore = calculateTanimoto(compoundValues.get(i), compoundValues.get(j), queries);
				tanimotoList.put(compoundNames.get(i) + " - " + compoundNames.get(j), tanimotoScore);
			}
		} 
		
		return tanimotoList;
	}
	
	public static void generateCsv(HashMap<String, Double> tanimotoList, String outputname){
//		System.out.println("In generator");
		try
		{
//			System.out.println("Making File");
		    FileWriter writer = new FileWriter(outputname);
		    
		    Iterator it = tanimotoList.entrySet().iterator();
			while (it.hasNext()) {
			     Map.Entry pair = (Map.Entry)it.next();
			     writer.append(String.valueOf(pair.getKey()));
				 writer.append(',');
				 writer.append((CharSequence) String.valueOf(pair.getValue()));
				 writer.append('\n');	    
			}
				
		    writer.flush();
		    writer.close();
		    System.out.println("CSV populated");
		}
		catch(IOException e)
		{
		     e.printStackTrace();
		} 
	}

}
