
/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.codehaus.jackson.JsonParser.Feature;
import org.codehaus.jackson.map.ObjectMapper;
import org.json.simple.JSONValue;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biosystems.BioSystem;
import biotransformer.biosystems.BioSystem.BioSystemName;


public class ErdbTask {

	public ErdbTask() {
		// TODO Auto-generated constructor stub
	
	}

	
	public void filter() throws IOException{
		BufferedReader bReader = new BufferedReader (new FileReader("data/CYP_3D_annotated_metabolites_withAllFeats_21012017_new.tsv"));

		ArrayList<String> lipinskiCompounds = new ArrayList<String>();
		ArrayList<String> leadLikeCompounds = new ArrayList<String>();
		ArrayList<String> lipinskiOrleadLikeCompounds = new ArrayList<String>();
		ArrayList<String> filtered = new ArrayList<String>();
		
		String line;
		int counter = 0;
		int reactantNr = 0;
		boolean lipinski, leadLike;
				
		// nHBAcc 15
		// nHBDon 16
		// nB 19
		// nRotB 21
		// MLogP 26
		// TopoPSA 28
		// MW 29
		
		while((line = bReader.readLine()) != null){
			counter++;
			System.out.println(counter);
			
			if(counter > 1 && counter < 501){
				String[] sline = line.split("\t");
				System.out.println(sline[29].getClass());
				if(sline[4] == "R" || sline[5] == "R" || sline[6] == "R" || sline[7] == "R"|| sline[8] == "R"
						|| sline[9] == "R"|| sline[10] == "R" || sline[11] == "R" || sline[12] == "R"){
					filtered.add(line);
					reactantNr++;
				} else{
					lipinski = (Integer.valueOf(sline[29]) < 500 && Integer.valueOf(28) < 5 &&
							Integer.valueOf(16) < 5 && Integer.valueOf(15) < 10 );
					leadLike = (Integer.valueOf(sline[29]) < 300 && Integer.valueOf(28) < 3 &&
							Integer.valueOf(16) < 3 && Integer.valueOf(3) < 10 &&
							Integer.valueOf(21) < 3);
					
					if(lipinski && leadLike){
						lipinskiCompounds.add(line);
						leadLikeCompounds.add(line);
						lipinskiOrleadLikeCompounds.add(line);
						filtered.add(line);
					} 
					else if(lipinski){
						lipinskiCompounds.add(line);
						lipinskiOrleadLikeCompounds.add(line);
						filtered.add(line);
						} 
					else if(leadLike){
						leadLikeCompounds.add(line);
						lipinskiOrleadLikeCompounds.add(line);
						filtered.add(line);
					}
	
				}	
			}			
		}
		
		System.out.println("lipinskiCompounds: " + lipinskiCompounds.size() );
		System.out.println("leadLikeCompounds: " + leadLikeCompounds.size() );
		System.out.println("lipinskiOrleadLikeCompounds: " + lipinskiOrleadLikeCompounds.size() );
		System.out.println("filtered: " + filtered.size() );
		System.out.println("reactantNr: " + reactantNr);
		
		bReader.close();
	}
	
	public void generateReactantsSetFromTSV() throws IOException{
		BufferedReader bReader = new BufferedReader (new FileReader("data/CYP_3D_annotated_metabolites_withAllFeats_24022017_ERDB_v60_.tsv"));
		ArrayList<String> reactantsOrInhibitors = new ArrayList<String>();
		ArrayList<String> nonReactantsToAllCPs = new ArrayList<String>();
		
		ArrayList<Double> massesForRIs = new ArrayList<Double>();
		ArrayList<Double> massesForNonReactantsToAllCPs = new ArrayList<Double>();
		String maxMassRIsSmiles = "";
		
		ArrayList<Double> tpsaForRIs = new ArrayList<Double>();
		ArrayList<Double> tpsaForNonReactantsToAllCPs = new ArrayList<Double>();
		String maxTpsaRIsSmiles = "";
		
		
		ArrayList<Double> asaForRIs = new ArrayList<Double>();
		ArrayList<Double> asaForNonReactantsToAllCPs = new ArrayList<Double>();
		String maxAsaRIsSmiles = "";
		
		ArrayList<Double> mLogForRIs = new ArrayList<Double>();
		ArrayList<Double> mLogPForNonReactantsToAllCPs = new ArrayList<Double>();
		String maxMLogPaRIsSmiles = "";		
		
		String line;
		double mass=0.0;
		double tpsa=0.0;
		double asa=0.0;
		double mlogp = 0.0;
		
		int counter = 0;
		
		while((line = bReader.readLine()) != null){
			counter++;
			
			if(counter>1){
			String[] sline = line.split("\t");
			
//			if(counter>1){
//			System.out.println(Double.valueOf(sline[28]));
//			if(Double.valueOf(sline[29])>mass){				
//				mass = Double.valueOf(sline[29]);
//				System.out.println(Double.valueOf(sline[29]));
//			}
//			
//			if(Double.valueOf(sline[28])>tpsa){				
//				tpsa = Double.valueOf(sline[28]);
//				System.out.println(Double.valueOf(sline[28]));
//			}
//				System.out.println(Double.valueOf(sline[28]));
//				if(Double.valueOf(sline[29])>mass){				
//					mass = Double.valueOf(sline[29]);
//					System.out.println(Double.valueOf(sline[29]));
//				}
//				
//				if(Double.valueOf(sline[28])>tpsa){				
//					tpsa = Double.valueOf(sline[28]);
//					System.out.println(Double.valueOf(sline[28]));
//				}
//			}
			
			if(sline[4].contains("N") && sline[5].contains("N") && sline[6].contains("N") && sline[7].contains("N") &&sline[8].contains("N")
					&& sline[9].contains("N") && sline[10].contains("N") && sline[4].contains("N") && sline[12].contains("N")){
				nonReactantsToAllCPs.add(line);
				tpsaForNonReactantsToAllCPs.add(Double.valueOf(sline[28]));
				mLogPForNonReactantsToAllCPs.add(Double.valueOf(sline[26]));
				

				if((!sline[50].contains("NaN")) && (sline[50] != null)){
					asaForNonReactantsToAllCPs.add(Double.valueOf(sline[50]));
				}
				
				if(sline[29] != null && Double.valueOf(sline[29]) > 0.0){
					massesForNonReactantsToAllCPs.add(Double.valueOf(sline[29]));
				} else
					System.out.println(sline[29] + " || " + Double.valueOf(sline[29]));

			} else
			
			if(sline[4].contains("R") ||sline[5].contains("R") || sline[6].contains("R") || sline[7].contains("R")|| sline[8].contains("R")
					|| sline[9].contains("R")|| sline[10].contains("R") || sline[11].contains("R") || sline[12].contains("R") ||
					
					sline[4].contains("I") ||sline[5].contains("I") || sline[6].contains("I") || sline[7].contains("I")|| sline[8].contains("I")
					|| sline[9].contains("I")|| sline[10].contains("I") || sline[11].contains("I") || sline[12].contains("I")){
				
				reactantsOrInhibitors.add(line);
				tpsaForRIs.add(Double.valueOf(sline[28]));
				mLogForRIs.add(Double.valueOf(sline[26]));
				
				if((!sline[50].contains("NaN")) && (sline[50] != null)){
					asaForRIs.add(Double.valueOf(sline[50]));
				}
				
				
				if(sline[29] != null && Double.valueOf(sline[29]) > 0.0){
					massesForRIs.add(Double.valueOf(sline[29]));	
				} else {
					System.out.println(sline[29] + " || " + Double.valueOf(sline[29]));
				}
				
				
//				System.out.print(sline[4] + "\t" + sline[4] +  "\t" +sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] +  "\t" +sline[9] + "\t" + sline[10] + "\t" + sline[11]+ "\t" + sline[12] + "\n");
//				System.out.println(Double.valueOf(sline[28]));
			
				if(Double.valueOf(sline[26])>mlogp){	
//				if((!sline[14].contains("NaN")) && (sline[14] != null) && Double.valueOf(sline[14])>mlogp){				
					mlogp = Double.valueOf(sline[26]);
					maxMLogPaRIsSmiles = sline[14] + "\n" + sline[4] +  "\t" +sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] +  "\t" +sline[9] + "\t" + sline[10] + "\t" + sline[11]+ "\t" + sline[12];
//					System.out.println(Double.valueOf(sline[29]));
				}
				
				if(Double.valueOf(sline[29])>mass){				
					mass = Double.valueOf(sline[29]);
					maxMassRIsSmiles = sline[14] + "\n" + sline[4] +  "\t" +sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] +  "\t" +sline[9] + "\t" + sline[10] + "\t" + sline[11]+ "\t" + sline[12];
//					System.out.println(Double.valueOf(sline[29]));
				}
				
				if(Double.valueOf(sline[28])>tpsa){				
					tpsa = Double.valueOf(sline[28]);
					maxTpsaRIsSmiles = sline[14] + "\n" + sline[4] + "\t" +sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] +  "\t" +sline[9] + "\t" + sline[10] + "\t" + sline[11]+ "\t" + sline[12];
//					System.out.println(Double.valueOf(sline[28]));
				}
				
				if((!sline[50].contains("NaN")) && (sline[50] != null) && Double.valueOf(sline[50])>asa){				
					asa = Double.valueOf(sline[50]);
					maxAsaRIsSmiles = sline[14] + "\n" + sline[4] + "\t" + sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] +  "\t" +sline[9] + "\t" + sline[10] + "\t" + sline[11]+ "\t" + sline[12];;
//					System.out.println(Double.valueOf(sline[28]));
				}
				
			}
			
		}
	}
		System.out.println("MAX mass: " + mass);
		System.out.println("MAX tpsa: " + tpsa);
		System.out.println("massesForReactantsOrInhibitors (MIN): " + Collections.min(massesForRIs));
		System.out.println("massesForReactantsOrInhibitorss (MAX): " + Collections.max(massesForRIs));	
		System.out.println("massesForNonReactantsToAllCPs (MIN): " + Collections.min(massesForNonReactantsToAllCPs));
		System.out.println("massesForNonReactantsToAllCPs (MAX): " + Collections.max(massesForNonReactantsToAllCPs));
		System.out.println("\n");
		System.out.println("tpsaForReactantsOrInhibitors (MIN): " + Collections.min(tpsaForRIs));
		System.out.println("tpsaForReactantsOrInhibitorss (MAX): " + Collections.max(tpsaForRIs));	
		System.out.println("tpsaForNonReactantsToAllCPs (MIN): " + Collections.min(tpsaForNonReactantsToAllCPs));
		System.out.println("tpsaForNonReactantsToAllCPs (MAX): " + Collections.max(tpsaForNonReactantsToAllCPs));
		System.out.println("\n");
		System.out.println("asaForReactantsOrInhibitors (MIN): " + Collections.min(asaForRIs));
		System.out.println("asaForReactantsOrInhibitorss (MAX): " + Collections.max(asaForRIs));	
		System.out.println("asaForNonReactantsToAllCPs (MIN): " + Collections.min(asaForNonReactantsToAllCPs));
		System.out.println("asaForNonReactantsToAllCPs (MAX): " + Collections.max(asaForNonReactantsToAllCPs));
		System.out.println("\n");
		System.out.println("mLogPForReactantsOrInhibitors (MIN): " + Collections.min(mLogForRIs));
		System.out.println("mLogPForReactantsOrInhibitorss (MAX): " + Collections.max(mLogForRIs));	
		System.out.println("mLogPForNonReactantsToAllCPs (MIN): " + Collections.min(mLogPForNonReactantsToAllCPs));
		System.out.println("mLogPForNonReactantsToAllCPs (MAX): " + Collections.max(mLogPForNonReactantsToAllCPs));
		System.out.println("\n");
		
		System.out.println("maxMassRIsSmiles: " + maxMassRIsSmiles);
		System.out.println("maxTpsaRIsSmiles: " + maxTpsaRIsSmiles);
		System.out.println("maxAsaRIsSmiles: " + maxAsaRIsSmiles);
		System.out.println("maxMLogPaRIsSmiles: " + maxMLogPaRIsSmiles);
		
		bReader.close();
	}
	
	
	public void generateBioTransformerDB() throws Exception{
		BufferedReader bReader = new BufferedReader (new FileReader("/Users/yandj/Programming/Projects/Metabolism/cyp450_metabolites.tsv"));
		BufferedReader bReader_2 = new BufferedReader (new FileReader("/Users/yandj/Programming/Projects/Metabolism/phase2_metabolites.tsv"));
		BufferedReader bReader_3 = new BufferedReader (new FileReader("/Users/yandj/Programming/Projects/Metabolism/GutMicrobialBiotransformationsAndMetabolites.tsv"));
		BufferedWriter bw0 = new BufferedWriter(new FileWriter("/Users/yandj/Programming/Projects/Metabolism/btdb_metabolites.json"));
		
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		LinkedHashMap<String, ArrayList<LinkedHashMap<String, Object>>> lmao = new LinkedHashMap<String, ArrayList<LinkedHashMap<String, Object>>>();
		LinkedHashMap<String, String> inchikeyToID = new LinkedHashMap<String, String>();
		LinkedHashMap<String, LinkedHashMap<String, Object>> cpdDict = new LinkedHashMap<String, LinkedHashMap<String, Object>>();
		
		LinkedHashMap<String, LinkedHashMap<String, Object>> lmaof = new LinkedHashMap<String, LinkedHashMap<String, Object>>();
		ObjectMapper mapper = new ObjectMapper();
		mapper.configure(Feature.ALLOW_COMMENTS, true);
		mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);
		BioSystem bsys = new BioSystem(BioSystemName.HUMAN,mapper);
		BioSystem bsys2 = new BioSystem(BioSystemName.GUTMICRO,mapper);
		
//		System.err.println(bsys.getReactionsHash());
//		System.err.println(bsys.getReactionsHash().get("HYDROXYLATION_OF_ALIPHATIC_SECONDARY_ANTEPENULTIMATE_CARBON_PATTERN1").commonName);
		
		 Pattern drugBankPattern = Pattern.compile("^DB[0-9]|^DBMET[0-9]");
		 Pattern hmdbPattern = Pattern.compile("^^HMDB[0-9]");
		
		int id = 1;
		// GET CYP450 TRANSFORMATIONS
		int btCounter = 0;
		String line	 = null;
		String line2 = null;
		String line3 = null;
		
		while((line = bReader.readLine()) != null && (line = bReader.readLine()).trim() != "" 
				&& (line = bReader.readLine()).trim() != "\n"){			
			String[] sline = line.split("\t");
//			System.out.println(line);
			
			if(sline.length>12 && sline[1] != null && sline[1].trim() != "" && (!sline[1].contains("Structure"))
					&& sline[7] != null && sline[7].length()>0 && sline[8] != null && sline[9] != null && sline[9].contains("CYP")
					&& sline[10] !=null && sline[10].trim().length()>0 && sline[11] != null && sline[11].trim().length()>0
					){
				
				String[] pname = sline[10].split(";");
				String[] pstruct = sline[11].split(Pattern.quote("."));
				System.err.println("sline[0] : " + sline[0]);
				System.err.println("sline[10] : " + sline[10]);
				System.err.println("sline[11] : " + sline[11]);	
				
				if(pname != null && pstruct != null && pstruct.length > 0 && pname.length == pstruct.length){

					System.out.println("GOOD LINE");
					IAtomContainer atc = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(sline[1].trim(), true).get("atomContainer");
					atc= AtomContainerManipulator.suppressHydrogens(atc); 
					String inchikey = inchiGenFactory.getInChIGenerator(atc).getInchiKey();
					
					if(!cpdDict.containsKey(inchikey)){
						cpdDict.put(inchikey, new LinkedHashMap<String, Object>());
						cpdDict.get(inchikey).put("Name", sline[0].trim());
						cpdDict.get(inchikey).put("SMILES", ChemStructureExplorer.smiGen.create(atc));
						cpdDict.get(inchikey).put("InChIKey", inchikey);
						cpdDict.get(inchikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
						cpdDict.get(inchikey).put("PubChem CID","NULL");
						cpdDict.get(inchikey).put("DrugBank ID", "NULL");
						cpdDict.get(inchikey).put("HMDB_ID", "NULL");
						
						LinkedHashMap<String,ArrayList<String>> syn = ChemdbRest.getSynonymsObjectViaInChIKey(inchikey);

						if(syn != null && syn.get("CID") != null){						
							cpdDict.get(inchikey).put("PubChem CID", syn.get("CID").get(0));

							for(String s : syn.get("Synonyms")){
								Matcher m = drugBankPattern.matcher(s);
								Matcher n = hmdbPattern.matcher(s);								
								if(m.find()){
									cpdDict.get(inchikey).put("DrugBank ID", s);
								}
								if(n.find()){
									cpdDict.get(inchikey).put("HMDB_ID", s);
								}					
							}						
						}
						
						id++;
					}
									
					LinkedHashMap<String, Object> substrate = new LinkedHashMap<String, Object>();
					substrate = cpdDict.get(inchikey);
					
					for(String i : sline[7].split(";")){

						LinkedHashMap<String, Object> bt = new LinkedHashMap<String, Object>();
						bt.put("Substrate",substrate);
						bt.put("Enzyme(s)", sline[9].replace(" (minor)", "").replace(" (major)", ""));

						System.out.println("Reaction");
						System.out.println('"' + i.trim() + '"');
						bt.put("Reaction Type", bsys.getReactionsHash().get(i.trim()).commonName);
						bt.put("BioTransformer Reaction ID (BTMRID)", bsys.getReactionsHash().get(i.trim()).reactionsBTMRID);
						bt.put("Biotransformation type", "Human Phase I");
						bt.put("Biosystem", "Human");
								
//						if(pname != null && pstruct != null && pstruct.length > 0 && pname.length == pstruct.length){
							btCounter++;
							
							for(int k = 0; k < pname.length; k++){
								ArrayList<LinkedHashMap<String, Object>> products = new ArrayList<LinkedHashMap<String, Object>>();
								System.err.println("K : " + k);						
								System.err.println("pstruct[k] : " + pstruct[k].trim());
								IAtomContainer a = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(pstruct[k].trim(), true).get("atomContainer");
								a= AtomContainerManipulator.suppressHydrogens(a); 
								String ikey = inchiGenFactory.getInChIGenerator(a).getInchiKey();
								
								if(cpdDict.containsKey(ikey)){
									products.add(cpdDict.get(ikey));
								}else{
									cpdDict.put(ikey, new LinkedHashMap<String, Object>());
									cpdDict.get(ikey).put("Name", pname[k].trim());
									cpdDict.get(ikey).put("SMILES", ChemStructureExplorer.smiGen.create(a));
									cpdDict.get(ikey).put("InChIKey", ikey);
									cpdDict.get(ikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
									cpdDict.get(ikey).put("PubChem CID","NULL");
									cpdDict.get(ikey).put("DrugBank ID", "NULL");
									cpdDict.get(ikey).put("HMDB_ID", "NULL");

									LinkedHashMap<String,ArrayList<String>> synp = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
									
									if(synp != null && synp.get("CID") != null){						
										cpdDict.get(ikey).put("PubChem CID", synp.get("CID").get(0));

										for(String s : synp.get("Synonyms")){
											Matcher m = drugBankPattern.matcher(s);
											Matcher n = hmdbPattern.matcher(s);								
											if(m.find()){
												cpdDict.get(ikey).put("DrugBank ID",s);
											}
											if(n.find()){
												cpdDict.get(ikey).put("HMDB_ID",s);
											}					
										}						
									}
									
									products.add(cpdDict.get(ikey));
									id++;
								}
								
								bt.put("Products", products);
							}

							
							if(sline.length>=14){
	//							bt.put("References", sline[15].replace("(R1)","").replace("(R2)","").replace("(R3)","").replace("(R4)","").trim().replace(" || ", "\n").replace("||", "\n"));						
								bt.put("References", sline[15]);
	//							bt.put("References", String.join("\n", sline[15].replace("(R1)","").replace("(R2)","").replace("(R3)","").replace("(R4)","").trim().split("||")));
								
							}
							
							lmaof.put("BIOTID" + String.format("%04d", btCounter), bt);
	//						btCounter++;
	//						System.err.println(v);	
							
//						}		
					}
									
				} else{
					System.err.println("Check the number of product names and structures on line\n" + line + "\n"+ String.valueOf(pname.length) + " names vs. " + String.valueOf(pstruct.length) + " structures." );
					break;
				}
			}
		}
		
		while( (line2 = bReader_2.readLine()) != null && (line2 = bReader_2.readLine()).trim() != "" 
				&& (line2 = bReader_2.readLine()).trim().length()>0){
			String[] sline2 = line2.split("\t");
			
			System.err.println(line2);
			
			if(sline2.length>=12 && sline2[1] != null && sline2[2] != null && sline2[9] != null && sline2[6] != null && sline2[7] != null){
				LinkedHashMap<String, Object> bt = new LinkedHashMap<String, Object>();
				LinkedHashMap<String, Object> substrate = new LinkedHashMap<String, Object>();
				btCounter++;
				String drugbankID = "NULL";
				String hmdbID = "NULL";
				
				if(cpdDict.containsKey(sline2[2])){				
					substrate = cpdDict.get(sline2[2]);										
				} else{
					IAtomContainer atc = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(sline2[1].trim(), true).get("atomContainer");
					atc= AtomContainerManipulator.suppressHydrogens(atc); 
					String inchikey = inchiGenFactory.getInChIGenerator(atc).getInchiKey();
					
					cpdDict.put(inchikey, new LinkedHashMap<String, Object>());
					cpdDict.get(inchikey).put("Name", sline2[0].trim());
					cpdDict.get(inchikey).put("SMILES", ChemStructureExplorer.smiGen.create(atc));
					cpdDict.get(inchikey).put("InChIKey", inchikey);
					cpdDict.get(inchikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
					cpdDict.get(inchikey).put("PubChem CID","NULL");
					cpdDict.get(inchikey).put("DrugBank ID","NULL");
					cpdDict.get(inchikey).put("HMDB_ID","NULL");
					
					LinkedHashMap<String,ArrayList<String>> syn = ChemdbRest.getSynonymsObjectViaInChIKey(sline2[2]);

					if(syn != null && syn.get("CID") != null){						
						cpdDict.get(inchikey).put("PubChem CID",syn.get("CID").get(0));
						for(String s : syn.get("Synonyms")){
							Matcher m = drugBankPattern.matcher(s);
							Matcher n = hmdbPattern.matcher(s);								
							if(m.find()){
								cpdDict.get(inchikey).put("DrugBank ID",s);
							}
							if(n.find()){
								cpdDict.get(inchikey).put("HMDB_ID",s);
							}					
						}						
					}
					
					substrate = cpdDict.get(inchikey);
					id++;				
				}
				

				
				bt.put("Substrate", substrate);
				bt.put("Enzyme(s)", sline2[11].replace(" (minor)", "").replace(" (major)", ""));
				bt.put("Reaction Type", sline2[6]);
				bt.put("BioTransformer Reaction ID (BTMRID)", sline2[7]);
				bt.put("Biotransformation type", "Human Phase II");
				bt.put("Biosystem", "Human");				
				
				ArrayList<LinkedHashMap<String, Object>> products = new ArrayList<LinkedHashMap<String, Object>>();
				System.err.println(sline2[9]);
				
				IAtomContainer a = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(sline2[9].trim(), true).get("atomContainer");
				a = AtomContainerManipulator.suppressHydrogens(a);
				String ikey = inchiGenFactory.getInChIGenerator(a).getInchiKey();				

				if(cpdDict.containsKey(ikey)){
					products.add(cpdDict.get(ikey));
				}else{
					cpdDict.put(ikey, new LinkedHashMap<String, Object>());
					cpdDict.get(ikey).put("Name",sline2[8].trim());
					cpdDict.get(ikey).put("SMILES", ChemStructureExplorer.smiGen.create(a));
					cpdDict.get(ikey).put("InChIKey", ikey);
					cpdDict.get(ikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
					cpdDict.get(ikey).put("PubChem CID","NULL");
					cpdDict.get(ikey).put("DrugBank ID","NULL");
					cpdDict.get(ikey).put("HMDB_ID","NULL");					

					LinkedHashMap<String,ArrayList<String>> syn2 = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
					
					syn2 = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
					if(syn2 != null && syn2.get("CID") != null){						
						cpdDict.get(ikey).put("PubChem CID", syn2.get("CID").get(0));

						for(String s : syn2.get("Synonyms")){
							Matcher m = drugBankPattern.matcher(s);
							Matcher n = hmdbPattern.matcher(s);								
							if(m.find()){
								cpdDict.get(ikey).put("DrugBank ID",s);
							}
							if(n.find()){
								hmdbID = s;
							}					
						}						
					}
										
					products.add(cpdDict.get(ikey));
					id++;
				}
				
				bt.put("Products", products);

				if(sline2.length>=13){
					bt.put("References", sline2[12]);
//							bt.put("References", String.join("\n", sline[15].replace("(R1)","").replace("(R2)","").replace("(R3)","").replace("(R4)","").trim().split("||")));
					
				}
				
				String v = "BIOTID" + String.format("%04d", btCounter);
				lmaof.put(v, bt);
				
				System.err.println(v);					
				
			}
		}
		
		
		while((line3 = bReader_3.readLine()) != null && (line3 = bReader_3.readLine()).trim() != "" 
				&& (line3 = bReader_3.readLine()).trim() != "\n"){			
			String[] sline3 = line3.split("\t");
//			System.out.println(line3);
			
			if(sline3.length>12 && sline3[1] != null && sline3[1].trim() != "" && (!sline3[1].contains("Structure"))
					&& sline3[7] != null && sline3[7].length()>0 && sline3[8] != null && sline3[9] != null
					&& sline3[10] !=null && sline3[10].trim().length()>0 && sline3[11] != null && sline3[11].trim().length()>0
					){
				
				String[] pname = sline3[10].split(";");
				String[] pstruct = sline3[11].split(Pattern.quote("."));
				System.err.println("sline[0] : " + sline3[0]);
				System.err.println("sline[10] : " + sline3[10]);
				System.err.println("sline[11] : " + sline3[11]);	
				
				if(pname != null && pstruct != null && pstruct.length > 0 && pname.length == pstruct.length){

					System.out.println("GOOD LINE");
					IAtomContainer atc = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(sline3[1].trim(), true).get("atomContainer");
					atc= AtomContainerManipulator.suppressHydrogens(atc); 
					String inchikey = inchiGenFactory.getInChIGenerator(atc).getInchiKey();

					
					if(!cpdDict.containsKey(inchikey)){
						cpdDict.put(inchikey, new LinkedHashMap<String, Object>());
						cpdDict.get(inchikey).put("Name", sline3[0].trim());
						cpdDict.get(inchikey).put("SMILES", ChemStructureExplorer.smiGen.create(atc));
						cpdDict.get(inchikey).put("InChIKey", inchikey);
						cpdDict.get(inchikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
						cpdDict.get(inchikey).put("PubChem CID","NULL");
						cpdDict.get(inchikey).put("DrugBank ID","NULL");
						cpdDict.get(inchikey).put("HMDB_ID","NULL");
						
						LinkedHashMap<String,ArrayList<String>> syn = ChemdbRest.getSynonymsObjectViaInChIKey(inchikey);
						if(syn != null && syn.get("CID") != null){						
							cpdDict.get(inchikey).put("PubChem CID", syn.get("CID").get(0));

							for(String s : syn.get("Synonyms")){
								Matcher m = drugBankPattern.matcher(s);
								Matcher n = hmdbPattern.matcher(s);								
								if(m.find()){
									cpdDict.get(inchikey).put("DrugBank ID",s);
								}
								if(n.find()){
									cpdDict.get(inchikey).put("HMDB_ID",s);
								}					
							}						
						}
						
						id++;
					}
					LinkedHashMap<String, Object> substrate = new LinkedHashMap<String, Object>();
					substrate = cpdDict.get(inchikey);
					
					for(String i : sline3[7].split(";")){

						


						LinkedHashMap<String, Object> bt = new LinkedHashMap<String, Object>();
						bt.put("Substrate",substrate);
						bt.put("Enzyme(s)", sline3[9].replace(" (minor)", "").replace(" (major)", ""));

						System.out.println("Reaction");
						System.out.println('"' + i.trim() + '"');
						bt.put("Reaction Type", bsys2.getReactionsHash().get(i.trim()).commonName);
						bt.put("BioTransformer Reaction ID (BTMRID)", bsys2.getReactionsHash().get(i.trim()).reactionsBTMRID);
						bt.put("Biotransformation type", "Human Gut Microbial");
						bt.put("Biosystem", "Human");
								
//						if(pname != null && pstruct != null && pstruct.length > 0 && pname.length == pstruct.length){
							btCounter++;
							
							for(int k = 0; k < pname.length; k++){
								ArrayList<LinkedHashMap<String, Object>> products = new ArrayList<LinkedHashMap<String, Object>>();
								System.err.println("K : " + k);						
								System.err.println("pstruct[k] : " + pstruct[k].trim());
								IAtomContainer a = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmiles(pstruct[k].trim(), true).get("atomContainer");
								a= AtomContainerManipulator.suppressHydrogens(a);
								String ikey = inchiGenFactory.getInChIGenerator(a).getInchiKey();
								
								if(cpdDict.containsKey(ikey)){
									products.add(cpdDict.get(ikey));
								}else{
									cpdDict.put(ikey, new LinkedHashMap<String, Object>());
									cpdDict.get(ikey).put("Name", pname[k].trim());
									cpdDict.get(ikey).put("SMILES", ChemStructureExplorer.smiGen.create(a));
									cpdDict.get(ikey).put("InChIKey", ikey);
									cpdDict.get(ikey).put("BTMDB_ID", "BTM" + String.format("%04d", id));
									cpdDict.get(ikey).put("PubChem CID","NULL");
									cpdDict.get(ikey).put("DrugBank ID","NULL");
									cpdDict.get(ikey).put("HMDB_ID","NULL");

									LinkedHashMap<String,ArrayList<String>> syn = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
									if(syn != null && syn.get("CID") != null){						
										cpdDict.get(ikey).put("PubChem CID", syn.get("CID").get(0));

										for(String s : syn.get("Synonyms")){
											Matcher m = drugBankPattern.matcher(s);
											Matcher n = hmdbPattern.matcher(s);								
											if(m.find()){
												cpdDict.get(ikey).put("DrugBank ID",s);
											}
											if(n.find()){
												cpdDict.get(ikey).put("HMDB_ID", s);
											}					
										}						
									}
									
									products.add(cpdDict.get(ikey));
									id++;
								}
								
								bt.put("Products", products);
							}
							
							if(sline3.length>=14){
								bt.put("References", sline3[15]);
							}
							
							lmaof.put("BIOTID" + String.format("%04d", btCounter), bt);
							
//						}		
					}
									
				} else{
					System.err.println("Check the number of product names and structures on line\n" + line + "\n"+ String.valueOf(pname.length) + " names vs. " + String.valueOf(pstruct.length) + " structures." );
					break;
				}
			}
		}
		
//		System.out.println(lmaof);
		String jsonText = JSONValue.toJSONString(lmaof);
		System.out.print(jsonText);
		bw0.write(jsonText);
		bw0.close();
		
		bReader.close();

	}
	
	
	
	public static void main(String[] args) throws Exception{
		
//		LinkedHashMap<String,IAtomContainer> map=new LinkedHashMap<String,IAtomContainer>();
//		FileWriter fw = new FileWriter("data/addedSetAnnotated.tsv");
//		BufferedWriter bw = new BufferedWriter(fw);
//		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
//		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader("data/AddedSet.sdf"), bldr);
//				
//		bw.write("InChiKey\tPubChemID\tHMDB\tDrugBank\t1A2\t2A6\t2B6\t2C8\t2C9\t2C19\t2D6\t2|E1\t3A4\tName\tIsomericSmiles");
//		bw.write("\tnHBAcc\tnHBDon\tnaAromAtom\tnAtomP\tnB\tnAromBond\tnRotB\tALogP\tALogp2\tAMR\tXLogP\tMLogP\tapol\tTopoPSA\tMW\tbpol\tATSc1\tATSc2\tATSc3\tATSc4\tATSc5\tATSm1\tATSm2\tATSm3\tATSm4\tATSm5\tnAcid\tnBase\tMOMI-X\tMOMI-Y\tMOMI-Z\tMOMI-XY\tMOMI-XZ\tMOMI-YZ\tMOMI-R\tAllSurfaceArea");
//
//		LinkedHashMap<String, String> fpatterns = ChemStructureFingerprinter.getRINFingerprintPatterns();
//		String[] labels = fpatterns.keySet().toArray(new String[fpatterns.size()]);
//		bw.write("\t" + StringUtils.join(labels,"\t"));
//		
//		for(int h = 0; h < 881; h++){
//			
//			bw.write(String.format("\tpubchem_f%03d", h+1));
//		}
//
//		for(int h = 0; h < 166; h++){
//			
//			bw.write(String.format("\tmaccs_k%03d", h+1));
//		}
//
//		bw.newLine();
//		
//		SmilesGenerator sg = new SmilesGenerator().isomeric();
//		FeatureGenerator fgen = new FeatureGenerator();
//		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
//		PubchemFingerprinter pbf 	= new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
//		MACCSFingerprinter maccs 	=  new MACCSFingerprinter(SilentChemObjectBuilder.getInstance());
//		ChemStructureManipulator bt = new ChemStructureManipulator();
//
//		
//		while (sdfr.hasNext()) {
//			IAtomContainer container = sdfr.next();
//			ArrayList<String> line = new ArrayList<String>();
//			line.add(container.getProperty("InChIKey").toString());
//			line.add("");
//			line.add("");
//			line.add("");
//			line.add(container.getProperty("1A2").toString());
//			line.add(container.getProperty("2A6").toString());
//			line.add(container.getProperty("2B6").toString());
//			line.add(container.getProperty("2C8").toString());
//			line.add(container.getProperty("2C9").toString());
//			line.add(container.getProperty("2C19").toString());
//			line.add(container.getProperty("2D6").toString());
//			line.add(container.getProperty("2E1").toString());
//			line.add(container.getProperty("3A4").toString());
//			line.add(container.getProperty(CDKConstants.TITLE).toString());
//			
//			
//			try{
//				String isomericSmiles = sg.create(container);
//				IAtomContainer prepContainer = bt.preprocessContainer(container);
//				line.add(isomericSmiles);
//				String extendedFeatures = StringUtils.join(fgen.generateExtendedMolecularFeatures(prepContainer).split(","), "\t");
//				
//				ArrayList<Double> bioTransformerFingerprint_bits = cs.generateClassyfireFingerprintAsDouble(prepContainer, fpatterns).getBitValues();
//				for(int x = 0; x < bioTransformerFingerprint_bits.size(); x++){
//					extendedFeatures =  extendedFeatures + "\t" + String.valueOf(bioTransformerFingerprint_bits.get(x));
//					
//				}
//				
//				// Adding PubchemFingerprints
//				
//				
//				ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
//				IBitFingerprint fingerp		= pbf.getBitFingerprint(prepContainer);
//	
//				int[] onbits = fingerp.getSetbits();
//	
//				for(int kp = 0; kp < 881; kp++){
//					fingerprint_bits.add(0.0);
//				}
//				for(int o = 0; o < onbits.length; o++){
//					fingerprint_bits.set(onbits[o], 1.0);
//				}
//				extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(fingerprint_bits,"\t");
//				
//				// Adding MACCS Fingerprints
//				
//				ArrayList<Double> maccs_fingerprint_bits = new ArrayList<Double>();
//				IBitFingerprint maccs_fingerp		= maccs.getBitFingerprint(prepContainer);
//				
//				int[] maccs_onbits = maccs_fingerp.getSetbits();
//				
//				for(int kp = 0; kp < 166; kp++){
//					maccs_fingerprint_bits.add(0.0);
//				}
//				for(int o = 0; o < maccs_onbits.length; o++){
//					maccs_fingerprint_bits.set(maccs_onbits[o], 1.0);
//				}			
//				
//				extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(maccs_fingerprint_bits,"\t");
//				
//				line.add(extendedFeatures);
//				
//				String joinedLine = "";
//				boolean first = true;
//				final StringBuilder sb = new StringBuilder();
//				
//				for(String s : line){
//				
//					if(!first){
//						sb.append("\t");
//					}
//					else{
//						first = false;
//					}
//					sb.append(s);
//				}
//				System.out.println(sb.toString());
//				bw.write(sb.toString());
//			}
//			catch (Exception exp ){
//				System.err.println("Error here: Could not process compound with cinchikey " + container.getProperty("InChIKey").toString());
//				System.err.println(exp.getMessage());
//				
//			}
//			bw.newLine();
//		}
//		
		 ErdbTask e = new  ErdbTask();	 
//		 e.generateReactantsSetFromTSV();
//		 e.generateBioTransformerDB();
		
		 String[] cyps = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"};
		 IAtomContainerSet containers = FileUtilities.parseSdf("/Users/yandj/Projects/Metabolism/data_set/CYP_DBs/merged.sdf");
		 InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
//		 SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(sdfFileName) + 
//					"/" + FilenameUtils.getBaseName(sdfFileName)  + "_part_"+ part_nr + ".sdf" ));
		 SDFWriter sdfWriter = new SDFWriter(new FileOutputStream("/Users/yandj/Projects/Metabolism/data_set/CYP_DBs/biotransformer_db_may92018.sdf"));
		 
		 for(IAtomContainer atc : containers.atomContainers()){
			LinkedHashMap<Object, Object> lm = new LinkedHashMap<Object, Object>();
			String inchikey = inchiGenFactory.getInChIGenerator(atc).getInchiKey();
			LinkedHashMap<String,ArrayList<String>> c = ChemdbRest.getSynonymsObjectViaInChIKey(inchikey);
//			System.out.println(c);
			lm.put("InChIKey", inchikey);
			if(c!= null){
				lm.put("PubChem CID", c.get("CID").get(0));
			}
			else{
				lm.put("PubChem CID", "");
			}
			lm.put(CDKConstants.TITLE, (String) atc.getProperty(CDKConstants.TITLE));
			lm.put("SOM_1A2", "");
			lm.put("SOM_2A6", "");
			lm.put("SOM_2B6", "");
			lm.put("SOM_2C8", "");
			lm.put("SOM_2C9", "");
			lm.put("SOM_2C19", "");
			lm.put("SOM_2D6", "");
			lm.put("SOM_2E1", "");
			lm.put("SOM_3A4", "");
			
			Set set = new HashSet(Arrays.asList(((String) atc.getProperty("Citation")).split("\n")));
			lm.put("References",StringUtils.join(set,"\n"));
			lm.put("FLAG_COMMENTS", null);
			
//			for(int i = 0; i<cyps.length; i++){
//				String cyp = cyps[i];
//				ArrayList<String> l = new ArrayList<String>();
//				if(atc.getProperty("PRIMARY_SOM_"+cyp) != null){
//					String[] a = atc.getProperty("PRIMARY_SOM_"+cyp).toString().trim().split(" ");
//					for(int k = 0; k < a.length; k++){
//						if( !l.contains(a[k].trim())){
//							l.add(a[k].trim());
//						}
//						
//					}
//				}
//				if(atc.getProperty("SECONDARY_SOM_"+cyp) != null){
//					String[] a = atc.getProperty("SECONDARY_SOM_"+cyp).toString().trim().split(" ");
//					for(int k = 0; k < a.length; k++){
//						if( !l.contains(a[k].trim())){
//							l.add(a[k].trim());
//						}
//					}
//				}				
//				if(atc.getProperty("TERTIARY_SOM_"+cyp) != null){
//					String[] a = atc.getProperty("SECONDARY_SOM_"+cyp).toString().trim().split(" ");
//					for(int k = 0; k < a.length; k++){
//						if( !l.contains(a[k].trim())){
//							l.add(a[k].trim());
//						}
//					}
//				}				
////				if(atc.getProperty("SECONDARY_SOM_"+cyp) != null){
////					l.addAll( Arrays.asList(atc.getProperty("SECONDARY_SOM_"+cyp).toString().trim().split("\t")));
////				}
////				if(atc.getProperty("TERTIARY_SOM_"+cyp) != null){
////					l.addAll( Arrays.asList(atc.getProperty("TERTIARY_SOM_"+cyp).toString().trim().split("\t")));
////				}
//			
//				lm.put(("SOM_"+ cyp), StringUtils.join(l,"\n"));
//			
//			}
			atc.setProperties(lm);
			sdfWriter.write(atc);			 
		 }
		 
		sdfWriter.close();
		
	}
	
}
