/**
 * This class implements the class of metabolic enzymes.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.biomolecule;

import java.util.ArrayList;

import org.openscience.cdk.silent.SilentChemObjectBuilder;

import ambit2.smarts.SMIRKSManager;
import biotransformer.transformation.MetabolicReaction;

public class Enzyme {

	private String						name;
	private String						acceptedName;
	private ArrayList<String>				uniProtIds; //specific to a BioSystem
	private String 							primaryUniProtId;
	private ArrayList<MetabolicReaction>	reactionSet;
	private ArrayList<String> 				cellularLocations;
	private String							description;

	public Enzyme(String enzName) {
		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
		name = enzName;
//		reactionSet = generateReactionObjects(name, smrkMan);
	}


	public Enzyme(String enzName, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {	
//		Enzyme(enzName, null ,uniprot_ids, cellularLocations, reactions);
		this.name = enzName;
		this.uniProtIds = uniprot_ids;
		this.primaryUniProtId = uniprot_ids.get(0);
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		
	}

	
	
	public Enzyme(String enzName, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions, String acceptedcName) {	
//		Enzyme(enzName, null ,uniprot_ids, cellularLocations, reactions);
		this.name = enzName;
		this.acceptedName = acceptedcName;
		this.uniProtIds = uniprot_ids;
		this.primaryUniProtId = uniprot_ids.get(0);
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;

	}	
	public Enzyme(String enzName, String description, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {
		this.name = enzName;
		
		if(!(this.uniProtIds == null || this.uniProtIds.isEmpty())){
			this.uniProtIds = uniprot_ids;
			this.primaryUniProtId = uniprot_ids.get(0);
		}
		
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
	}
	
	public Enzyme(String enzName, String description, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions, String systematicName) {
		this.name = enzName;
		this.acceptedName = acceptedName;
		
		if(!(this.uniProtIds == null || this.uniProtIds.isEmpty())){
			this.uniProtIds = uniprot_ids;
			this.primaryUniProtId = uniprot_ids.get(0);
		}
		
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
	}
	
	public Enzyme(String enzName, String description, String uniprot_id, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {
		this.name = enzName;
		this.uniProtIds = new ArrayList<String>();
		this.uniProtIds.add(uniprot_id);
		this.primaryUniProtId = uniprot_id;
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
		
	}

	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet, ArrayList<String>cellLocations) {
		name = enzName;
		this.reactionSet = reactionSet;
		this.cellularLocations = cellLocations;

	}
	
	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet, SMIRKSManager smrkMan) {
		name = enzName;
		this.reactionSet = reactionSet;
	}	
	
	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet) {
		name = enzName;
		this.reactionSet = reactionSet;
	}
	
	
//	public Enzyme(String enzName) {
//		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//		name = EnzymeName.valueOf(enzName);
////		reactionSet = generateReactionObjects(name, smrkMan);
//	}
//
//	public Enzyme(EnzymeName enzName, SMIRKSManager smrkMan) {
//		name = enzName;
////		reactionSet = generateReactionObjects(name, smrkMan);
//	}
//
//	public Enzyme(String enzName, SMIRKSManager smrkMan) {
//		name = EnzymeName.valueOf(enzName);
////		reactionSet = generateReactionObjects(name, smrkMan);
//	}



	
	public void addReaction(MetabolicReaction reaction) {
		reactionSet.add(reaction);
	}

	public String getName() {
		return name.toString();
	}

	public String getPrimaryUniprotId() {
		return primaryUniProtId;
	}

	public ArrayList<String> getUniprotIds() {
		return uniProtIds;
	}
	
	public ArrayList<MetabolicReaction> getReactionSet() {
		return reactionSet;
	}

	public void setprimaryUniprotId(String uniprotID) {
		primaryUniProtId = uniprotID;
	}
	
	public void setUniprotIds(ArrayList<String> uniprotIds) {
		uniProtIds =  uniprotIds;
	}


	
	
//	public enum EnzymeName {
//		CYP1A1, CYP1A2, CYP1B1, CYP2A6, CYP2D4, CYP2A13, CYP2B1, CYP2B6, CYP2C8, CYP2C9, CYP2C18,
//		CYP2C19, CYP2D6, CYP2E1, CYP2F1, CYP3A4, CYP3A5, CYP3A7, CYP4A11, CYP4B1,
//		ABKAR1, ABKAR2,
//		FMO1, FMO2, FMO3,
//		EC_1_1_1_1, EC_1_1_1_35, EC_1_1_1_211, EC_1_1_1_213, EC_1_1_1_53, EC_1_1_1_90, EC_1_1_1_100, EC_1_1_1_145, EC_1_1_1_146,
//		EC_1_1_1_147, EC_1_1_1_149, EC_1_1_1_150, EC_1_1_1_151, EC_1_1_1_159,
//		EC_1_1_1_170, EC_1_1_1_184, EC_1_1_1_239, EC_1_1_1_330, EC_1_1_1_345, EC_1_1_1_357,
//		EC_1_1_1_375, EC_1_1_3_7, EC_1_1_3_15, EC_1_1_3_20, EC_1_2_1_7, EC_1_2_1_28,
//		EC_1_2_1_29, EC_1_2_1_30, EC_1_2_1_48, EC_1_2_1_50,
//		EC_1_2_3_1, EC_1_2_3_9, EC_1_2_99_7, EC_1_3_1_22, EC_1_3_1_34, EC_1_3_1_38, EC_1_3_1_70, EC_1_3_1_77,
//		EC_1_3_1_93, EC_1_3_8_1, EC_1_3_8_7, EC_1_3_3_6, EC_1_3_8_8, EC_1_3_8_9, EC_1_3_99_4,
//		EC_1_3_99_5, EC_1_3_99_6, EC_1_5_3_20, EC_1_6_5_2, EC_1_14_13_70, EC_1_14_13_142, 
//		EC_1_14_14_14, EC_1_14_14_16,EC_1_14_14_19, EC_1_14_14_29, 
//		EC_1_14_15_4, EC_1_14_19_5, EC_1_14_19_17, EC_1_14_19_18, EC_1_14_19_19,EC_1_14_14_32, EC_1_14_19_20,
//		EC_1_14_19_29, EC_2_1_1_6, EC_2_1_1_9, EC_2_1_1_15, EC_2_1_1_17, EC_2_1_1_42, 
//		EC_2_1_1_45, EC_2_1_1_67, EC_2_1_1_71, EC_2_1_1_96, EC_2_1_1_97, EC_2_1_1_317, 
//		EC_2_3_1_5, EC_2_3_1_15, EC_2_3_1_16, EC_2_3_1_20, EC_2_3_1_22, EC_2_3_1_23, EC_2_3_1_24,
//		EC_2_3_1_25, EC_2_3_1_26, EC_2_3_1_42, EC_2_3_1_43, EC_2_3_1_51, EC_2_3_1_52,
//		EC_2_3_1_62, EC_2_3_1_63, EC_2_3_1_65, EC_2_3_1_67, EC_2_3_1_75, EC_2_3_1_76, 
//		EC_2_3_1_80, EC_2_3_1_85, EC_2_3_1_87, EC_2_3_1_88, EC_2_3_1_102, EC_2_3_1_116, 
//		EC_2_3_1_118, EC_2_3_1_121, EC_2_3_1_125, EC_2_3_1_135,  EC_2_3_1_152, EC_2_3_1_153, EC_2_3_1_171, EC_2_3_1_199, 
//		EC_2_3_1_215, EC_2_4_1_17, EC_2_4_1_46, EC_2_4_1_47, EC_2_4_1_62, EC_2_4_1_80, EC_2_4_1_83, EC_2_4_1_91, EC_2_4_1_92, EC_2_4_1_115, EC_2_4_1_117, EC_2_4_1_157, 
//		EC_2_4_1_159, EC_2_4_1_173, EC_2_4_1_198, EC_2_4_1_239, EC_2_4_1_240, EC_2_4_1_241, EC_2_4_1_274, EC_2_4_99_9,
//		EC_2_5_1_26,
//		EC_2_6_1_57, EC_2_7_1_60, EC_2_7_1_67, EC_2_7_1_68, EC_2_7_1_94, EC_2_7_1_107,
//		EC_2_7_1_137, EC_2_7_1_138, EC_2_7_1_149, EC_2_7_1_150, EC_2_7_1_153,
//		EC_2_7_1_154, EC_2_7_7_41, EC_2_7_7_43, EC_2_7_8_1, EC_2_7_8_2, EC_2_7_8_8, EC_2_7_8_11, EC_2_7_8_27, EC_2_7_8_29,
//		EC_2_7_8_41, EC_2_8_2_1, EC_2_8_2_2, EC_2_8_2_11, EC_2_8_2_15, EC_3_1_1_1, EC_3_1_1_2,
//		EC_3_1_1_3, EC_3_1_1_4, EC_3_1_1_5, EC_3_1_1_8, EC_3_1_1_13, EC_3_1_1_23,
//		EC_3_1_1_32, EC_3_1_1_34, EC_3_1_1_47, EC_3_1_2_2, EC_3_1_2_6, EC_3_1_3_1, EC_3_1_3_2, EC_3_1_3_4,
//		EC_3_1_3_5, EC_3_1_3_6, EC_3_1_3_27, EC_3_1_3_29, EC_3_1_3_36, EC_3_1_3_64,
//		EC_3_1_3_66, EC_3_1_3_67, EC_3_1_3_78, EC_3_1_3_81, EC_3_1_3_86, EC_3_1_3_95,
//		EC_3_1_4_4, EC_3_1_4_11, EC_3_1_4_12, EC_3_1_4_16, EC_3_1_4_17, EC_3_1_4_39, EC_3_1_4_41,
//		EC_3_1_4_42, EC_3_1_4_43, EC_3_1_4_44, EC_3_1_4_50, EC_3_1_4_54, EC_3_1_6_1,
//		EC_3_1_6_2, EC_3_1_6_8, EC_3_1_8_1, EC_3_2_1_22,  EC_3_2_1_23, EC_3_2_1_45, EC_3_2_1_46, 
//		EC_3_2_1_62, EC_3_2_1_147, EC_3_2_2_8, EC_3_3_2_2, EC_3_3_2_5, EC_3_1_4_37, EC_3_1_4_46,
//		EC_3_3_3_2, EC_3_3_3_5, EC_3_5_1_14, EC_3_5_1_15, EC_3_5_1_17, EC_3_5_1_23,
//		EC_3_5_1_55, EC_3_5_1_60, EC_3_5_1_89, EC_3_5_1_114, EC_3_6_1_6, EC_3_6_1_7,
//		EC_3_6_1_9, EC_3_6_1_15, EC_3_6_1_19, EC_3_6_1_20, EC_3_7_1_5, EC_4_1_1_65,
//		EC_4_2_1_17, EC_4_2_1_134, EC_4_3_2_5, EC_5_3_3_1, EC_6_2_1_3, EC_6_2_1_20,  EC_3_7_1_,
//		EC_3_3_2_, EC_5_3_2_,
//		EC_1_6_99_1,
//		EC_2_3_1_13,
//		EC_1_14_18_5,
//		EC_2_7_8_5,
//		EC_3_1_1_26, 
//
//		
//		// X means, there is no number associated at the 4th level
//		EC_3_1_1_X,
//		
//		EC_3_2_1_20, EC_3_2_1_21, EC_3_2_1_24, EC_3_2_1_31, EC_3_2_1_50, EC_3_2_1_52,
//		
//		EC_3_2_1_X, EC_3_2_1_40,
//		EC_2_3_1_18,
//		
//		EC_3_5_2_X,
//		
//		EC_4_99_1_M4,
//		
//		
//		 // EC 3.2.1.22 , EC 2.3.1.158, EC 3.1.1.34
//		
//		// EC_3_1_4_40, 
//		/*
//		 * pseudo enzymes for gut microbial and environmental microbial this is more because of the rules are not 
//		 * always associated with a specific enzymes. The more we know about the enzyme(s) associated to a reaction, 
//		 * we can add that enzyme (if not available and then associate it with the given rule.
//		 */
//		UNSPECIFIED_GUT_BACTERIAL_ENZYME, UNSPECIFIED_ENVIRONMENTAL_BACTERIAL_ENZYME, DECARBOXYLASE, DEMETHYLASE, DEHYDROXYLASE, BACTERIAL_LACTONASE,
//		HYDROXYCINNAMATE_DECARBOXYLASE, VINYL_PHENOL_REDUCTASE, 
//		UDP_GLUCURONOSYLTRANSFERASE,
//		SULFOTRANSFERASE,
//		ACETYLTRANSFERASE,
//		BACTERIAL_BILE_SALT_3_HYDROXYSTEROID_DEHYDROGENASE,
//		BACTERIAL_BILE_SALT_7_HYDROXYSTEROID_DEHYDROGENASE,
//		BACTERIAL_BILE_SALT_12_HYDROXYSTEROID_DEHYDROGENASE,
//		EC_3_5_1_24,
//		UNSPECIFIED_BACTERIAL_ISOFLAVONE_REDUCTASE,
//		
//		// Add reactions for these enzymes
//		EC_2_4_1_37, EC_2_4_1_65, EC_2_4_1_79, EC_2_4_1_149, EC_2_4_1_228, 
//		EC_1_1_1_101, EC_3_2_1_18, EC_1_14_13_205,
//		EC_3_1_1_64, EC_3_1_1_20,
//		EC_2_5_1_18,
//		EC_3_1_1_81,
//		
//		EC_1_1_1_102,
//		EC_2_4_1_38, EC_2_4_1_40, EC_2_4_1_69, EC_2_4_1_86, EC_2_4_1_87, EC_2_4_1_88, EC_2_4_1_90, EC_2_4_1_152, EC_2_4_1_206, 
//		EC_2_4_1_275,
//		EC_2_4_99_1, EC_2_4_99_2, EC_2_4_99_4, EC_2_4_99_6, EC_2_4_99_8,
//		// https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-5921&orgids=ECOLI
//		// catalyzing carnitine ad choline to TMA
//		EC_1_14_13_M4,
//		
//		EC_3_7_1_4,
//		
//		FLAVONOID_C_GLYCOSIDASE,
//		BACTERIAL_NITROREDUCTASE,
//		EC_1_6_2_4,
//		CYPB5_CYPBR5,
//		EC_2_7_8_X,
//		EC_3_1_1_28,
//		EC_3_5_1_13,
//		EC_1_7_1_6,
//		EC_1_8_99_3,
//		EC_4_2_1_84,
//		EC_3_5_5_1,
//		EC_3_5_5_7,
//		EC_3_1_1_73,
//		EC_3_1_1_45,
//		EC_3_8_1_2,
//		EC_3_8_1_5,
//		EC_3_5_2_2,
//		EC_6_2_1_7,
//		UNSPECIFIED_BACTERIAL_SULFOXIDE_REDUCTASE,
//		PHENOLIC_ACID_DECARBOXYLASE,
//		EC_4_1_1_28,
//		EC_1_14_15_15,
//		EC_1_3_1_31,
//		EC_4_3_99_4,
//		EC_3_1_6_19
//	}
	

	public ArrayList<String> getReactionsNames() {
		ArrayList<String> names = new ArrayList<String>();
		for (int i = 0; i < reactionSet.size(); i++) {
			names.add(reactionSet.get(i).getReactionName());
		}
		return names;
	}
	
	@Override
	public String toString(){
		String representation = "";
		representation += String.format("%-20s\t%-25s\n","Name:",this.name);
		representation += String.format("%-30s\t%-35s\n","Systematic name:",this.acceptedName);
		representation += String.format("%-20s\n","Reactions:");
		
		for(MetabolicReaction m : this.getReactionSet()) {
			representation += String.format("%-60s\n",m.name);
		}
		
		return representation;
	}

}
