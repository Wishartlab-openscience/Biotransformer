/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */


package biotransformer.transformation;

import biotransformer.biosystems.BioSystem;


/**
 * This class implements metabolic pathways. A metabolic pathway is created as a list of enzymes. Based on that list,
 * a limited number of reactions they catalyze will be applicable, and so, the number of possible substrates.
 */

import biotransformer.biosystems.BioSystem.BioSystemName;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import biotransformer.biomolecule.Enzyme;


public class MetabolicPathway {
	
	private String  name;
	private ArrayList<String> enzymeList;
	private LinkedHashMap<String,Enzyme> enzymes;
	private BioSystem bSystem;
	
	
	public enum MPathwayName {
			BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS,BIOSYNTHESIS_OF_AMINO_ACIDS, 
			PRIMARY_BILE_ACID_BIOSYNTHESIS, SECONDARY_BILE_ACID_BIOSYNTHESIS,
//			FATTY_ACID_BIOSYNTHESIS,
			FATTY_ACID_DEGRADATION, FATTY_ACID_ELONGATION, ETHER_LIPID_METABOLISM,
			GLYCEROLIPID_METABOLISM, GLYCEROPHOSPHOLIPID_METABOLISM, SPHINGOLIPID_METABOLISM,			
			UNSATURATED_FATTY_ACID_BIOSYNTHESIS, INOSITOL_PHOSPHATE_METABOLISM,
			MICROBIAL_POLYPHENOL_METABOLISM, BILE_ACID_METABOLISM
			
	}

	
//	public MetabolicPathway(String name, BioSystemName bsys, ArrayList<Enzyme> enzymeArray) {
//		// TODO Auto-generated constructor stub	
//		this.name = name;
//		this.bSystem = new BioSystem(bsys);
//		this.enzymeList = enzymeArray;	
//	}
	
	public MetabolicPathway(MPathwayName pathwayName, BioSystem bsys) {
		// TODO Auto-generated constructor stub	
	}	

	public MetabolicPathway(String pathwayName, BioSystem bsys) {
		// TODO Auto-generated constructor stub	
	}
	
	public void filterEnzymesAndReactions(){
	//	takes the pathway's list of enzymes from the database, and remove the ones that do not occur in the current biosystem,
		// so that this pathway reflects the selected organis only.
	}
	
}
