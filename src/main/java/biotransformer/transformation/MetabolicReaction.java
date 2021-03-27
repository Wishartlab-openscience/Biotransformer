package biotransformer.transformation;

/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

//import java.util.Arrays;
//import java.util.HashSet;
//import java.util.LinkedHashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//import java.util.stream.Collectors;
import java.util.ArrayList;



//import com.google.common.base.Joiner;

import ambit2.smarts.SMIRKSReaction;
import biotransformer.biomolecule.Enzyme;
import biotransformer.transformation.MRPatterns.ReactionName;
import ambit2.smarts.SMIRKSManager;


/**
 * The first reactantsSMARTS is the main one, for which indices might be
 * retrieved to apply position-specific transformations. The subsequent SMARTS
 * (excludedReactantsSMARTS) will be checked just to match further constraints.
 * For instance, long-chain fatty acids must NOT match the pattern of
 * very-long-chain fatty acids. This combination of SMARTS patterns could be
 * replaced in favor of MARKUSH formats. Because we are using CDK, we rather
 * stick up the the combination of SMARTS
 *
 */

public class MetabolicReaction {

	public String				name;
	public String				commonName;
	public String				reactionsBTMRID;
	private ArrayList<String>	reactantsSMARTS			= new ArrayList<String>();
	private ArrayList<String>	excludedReactantsSMARTS	= new ArrayList<String>();
	private String				productsSMARTS;
	private String				reactionSMIRKS;
	private String				reactionSMIRKS_text;
	private String				reactionEquation;
	private SMIRKSReaction		smirksReaction;
	private ArrayList<Enzyme>	catalyzingEnzymes 		= new ArrayList<Enzyme>();

	public static void main (String[] args){

	}
	
	
//	public MetabolicReaction(String r_name) {
//		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//		this.name = r_name.toString();
//		this.reactionSMIRKS = MRPatterns.setOfSMIRKS.get(name);
//
//		for (String smarts : MRPatterns.setOfReactantSMARTS.get(name)) {
//			if (!(smarts.trim().isEmpty())) {
//				this.reactantsSMARTS.add(smarts);
//			}
//		}
//
//		if (!(MRPatterns.setOfExcludedReactantSMARTS.get(r_name.toString()) == null)) {
//			for (String n_smarts : MRPatterns.setOfExcludedReactantSMARTS.get(name)) {
//				if (!(n_smarts.trim().isEmpty())) {
//					this.excludedReactantsSMARTS.add(n_smarts);
//				}
//			}
//		}
//		this.smirksReaction = smrkMan.parse(reactionSMIRKS);
//
//	}
//	
//	public MetabolicReaction(ReactionName r_name) {
//		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//		this.name = r_name.toString();
//		this.reactionSMIRKS = MRPatterns.setOfSMIRKS.get(name);
//
//		for (String smarts : MRPatterns.setOfReactantSMARTS.get(name)) {
//			if (!(smarts.trim().isEmpty())) {
//				this.reactantsSMARTS.add(smarts);
//			}
//		}
//
//		if (!(MRPatterns.setOfExcludedReactantSMARTS.get(r_name.toString()) == null)) {
//			for (String n_smarts : MRPatterns.setOfExcludedReactantSMARTS.get(name)) {
//				if (!(n_smarts.trim().isEmpty())) {
//					this.excludedReactantsSMARTS.add(n_smarts);
//				}
//			}
//		}
//		this.smirksReaction = smrkMan.parse(reactionSMIRKS);
//
//	}

	
	
	public MetabolicReaction(String r_name, String commmonName, String reactionsBTMRID, String reactionSMIRKS, ArrayList<String> reactantsSMARTS, ArrayList<String> excludedReactantsSMARTS, SMIRKSManager smrkMan) {
		this.name = r_name.toString();
		this.commonName = commmonName;
		this.reactionsBTMRID = reactionsBTMRID;
		this.reactionSMIRKS = reactionSMIRKS;
		this.reactantsSMARTS = reactantsSMARTS;
		this.excludedReactantsSMARTS = excludedReactantsSMARTS;
		this.smirksReaction = smrkMan.parse(reactionSMIRKS);
	}

	
	public MetabolicReaction(String r_name, String commmonName, String reactionSMIRKS, ArrayList<String> reactantsSMARTS, ArrayList<String> excludedReactantsSMARTS, SMIRKSManager smrkMan) {
		this.name = r_name.toString();
		this.commonName = commmonName;
		this.reactionSMIRKS = reactionSMIRKS;
		this.reactantsSMARTS = reactantsSMARTS;
		this.excludedReactantsSMARTS = excludedReactantsSMARTS;
		this.smirksReaction = smrkMan.parse(reactionSMIRKS);
	}
	
	public MetabolicReaction(ReactionName r_name, String reactionSMIRKS, ArrayList<String> reactantsSMARTS, ArrayList<String> excludedReactantsSMARTS, SMIRKSManager smrkMan) {
		this.name = r_name.toString();
		this.reactionSMIRKS = reactionSMIRKS;
		this.reactantsSMARTS = reactantsSMARTS;
		this.excludedReactantsSMARTS = excludedReactantsSMARTS;
		this.smirksReaction = smrkMan.parse(reactionSMIRKS);
	}
	
//	public MetabolicReaction(String smirks) {
//		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//		reactionSMIRKS = smirks;
//		smirksReaction = smrkMan.parse(reactionSMIRKS);
//	}

//	public MetabolicReaction(ReactionName r_name, SMIRKSManager smrkMan) {
//		name = r_name.toString();
//		reactionSMIRKS = MRPatterns.setOfSMIRKS.get(r_name.toString());
//
//		for (String smarts : MRPatterns.setOfReactantSMARTS.get(r_name.toString())) {
//			if (!(smarts.trim().isEmpty())) {
//				reactantsSMARTS.add(smarts);
//			}
//		}
//
//		if (!(MRPatterns.setOfExcludedReactantSMARTS.get(r_name.toString()) == null)) {
//			for (String n_smarts : MRPatterns.setOfExcludedReactantSMARTS.get(r_name.toString())) {
//				if (!(n_smarts.trim().isEmpty())) {
//					excludedReactantsSMARTS.add(n_smarts);
//				}
//			}
//		}
//
//		smirksReaction = smrkMan.parse(reactionSMIRKS);
//
//	}

	public MetabolicReaction(String smirks, SMIRKSManager smrkMan) {
		reactionSMIRKS = smirks;
		smirksReaction = smrkMan.parse(reactionSMIRKS);
	}

	public ArrayList<String> getReactantSMARTS() {
		return this.reactantsSMARTS;
	}

	public ArrayList<String> getExcludedReactantsSMARTS() {
		return this.excludedReactantsSMARTS;
	}	
	
	public String getBTRMID(){
		return this.reactionsBTMRID;
	}
	
	public String getComonName(){
		return this.commonName;
	}
	
	public String getReactionSMIRKS() {
		return this.reactionSMIRKS;
	}

	public SMIRKSReaction getSmirksReaction() {
		return this.smirksReaction;
	}

	public String getReactionName() {
		return this.name;
	}

	public String getReactionEquation() {
		return reactionEquation;
	}

	public void setReactionEquation(String rEquation) {
		reactionEquation = rEquation;
	}

	public String transformationDataToString() {
//		String description = "Name: " + name + "\n"
//				+ smirksReaction.transformationDataToString();
		String description = "";
		description.concat(String.format("%-20s\t%-25s","Name:",this.name));
		description.concat(String.format("%-20s\t%-150s","SMIRKS:",this.reactionSMIRKS));
//		description.concat(String.format("%-20s\t","Catalyzing enzymes:"));		
//		for(Enzyme enz : this.catalyzingEnzymes){
//			System.out.print(String.format("%-7s, ",enz.getName()));
//		}		
		
		System.out.print("\n");
		return description;
	}
	
	@Override
	public String toString(){
		String representation = "";
		representation += String.format("%-20s\t%-25s\n","Name:", this.name);
		representation += String.format("%-30s\t%-35s\n","Common name:", this.commonName);
		representation += String.format("%-20s\t%-300s\n","SMIRKS:", this.reactionSMIRKS);
		representation += String.format("%-20s", "reactantsSMARTS:") + this.reactantsSMARTS + "\n";
		representation += String.format("%-30s", "excludedReactantsSMARTS:") + this.excludedReactantsSMARTS + "\n";
		if(this.reactionsBTMRID !=null){
			representation += String.format("%-20s\t%-25s\n","btmrID:", this.reactionsBTMRID);
		}
//		System.out.print(String.format("%-20s\t","Catalyzing enzymes:"));		
//		for(Enzyme enz : this.catalyzingEnzymes){
//			System.out.print(String.format("%-7s, ",enz.getName()));
//		}		
		
		return representation;
	}
	
	
}
