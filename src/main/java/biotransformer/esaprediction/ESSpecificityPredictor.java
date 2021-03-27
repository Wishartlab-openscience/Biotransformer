/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.esaprediction;

import java.io.IOException;
import java.util.ArrayList;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import ambit2.smarts.query.SMARTSException;
import biotransformer.biomolecule.Enzyme;
//import biotransformer.biomolecule.Enzyme.EnzymeName;
import biotransformer.biosystems.BioSystem;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemicalClassFinder;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;
import reactantpredictor.BioTransformerAPIs;

/**
 * Predicts whether a compound is a substrate of a given enzyme within a specific biosystem.
 * @author Djoumbou Feunang, Yannick
 *
 */

public class ESSpecificityPredictor {
	public static final String[] cyp450EnzymeList = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"};
	BioSystem bSys;

	/**
	 * Creates a specific 
	 * @param bsys
	 */
	public ESSpecificityPredictor(BioSystem bsys) {
		// TODO Auto-generated constructor stub
		this.bSys = bsys;
	
	}
	
	public boolean isValidCyp450Substrate(IAtomContainer substrate, String enz) throws Exception{
	
		// ADD testing condition for scenario where an enzyme is not in the biosystem.	
		ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(substrate);
		return isValidCyp450Substrate(substrate, enz, chemClasses);
	}
	
	
	public boolean isValidCyp450Substrate(IAtomContainer substrate, String enz, ArrayList<ChemicalClassName> chemClasses) throws Exception{
		
		// ADD testing condition for scenario where an enzyme is not in the biosystem.	
		if(!(enz.toString().contains("CYP1A2") || enz.toString().contains("CYP2A6") || enz.toString().contains("CYP2B6")
				|| enz.toString().contains("CYP2C8") || enz.toString().contains("CYP2C9") || enz.toString().contains("CYP2C19")
				|| enz.toString().contains("CYP2D6") || enz.toString().contains("CYP2E1") || enz.toString().contains("CYP3A4"))){
			
			throw new IllegalArgumentException(enz.toString() + " is not a valid CYP isozyme for this system. The selected isozyme must"
					+ "be either of the following: CYP1A2, CYP2A6, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP2E1, or CYP3A4.");
		} else if(!ChemStructureExplorer.isMixture(substrate)){
			boolean validCyp450 = false;
//			System.err.println(chemClasses);
			
			if(!( ChemStructureExplorer.getMajorIsotopeMass(substrate) > 1500.0 || 
					chemClasses.contains(ChemicalClassName.GLYCOSYLATED_COMPOUND) ||
					chemClasses.contains(ChemicalClassName.GLUTATHIONE_CONJUGATE) ||
					chemClasses.contains(ChemicalClassName.SULFATED_COMPOUND) ||
					chemClasses.contains(ChemicalClassName.GLYCINATED_COMPOUND) ||
					chemClasses.contains(ChemicalClassName.TAURINATED_COMPOUND) ||
					chemClasses.contains(ChemicalClassName.GLUTAMATE_CONJUGATE) ||
					chemClasses.contains(ChemicalClassName.ACYLCARNITINE_CONJUGATE) ||
					chemClasses.contains(ChemicalClassName.CYSTEINYLGLYCINE_S_CONJUGATE) ||					
					chemClasses.contains(ChemicalClassName.ACYL_CoA_CONJUGATE) ||
					chemClasses.contains(ChemicalClassName.TETRAPYRROLE) ||
					chemClasses.contains(ChemicalClassName.SACCHARIDE) ||
					chemClasses.contains(ChemicalClassName.ETHER_LIPID) ||
					chemClasses.contains(ChemicalClassName.GLYCEROLIPID) ||
					chemClasses.contains(ChemicalClassName.GLYCEROPHOSPHOLIPID) ||
					chemClasses.contains(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL) ||
					chemClasses.contains(ChemicalClassName.SPHINGOLIPID) )) {
				
				if(this.bSys.name.toString() == "HUMAN"){
					BioTransformerAPIs reactantPred = new BioTransformerAPIs();
					validCyp450 = reactantPred.predictReactant(substrate, cyp450EnzymeList);
				}
			}
			// use CYP450 predictor
			return validCyp450;
		}
		else{
			return false;
		}
	}
	
	
	public boolean isValid_EC_2_8_2_2_substrate(IAtomContainer substrate, boolean preprocess) throws Exception{
		boolean validEC_2_8_2_2 = false;
		
		IAtomContainer atc = substrate.clone();
		
//		if(preprocess){
//			
//		}
		
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(substrate);	
		
		// phenolic hydroxysteroid
		Pattern smp = SmartsPattern.create(
				"[#6]=,:1[#6]=,:[#6][#6]~2=,:[#6]([#6]=,:1)~[#6,#7]~[#6,#7]~[#6,#7]~1-,=[#6,#7]-,=3-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A]-,=[#6,#7]-,=3~[#6,#7,#16]~[#6,#7,#16]~[#6,#7]~2~1"
				, SilentChemObjectBuilder.getInstance());
		
		validEC_2_8_2_2 = (!smp.matches(atc)) && (isPotentialSubstrateByReactionPatternMatching(atc, "EC_2_8_2_2"));
		
		return validEC_2_8_2_2;
	}
	
	
	public boolean isValidSubstrate(IAtomContainer substrate, String enz) throws Exception {
		
//		boolean validSubstrate = false;
//		
//		
//		if(  (enz.toString().contains("CYP1A2") || enz.toString().contains("CYP2A6") || enz.toString().contains("CYP2B6")
//				|| enz.toString().contains("CYP2C8") || enz.toString().contains("CYP2C9") || enz.toString().contains("CYP2C19")
//				|| enz.toString().contains("CYP2D6") || enz.toString().contains("CYP2E1") || enz.toString().contains("CYP3A4"))){
//			
//			
//			return isValidCyp450Substrate(substrate, enz);
//		} 
//		
//		else if(enz.toString().contentEquals("EC_2_8_2_2")){
//			return isValid_EC_2_8_2_2_substrate(substrate,false);
//		}
//		
//		else{
//			
//			/* This will apply to the enzyme for which we do not have machine-learning prediction models. Specific substrate specificity
//			 * prediction methods can be implemented for specific enzymes, if the rules take into consideration more than the pattern matching.
//			 * For instance, for some molecules, one could incorporate the LogP.
//			 */
////			if(this.bSys.getEnzymeHash().containsKey(enz)){
////				validSubstrate = isPotentialSubstrateByReactionPatternMatching(substrate, enz);
////			}
//			validSubstrate =  true;
//		}
//		return validSubstrate;
		ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(substrate);
		return isValidSubstrate(substrate, enz, chemClasses);
	}

	public boolean isValidSubstrate(IAtomContainer substrate, String enz, ArrayList<ChemicalClassName> chemClasses) throws Exception {
		
		boolean validSubstrate = false;
		
		
		if(  (enz.toString().contains("CYP1A2") || enz.toString().contains("CYP2A6") || enz.toString().contains("CYP2B6")
				|| enz.toString().contains("CYP2C8") || enz.toString().contains("CYP2C9") || enz.toString().contains("CYP2C19")
				|| enz.toString().contains("CYP2D6") || enz.toString().contains("CYP2E1") || enz.toString().contains("CYP3A4"))){
			
			
			return isValidCyp450Substrate(substrate, enz, chemClasses);
		} 
		
		else if(enz.toString().contentEquals("EC_2_8_2_2")){
			return isValid_EC_2_8_2_2_substrate(substrate,false);
		}
		
		else{
			
			/* This will apply to the enzyme for which we do not have machine-learning prediction models. Specific substrate specificity
			 * prediction methods can be implemented for specific enzymes, if the rules take into consideration more than the pattern matching.
			 * For instance, for some molecules, one could incorporate the LogP.
			 */
//			if(this.bSys.getEnzymeHash().containsKey(enz)){
//				validSubstrate = isPotentialSubstrateByReactionPatternMatching(substrate, enz);
//			}
			validSubstrate =  true;
		}
		return validSubstrate;		
	}
	
	
	
	
	public boolean isPotentialSubstrateByReactionPatternMatching(IAtomContainer substrate, String enz) throws SMARTSException, CDKException, IOException{
		boolean isPotentialSubstrate = false;
		
		if(this.bSys.getEnzymeHash().containsKey(enz)){
			Enzyme e = this.bSys.getEnzymeHash().get(enz);
//			System.out.println(e.getName());
//			System.out.println(e.getReactionSet().size());
			for( MetabolicReaction mreact : e.getReactionSet()){
//				System.out.println(mreact.name);
				if(ChemStructureExplorer.compoundMatchesReactionConstraints(mreact, substrate)){
					isPotentialSubstrate = true;
//					System.out.println("This compound is a substrate of " + enz + " and matches at least " + mreact.name);
					break;
				}
			}
			return isPotentialSubstrate;
		} else {
			throw new IllegalArgumentException(enz.toString() + " is not associated with the biosystem " + this.bSys.name);
		}
	}

	


}
