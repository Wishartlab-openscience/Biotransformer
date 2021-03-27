package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesGenerator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Biotransformer.bType;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;

public class BiotransformerSequence {
	//	public enum sequenceMode {
	//	
	//}
	
	protected ArrayList<BiotransformerSequenceStep> sequence;
	protected double scoreThreshold;
	protected Biotransformer utilityBiotransformer = null;
	

	public BiotransformerSequence(ArrayList<BiotransformerSequenceStep> mySequence, 
			double scoreThreshold) {
		this.sequence 		= mySequence;
		this.scoreThreshold = scoreThreshold;
	}	
	public BiotransformerSequence(ArrayList<BiotransformerSequenceStep> mySequence) {
		this.sequence = mySequence;
		this.scoreThreshold = 0.0;
	}

	
	public BiotransformerSequence(String mySequence, double scoreThreshold) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence, scoreThreshold);
		this.scoreThreshold = scoreThreshold;
	}
	public BiotransformerSequence(String mySequence) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence, 0.0);
		this.scoreThreshold = 0.0;
	}	
	
	
	
	private ArrayList<BiotransformerSequenceStep> createSequenceFromString(String mySequence, double scoreThreshold) throws Exception{
		ArrayList<BiotransformerSequenceStep> final_seq = new ArrayList<BiotransformerSequenceStep>();
		try {
			ArrayList<BiotransformerSequenceStep> seq = new ArrayList<BiotransformerSequenceStep>();
			String[] steps = mySequence.split(";");
			for(String step : steps) {
//				System.out.println("step: " + step);
				String[] attr = step.split(":");
				if(attr.length==2) {
					String btype_str = attr[0].trim().toUpperCase();
//					System.out.println('"'+attr[1].trim()+'"');
					int nsteps = Integer.valueOf(attr[1].trim());
					
					try {
						bType btype_ = Biotransformer.bType.valueOf(btype_str);
					}
					catch (IllegalArgumentException ie) {
						throw new IllegalArgumentException("Illegal biotransformer sequence type '" + btype_str + "'.");
					}
					catch (Exception e) {
						System.err.println(e.getMessage());
					}
					
//					if (btype_ != null) {
//						System.out.println("btype = " + btype_);
//					}
//					else {
//						throw new IllegalArgumentException("Illegal biotransformer sequence type '" + btype_str + "'.");
//					}
									
					seq.add(new BiotransformerSequenceStep(Biotransformer.bType.valueOf(btype_str), nsteps, scoreThreshold));
				}
				else {
					throw new IllegalArgumentException("Illegal biotransformer sequence '" + mySequence + "\n\tMake sure your sequence is semi-colon separated, and that each step is colon-separated pair of one valid biotransformer type, and one integer.");					
					
				}
			}	
			final_seq = seq;
		}
		catch (Exception e){
			System.out.println(e);
		}

		return final_seq;
	}
	
	public String toString() {
		String str_representation = "";
		ArrayList<String> steps = new ArrayList<String>();
		int counter = 0;
		for(BiotransformerSequenceStep s : this.sequence) {
			counter++;
			steps.add("Step " + counter + ": " + s.toString());
		}
		str_representation = String.join("\n", steps);
		return str_representation;
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainer startingCompound, double scoreThreshold) throws Exception{
		return runSequence(startingCompound, scoreThreshold, 1);
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainer startingCompound, double scoreThreshold, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		LinkedHashMap<Biotransformer.bType, Object> btransformers = new LinkedHashMap<Biotransformer.bType, Object>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentMetabolites.addAtomContainer(startingCompound);
		
		for(BiotransformerSequenceStep step : this.sequence) {
			ArrayList<Biotransformation> currentBiots =  new ArrayList<Biotransformation>();
			
			if(step.btype == bType.ALLHUMAN) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new HumanSuperBioTransformer());
				}
				currentBiots = ((HumanSuperBioTransformer) btransformers.get(step.btype)).predictAllHumanBiotransformationChain(currentMetabolites, step.nOfIterations, this.scoreThreshold, cyp450Mode);
			}
			
			if(step.btype == bType.SUPERBIO) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new HumanSuperBioTransformer());
				}
				for(IAtomContainer atc : currentMetabolites.atomContainers()) {
					currentBiots = ((HumanSuperBioTransformer) btransformers.get(step.btype)).simulateHumanSuperbioMetabolism(atc, this.scoreThreshold ,cyp450Mode);
				}				
			}		
			else if(step.btype == bType.CYP450) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new Cyp450BTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((Cyp450BTransformer) btransformers.get(step.btype)).predictCyp450BiotransformationChainByMode(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold, cyp450Mode);			
			}
			else if(step.btype == bType.ECBASED) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new ECBasedBTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((ECBasedBTransformer) btransformers.get(step.btype)).simulateECBasedMetabolismChain(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);		
			}
			else if(step.btype == bType.ENV) {
//				System.out.println("Starting ENV");
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new EnvMicroBTransformer());
				}
				currentBiots = ((EnvMicroBTransformer) btransformers.get(step.btype)).applyEnvMicrobialTransformationsChain(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);
//				System.out.println("Number of env. biotransformations: " + currentBiots.size());
			}			
			else if(step.btype == bType.HGUT) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new HGutBTransformer());
				}
				currentBiots = ((HGutBTransformer) btransformers.get(step.btype)).simulateGutMicrobialMetabolism(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);						
			}
			else if(step.btype == bType.PHASEII) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new Phase2BTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((Phase2BTransformer) btransformers.get(step.btype)).applyPhase2TransformationsChainAndReturnBiotransformations(currentMetabolites, true, true, true, step.nOfIterations, this.scoreThreshold);										
			}
			else if(step.btype == bType.SUPERBIO) {
				throw new IllegalArgumentException("Invalid Argument: SUPERBIO cannot be used within a sequence, as it is already a customized sequence. Valid biotransformers within a sequence are ALLHUMAN, CYP450, ECBASED, ENVMICRO, HGUT, PHASEII.");
			}
			
			currentMetabolites.add(((Biotransformer) btransformers.get(step.btype)).extractProductsFromBiotransformations(currentBiots));
			biotransformations.addAll(currentBiots);
		}
		
//		this.utilityBiotransformer = ((Collection<Entry<bType, Object>>) btransformers.entrySet()).stream().reduce((first, second) -> second).orElse(null)biotransformations.gtValue();
//		for(Biotransformation b: biotransformations) {
//			b.display();		
//		}
		

		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	

	public ArrayList<Biotransformation> runSequence(IAtomContainerSet startingCompounds, double scoreThreshold) throws Exception{
		return runSequence(startingCompounds, scoreThreshold, 1);
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainerSet startingCompounds, double scoreThreshold, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		SmilesGenerator smiGen 		= new SmilesGenerator().isomeric(); 
		for(IAtomContainer starting_ac : startingCompounds.atomContainers()) {
//			System.out.println(smiGen.create(starting_ac));
			try {
				biotransformations.addAll(runSequence(starting_ac, scoreThreshold, cyp450Mode));
			}
			catch (Exception e) {
				System.out.println(e);
			}
			
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
}
