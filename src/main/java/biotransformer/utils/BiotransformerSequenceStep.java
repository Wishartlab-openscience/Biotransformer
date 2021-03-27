package biotransformer.utils;

import biotransformer.btransformers.Biotransformer;

public class BiotransformerSequenceStep {
	public Biotransformer.bType btype;
	public int nOfIterations;
	public double scoreThreshold;

	public BiotransformerSequenceStep(Biotransformer.bType btype, int nofIterations, double scoreThreshold) {
		this.btype 		= btype;
		this.nOfIterations = nofIterations;
		this.scoreThreshold = scoreThreshold;

	}	

	public BiotransformerSequenceStep(Biotransformer.bType btype) {
		this.btype 		= btype;
		this.nOfIterations = 1;
//		this.scoreThreshold = null;

	}
	
	public BiotransformerSequenceStep(Biotransformer.bType btype, double scoreThreshold) {
		this.btype 		= btype;
		this.nOfIterations = 1;
		this.scoreThreshold = scoreThreshold;

	}
	
	public BiotransformerSequenceStep(Biotransformer.bType btype, int nofIterations) {
		this.btype 		= btype;
		this.nOfIterations = nofIterations;
//		this.scoreThreshold = 0.0;

	}
	
}

