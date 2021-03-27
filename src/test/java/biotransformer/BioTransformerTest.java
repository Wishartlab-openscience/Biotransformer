package biotransformer;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.io.FilenameUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.utils.ChemStructureManipulator;


public class BioTransformerTest {

	public BioTransformerTest() {
		// TODO Auto-generated constructor stub
	}


	public static void main(String[] args) throws Exception {		
		Biotransformer humanBT = new Biotransformer(BioSystemName.HUMAN);
		Biotransformer bt = new Biotransformer(BioSystemName.GUTMICRO);

//		ArrayList<MetabolicReaction> customArray = new ArrayList<MetabolicReaction>(bt.bSystem.getReactionsHash().values()) ;
//		System.out.println("customArray: " + customArray.size());
		IAtomContainer ac = bt.getSmiParser().parseSmiles("[H]OC1=C(O)C2=C(OC(=C(O[H])C2=O)C2=C([H])C(O[H])=C(O)C(O)=C2[H])C([H])=C1[H]");
		
//		ReactantPred rp = new ReactantPred();
//		IAtomContainerSet inputMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		inputMolecules.addAtomContainer(ac);
//		ArrayList<HashMap<String,String>> predictedResult = new ArrayList<HashMap<String,String>>();
//		predictedResult = rp.initPreResults(predictedResult,inputMolecules.getAtomContainerCount());
//		ArrayList<HashMap<String,String>>  res = rp.makePrediction("CYP1A2", System.getProperty("user.dir")+"/supportfiles/CYP1A2/model/1A2_RT.model", inputMolecules, "supportfiles/CYP1A2/supportfile.csv", predictedResult);
//		System.out.println(res);
//		System.out.println(bt.smiGen.create(ac));
//		IAtomContainer stac = bt.standardizeMolecule(ac, true);
//		System.out.println(bt.smiGen.create(stac));
//		MetabolicReaction mreact = new MetabolicReaction(ReactionName.FLAVANON_3_OL_C_RING_FISSION);
//		
//		ArrayList<Biotransformation> acBiotransormations = bt.applyReactionAndReturnBiotransformations(stac, 
//				mreact, false);
//		
//		IAtomContainer ac2 = bt.smiParser.parseSmiles("Clc1c(OP(=S)(OCC)OCC)nc(Cl)c(Cl)c1");
//		ArrayList<Biotransformation> acBiotransormations = bt.applyReactionsChainAndReturnBiotransformations(ac2, customArray, true, true, 2);
//		
//		for(Biotransformation b : acBiotransormations){
//			b.display();
//			System.out.println("\n");
//		}
//		System.out.println("A total of " + acBiotransormations.size() + " transformations.");
//		IAtomContainerSet acMetabolites = bt.extractProductsFromBiotransformations(acBiotransormations);
//		System.out.println("A total of " + acMetabolites.getAtomContainerCount() + " metabolites.");
		
//		IAtomContainer a = bt.smiParser.parseSmiles("[H]OC1=C(O[H])C([H])=C2C(=O)OC3=C(C2=C1)C(=C([H])C(O[H])=C3O[H])C(O)=O");
//		IAtomContainer stac = bt.standardizeMolecule(a, true);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stac);
//		System.out.println(bt.smiGen.create(stac));
//		SMIRKSReaction se = bt.smrkMan.parse("[H][#8]-[#6;R0](=[O;R0])-[#6;R1:1]=,:1[#6;R1:6]=,:[#6;R1:5][#6;R1:4]=,:[#6;R2:3][#6;R2:2]=,:1>>[H][#6;R1:1]=,:1[#6;R1:6]=,:[#6;R1:5][#6;R1:4]=,:[#6;R2:3][#6;R2:2]=,:1");
//		MetabolicReaction mreact = new  MetabolicReaction(MRPatterns.ReactionName._4P_O_DEMETHYLATION_OF_FLAVONE);
		

//		ArrayList<Biotransformation> acBiotransormations = 	bt.applyReactionsChainAndReturnBiotransformations(stac, bt.reactionsList.get("gutMicroReactions"), true, true, 2, 0.0);
//		System.out.println(acBiotransormations.isEmpty());
//		
//			for(Biotransformation b : acBiotransormations){
//			b.display();
//			System.out.println("\n");
//		}
//					
//		IAtomContainerSet meta = bt.generateAllMetabolitesFromAtomContainer(stac, se, true);			
//		if(meta != null){			
//			System.out.println(meta.getAtomContainerCount());
//			for(IAtomContainer m : meta.atomContainers()){
//				System.out.println(bt.smiGen.create(m));
//			}
//		}
		
//		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//		String smirks = "[H][#6;A:1][#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16]):2]([#6;A:3])[#6;A:4]>>[#6;A:1]=O.[H][#7;A;X3:2]([#6;A:3])[#6;A:4]";
//		String smarts = "[H][#6;A][#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6;A])[#6;A]";
//		String rName = "N_DEALKYLATION_OF_ALIPHATIC_TERTIARY_AMINES_PATTERN1";
//		MetabolicReaction mm = new MetabolicReaction(ReactionName.valueOf(rName), smirks, new ArrayList<String>(Arrays.asList(smarts.split("; "))), 
//				new ArrayList<String>(), smrkMan);
//
//		
//		IAtomContainer ac = bt.smiParser.parseSmiles("C(N(C([H])([H])[H])C(C(C(=C1C2=C(C(=C(C(=C2C(C(C3=C(C(=C(C(=C13)[H])[H])[H])[H])([H])[H])([H])[H])[H])[H])[H])[H])[H])([H])[H])([H])[H])([H])([H])[H]");
//
//		if(ChemStructureExplorer.compoundMatchesReactionConstraints(mm, ac)){
//			IAtomContainerSet meta = bt.generateAllMetabolitesFromAtomContainer(ac, mm, true);
//			
//			for(IAtomContainer m : meta.atomContainers()){
//				System.out.println(bt.smiGen.create(m));			
//			}
//		}

//		IAtomContainer ac = bt.smiParser.parseSmiles("C(C(C([H])([H])[H])(C1(O[H])C(C(C(C([H])([H])[H])=C(C1([H])[H])[H])([H])[H])([H])[H])[H])([H])([H])[H]");
		
//		IAtomContainer ac = bt.smiParser.parseSmiles("C=1(C(C(C(C=2C(=C(C(=C(C2OC([C@](C(N(C(C(C([H])([H])[H])([H])[H])([H])[H])[H])([H])[H])(O[H])[H])([H])[H])[H])[H])[H])[H])=O)([H])[H])([H])[H])C(=C(C(=C(C1[H])[H])[H])[H])[H]");
//		IAtomContainer a = ChemStructureManipulator.preprocessContainer(ac);
		
//		MetabolicReaction mn = new MetabolicReaction(ReactionName.HYDROXYLATION_OF_ALIPHATIC_TERTIARY_PENULTIMATE_CARBON);
//		if(ChemStructureExplorer.compoundMatchesReactionConstraints(mn, a)){
//			IAtomContainerSet meta = bt.generateAllMetabolitesFromAtomContainer(a, mn, false);
//			System.out.println(meta.getAtomContainerCount());
//			for(IAtomContainer m : meta.atomContainers()){
//				System.out.println(bt.smiGen.create(m));			
//			}
//		}
		
//		ReactantPred rp = new ReactantPred();
//		IAtomContainerSet inputMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		inputMolecules.addAtomContainer(a);
//		ArrayList<HashMap<String,String>> predictedResult = new ArrayList<HashMap<String,String>>();
//		predictedResult = rp.initPreResults(predictedResult,inputMolecules.getAtomContainerCount());
//		
//		
//		
////		ArrayList<HashMap<String,String>>  res = rp.makePrediction("CYP1A2", "supportfiles/CYP1A2/model/1A2_RT.model", inputMolecules, "supportfiles/CYP1A2/supportfile.csv", predictedResult);
//		ArrayList<HashMap<String,String>>  res = rp.makeMultiPrediction("CYP1A2,CYP2A6, CYP2D6,CYP3A4", inputMolecules, predictedResult);
//		
//		System.out.println(res);

//		System.out.println("Reactions with preference rules: " + bt.bSystem.mrFilter.btPriorityHash.keySet().size());
//		int ruleCount = 0;
//		
//		for(Entry<ReactionName, ReactionName[]> r : bt.bSystem.mrFilter.btPriorityHash.entrySet()) {
//			System.out.println(r.getValue().length);
//			ruleCount = ruleCount + r.getValue().length;
//		}
//		
//		System.out.println("Total number of precedence rules: " + ruleCount);
//		
//		
//		Statistics stats = new Statistics();
//		System.out.println(stats.generalStatistics);
//
//		MDLRXNV2000Reader  rxnReader = new MDLRXNV2000Reader(new FileReader("data/rxn/26274.rxn"));
//		System.out.println(rxnReader.getFormat().getFormatName());
//		MDLRXNWriter writer = new MDLRXNWriter(new FileWriter(new File("data/output.mol")));
////		System.out.println(rxnReader.);
		
		
//		System
		
//		buildSdfFromTSV("/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_mini.txt");
		buildSdfTSVViaFromRandomSelection("../compounds_mini.txt", 50, 28771);
	}
	

	public static void buildSdfFromTSV(String tsvFileName) throws Exception{
		
		BufferedReader bRead = new BufferedReader(new FileReader(tsvFileName));
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(tsvFileName) + 
				"/" + FilenameUtils.getBaseName(tsvFileName)  + ".sdf" ));
		
		int counter = 0;
		String line = null;
		
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		
//		System.out.println(counter);
		
		while((line = bRead.readLine()) !=null){
			counter++;
			System.out.println(counter);
			String[] sline = line.split("\t");
			
			if(! (line.contains("moldb_smiles") ||  sline.length<4 || sline[3].contentEquals("NULL") )){
				
				System.out.println(counter + " " + sline[0] + " " + sline[1]);
				IAtomContainer atc= smiParser.parseSmiles(sline[3].trim());
				AtomContainerManipulator.suppressHydrogens(atc);
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
				atc.setProperty(CDKConstants.TITLE, sline[2]);
//				atc.setProperty("Name", sline[2]);
				
//				if(gen.getInchiKey().trim().contentEquals(sline[3].trim())){
//					atc.setProperty("InChIKey", gen.getInchiKey());
					atc.setProperty("InChIKey", sline[4].trim().replace("InChIKey=",""));
//				}
//				else{
//					System.err.println("Issue of inchikey incompatibility with " + sline[0]);
//					System.err.println(gen.getInchiKey());
//					System.err.println(sline[0].trim());
//	//				break;
//				}
//				if(sline.length>=5){
//					atc.setProperty("Origin", sline[4]);
//				}
//				else{
//					atc.setProperty("Origin", null);
//				}
				atc.setProperty("FooDB ID", sline[1]);
//				atc.setProperty("SoMs", "");
//				atc.setProperty("References","");
				sdfWriter.write( ChemStructureManipulator.preprocessContainer(atc) );
			}
		}
		
		sdfWriter.close();
		bRead.close();
	}
	
	
	public static void buildSdfTSVViaFromRandomSelection(String tsvFileName, int size, int originalSampleSize) throws Exception{
		
		BufferedReader bRead = new BufferedReader(new FileReader(tsvFileName));
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(tsvFileName) + "/"  +  
				FilenameUtils.getBaseName(tsvFileName) + "_random_" + 
				Integer.valueOf(size).toString() + ".sdf" ));
		
		int counter = 0;
		String line = null;
		
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		
//		System.out.println(counter);
		
		
		
//		List<Integer> numbersList = new ArrayList();
//		Collections.addAll(numbersList, numbers); 
		
		ArrayList<String[]> lines = new ArrayList<String[]>();
		
	    
		while((line = bRead.readLine()) !=null){
			counter++;
			System.out.println(counter);
			String[] sline = line.split("\t");
			if(!( line.contains("moldb_smiles") ||  sline.length<4 || sline[3].contentEquals("NULL") )){
				lines.add(sline);
			}
		}
		int newMaxSize=originalSampleSize;
		if(lines.size()<originalSampleSize){
			newMaxSize = lines.size();
		}
		
		System.out.println("new max. size: " + newMaxSize);
			
		int[] numbers  = sampleRandomNumbersWithoutRepetition(0, newMaxSize, size);

		for(int i = 0; i < numbers.length; i++){
			System.out.println(numbers[i] + " " +  lines.get(numbers[i])[0] + " " + lines.get(numbers[i])[1]);
			IAtomContainer atc= smiParser.parseSmiles( lines.get(numbers[i])[3].trim());
			AtomContainerManipulator.suppressHydrogens(atc);
			InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
			atc.setProperty(CDKConstants.TITLE,  lines.get(numbers[i])[2]);

			atc.setProperty("InChIKey",  lines.get(numbers[i])[4].trim().replace("InChIKey=",""));

			atc.setProperty("FooDB ID",  lines.get(numbers[i])[1]);
//				atc.setProperty("SoMs", "");
//				atc.setProperty("References","");
			sdfWriter.write( ChemStructureManipulator.preprocessContainer(atc) );				
						
		}
	
	sdfWriter.close();
	bRead.close();
	}	
	
	
	public static int[] sampleRandomNumbersWithoutRepetition(int start, int end, int count) {
	    Random rng = new Random();

	    int[] result = new int[count];
//	    int cur = 0;
//	    int remaining = end - start;
//	    for (int i = start; i < end && count > 0; i++) {
//	        double probability = rng.nextDouble();
//	        if (probability < ((double) count) / (double) remaining) {
//	            count--;
//	            result[cur++] = i;
//	        }
//	        remaining--;
//	    }
	    
	   ArrayList<Integer> l = new ArrayList<Integer>();
//	   int c = 0;
	   
	   while (l.size()<count){
		   int n = rng.nextInt((end - start) + 1) + start;
		   if(!l.contains(n)){
			   System.out.println("Size of List: " + l.size());
//			   c = c++;
			   l.add(n);
			   result[l.size()-1] = n;
		   }
		   
	   }

	   for(int j = 0; j <result.length; j++){
		   System.out.println(result[j]);
	   }
	   
	   return result;
	}
	
	
}
