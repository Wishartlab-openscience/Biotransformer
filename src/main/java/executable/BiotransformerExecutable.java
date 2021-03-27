package executable;

import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Biotransformer.bType;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.BiotransformerSequence;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HumanSuperBioTransformer;
import biotransformer.utils.MetaboliteFinder;
import biotransformer.utils.MetaboliteFinder.FinderOption;
import biotransformer.utils.UniversalBioTransformer;
import biotransformer.version.Version;



public class BiotransformerExecutable {
	private static LinkedHashMap<String, Biotransformer.bType> optionsToBtTypes = new LinkedHashMap<String, Biotransformer.bType>();

	public BiotransformerExecutable() {
		// TODO Auto-generated constructor stub
		setOptionsTobType();
		
	}
	
	
	private static void setOptionsTobType(){
		optionsToBtTypes.put("allhuman", Biotransformer.bType.ALLHUMAN);
		optionsToBtTypes.put("cyp450", Biotransformer.bType.CYP450);
		optionsToBtTypes.put("ecbased", Biotransformer.bType.ECBASED);		
		optionsToBtTypes.put("env", Biotransformer.bType.ENV);
		optionsToBtTypes.put("envmicro", Biotransformer.bType.ENV);
		optionsToBtTypes.put("hgut", Biotransformer.bType.HGUT);	
		optionsToBtTypes.put("phaseii", Biotransformer.bType.PHASEII);
		optionsToBtTypes.put("phase2", Biotransformer.bType.PHASEII);
		optionsToBtTypes.put("superbio", Biotransformer.bType.SUPERBIO);
	}
	private static Options generateOptions(){
		
		final Option taskOption = Option.builder("k")
				.required(true)
				.hasArg(true)
				.argName("BioTransformer Task")
				.longOpt("task")
				.desc("The task to be permed: pred for prediction, or cid for compound identification ")
				.build();
		
		final Option biotransformerOption = Option.builder("b")
				.required(false)
				.hasArg(true)
				.argName("BioTransformer Option")
				.longOpt("btType")
				.desc("The type of description: Type of biotransformer - EC-based  (ecbased), CYP450 (cyp450), Phase II (phaseII), "
						+ "Human gut microbial (hgut), human super transformer* (superbio, or allHuman), Environmental microbial (envimicro)**.\n"
						+ "If option -m is enabled, the only valid biotransformer types are allHuman, superbio and env.")
				.build();

		final Option biotransformerSequenceOption = Option.builder("q")
				.required(false)
				.hasArg(true)
				.argName("BioTransformer Sequence Option")
				.longOpt("bsequence")
				.desc("Define an ordered sequence of biotransformer/nr_of_steps to apply. Choose only from the following "
						+ "BioTranformer Types: allHuman, cyp450, ecbased, env, hgut, and phaseII. For instance, the "
						+ "following string representation describes a sequence of 2 steps of CYP450 metabolism, followed by 1 "
						+ "step of Human Gut metabolism, 1 step of Phase II, and 1 step of Environmental Microbial Degradation:\n"
						+ "'cyp450:2; hgut:1; phaseII:1; env:1'")
				.build();

		
		final Option nrOfStepsOption = Option.builder("s")
				.required(false)
				.hasArg(true)
				.argName("Number of steps")
				.longOpt("nsteps")
				.desc("The number of steps for the prediction. This option can be set by the user for the EC-based, CYP450, Phase II, and Environmental microbial biotransformers. The default value is 1.")
				.build();

		final Option smiInputOption = Option.builder("ismi")
				.required(false)
				.hasArg(true)
				.argName("SMILES Input")
				.longOpt("ismiles")
				.desc("The input, which can be a SMILES string")
				.build();

		final Option molInputOption = Option.builder("imol")
				.required(false)
				.hasArg(true)
				.argName("MOL Input")
				.longOpt("molinput")
				.desc("The input, which can be a Mol file")
				.build();
		
		final Option sdfInputOption = Option.builder("isdf")
				.required(false)
				.hasArg(true)
				.argName("Sdf Input")
				.longOpt("sdfinput")
				.desc("The input, which can be an SDF file.")
				.build();
		
		final Option csvOutputOption = Option.builder("ocsv")
				.required(false)
				.hasArg(true)
				.argName("Csv Output")
				.longOpt("csvoutput")
				.desc("Select this option to return CSV output(s). You must enter an output filename")
				.build();
		
		
		final Option sdfOutputOption = Option.builder("osdf")
				.required(false)
				.hasArg(true)
				.argName("Sdf Output")
				.longOpt("sdfoutput")
				.desc("Select this option to return SDF output(s). You must enter an output filename")
				.build();
		
		

		final Option annotateOption = Option.builder("a")
				.required(false)
				.hasArg(false)
				.argName("Annotate")
				.longOpt("annotate")
				.desc("Search PuChem for each product, and store with CID and synonyms, when available.")
				.build();
		
		final Option indentificationMassOption = Option.builder("m")
				.required(false)
				.hasArg(true)
				.argName("Masses")
				.longOpt("masses")
				.desc("Semicolon-separated list of masses of compounds to identify")
				.build();
		
		final Option indentificationFormulaMetadataOption = Option.builder("f")
				.required(false)
				.hasArg(true)
				.argName("Formulas")
				.longOpt("formulas")
				.desc("Semicolon-separated list of formulas of compounds to identify")
				.build();
				
		final Option massToleranceOption = Option.builder("t")
				.required(false)
				.hasArg(true)
				.argName("Mass Tolerance")
				.longOpt("mTolerance")
				.desc("Mass tolerance for metabolite identification (default is 0.01).")
				.build();

		final Option helpOption = Option.builder("h")
				.required(false)
				.hasArg(false)
				.argName("help")
				.longOpt("help")
				.desc("Prints the usage.")
				.build();
		
		final Option cypModeOpton = Option.builder("cm")
				.required(false)
				.hasArg(true)
				.argName("CYP450 Prediction Mode")
				.longOpt("cyp450mode")
				.desc("Specify the CYP450 predictoin Mode here: 1) CypReact + BioTransformer rules; 2) CyProduct only; 3) Combined: CypReact + BioTransformer rules + CyProducts.\nDefault mode is 1.")
				.build();
		
		final Options options = new Options();
		options.addOption(taskOption);
		options.addOption(biotransformerOption);
		options.addOption(biotransformerSequenceOption);
		options.addOption(nrOfStepsOption);
		options.addOption(smiInputOption);
		options.addOption(molInputOption);
		options.addOption(sdfInputOption);
		options.addOption(csvOutputOption);
		options.addOption(sdfOutputOption);
		options.addOption(indentificationMassOption);
		options.addOption(indentificationFormulaMetadataOption);	
		options.addOption(massToleranceOption);
		options.addOption(annotateOption);
		options.addOption(helpOption);
		options.addOption(cypModeOpton);

		return options;
	}
	
	public static CommandLine generateCommandLine(
			final Options options, final String[] commandLineArguments) throws ParseException{
		final CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine commandLine = null;
		;
		String header = "\nThis is the version " + Version.current + " of BioTransformer. BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota,"
				+ " as well as the soil/aquatic microbiota. BioTransformer also assists scientists in metabolite identification, based on the metabolism prediction. \n\n";
		
		String footer = "\n(* ) While the 'superbio' option runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.),"
				+ " the 'allHuman' option predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step."
				+ "(** ) For the environmental microbial biodegradation, all reactions (aerobic and anaerobic) are reported, and not only the aerobic biotransformations (as per default in the EAWAG BBD/PPS system)."
				+ "\n\n*********\n"
				+ "Examples:\n"
				+ "*********\n\n"
				+"1) To predict the biotransformation of a molecule from an SDF input using the human super transformer (option superbio) and annotate the metabolites with names and database IDs (from PubChem), run\n"
				+ "\n	java -jar biotransformer-" + Version.current +".jar -k pred -b superbio -isdf #{input file name} -osdf #{output file} -a."
				+ "\n\n2) To predict the 2-step biotransformation of Thymol (a monoterpene) using the human super transformer (option allHuman) using the SMILES input, cypMode = 3 (combined mode), and saving to a CSV file, run"
				+ "\n	java -jar biotransformer-" + Version.current +".jar  -k pred -b allHuman -ismi \"CC(C)C1=CC=C(C)C=C1O\" -ocsv #{replace with output file name} -s 2 -cm 3"
				+ "\n\n3) Identify all human metabolites (max depth = 2) of Epicatechin (\"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\") with masses 292.0946 Da and 304.0946 Da, with a mass tolerance of 0.01 Da."
				+ " Provide an annotation (Common name, synonyms, and PubChem CID), when available."
				+ "\n	java -jar biotransformer-" + Version.current +".jar  -k cid -b allHuman -ismi \"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\" -osdf #{replace with output file name} -s 2 -m \"292.0946;304.0946\" -t 0.01 -a"
				+ "\n	- DO NOT forget the quotes around the SMILES string or the list of masses"
				+ "\n"
				+ "4) Simulate an order sequence of metabolism of Atrazine (\"CCNC1=NC(=NC(=N1)Cl)NC(C)C\"), starting with two steps of Cyp450 oxidation, followed by one step of conjugation."
				+ "\n java -jar biotransformer-1.1.6.jar -ismi \"CCNC1=NC(=NC(=N1)Cl)NC(C)C\" -osdf ~/atrazine-sequence.sdf -k pred -q \"cyp450:2; phaseII:1\"\n"
				+ "\nTo report issues, provide feedback, or ask questions, please send an e-mail the following address: djoumbou@ualberta.ca\n\n"
				+ "BioTransformer is offered to the public as a freely acessible software package under the GNU License LGPL v3. Users are free"
				+ " to copy and redistribute the material in any medium or format. Moreover, they could modify, and build upon the material unfer "
				+ "the condition that they must give appropriate credit, provide links to the license, and indicate if changes were made. Furthermore, "
				+ "the above copyright notice and this permission notice must be included. Use and re-distribution of the these resources, in whole or in part, "
				+ "for commercial purposes requires explicit permission of the authors. We ask that all users of the BioTransformer software tool, the BioTransformer web server, "
				+ "or BioTransformerDB to cite the BioTransformer reference in any resulting publications, and to acknowledge the authors."
				+ "\n\n";

		HelpFormatter formatter = new HelpFormatter();

		
		try{
			commandLine = cmdLineParser.parse(options, commandLineArguments);
		}
		catch (MissingOptionException missingOptionException){
			
			if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
//				System.out.println("Version: " + Version.current);
				formatter.printHelp("\njava -jar biotransformer-" + Version.current +".jar", header, options, footer, true);
			}
			else {
				System.out.println(missingOptionException.getLocalizedMessage());
			}			
		}
		catch (ParseException parseException){
			System.out.println("Could not parse the command line arguments "
					+ Arrays.toString(commandLineArguments) + "\nfor the following reaons:" 
					+ parseException);		
		}
		return commandLine;
	}
	
	
	
	public static void main(String[] args) throws Exception{
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesParser	smiParser		= new SmilesParser(builder);

		Options options = generateOptions();
		CommandLine commandLine = generateCommandLine(options, args);
		IAtomContainer singleInput = null;
		String inputFileName = null;
		int nrOfSteps = 1;
		
		boolean annotate = false;
		String masses = null;
		String formulas = null;
		String metadata_input = null;
		Double defaultMassToleranceThreshold = 0.01;
		Double massToleranceThreshold = null;
		String task = null;
		FinderOption opt = null;
		String iFormat = null;
		String oFormat = null;
		String outputF = null;
		BiotransformerSequence biotransformerSeqeuence = null;
		double scoreThreshold = 0.5;
		int cyp450Mode = 1;
		
		int number_of_molecules = 0;
		int successful_predictions = 0;

		// Setting up the OptionsTobType hash map
		setOptionsTobType();
		
		if(Arrays.asList(args).contains("-a") || Arrays.asList(args)
				.contains("--annotate")){		
			annotate = true;
			
			System.out.println("\n\n=============================================================");
			System.out.println("Compounds will be annotated using the PubChem API. Make sure");
			System.out.println("to have a secure internet connection.");
			System.out.println("=============================================================\n");
		}
		
		if(Arrays.asList(args).contains("-ismi") || Arrays.asList(args)
				.contains("--ismiles")){
			iFormat = "smi";
		}
		else if(Arrays.asList(args).contains("-imol") || Arrays.asList(args)
				.contains("--sdfinput")){
			iFormat = "mol";
		}
		else if(Arrays.asList(args).contains("-isdf") || Arrays.asList(args)
				.contains("--sdfinput")){
			iFormat = "sdf";
		}
				
		if(Arrays.asList(args).contains("-ocsv") || Arrays.asList(args)
				.contains("--csvoutput")){
			oFormat = "csv";
			outputF = commandLine.getOptionValue("ocsv").trim();
		}
		else if(Arrays.asList(args).contains("-osdf") || Arrays.asList(args)
				.contains("sdfoutput")){
			oFormat = "sdf";
			outputF = commandLine.getOptionValue("osdf").trim();
		}		
		
		if(commandLine !=null){
			String mode = commandLine.getOptionValue("cm");
			if(mode != null){
//				System.out.println("mode.trim() = " + mode.trim() + "\n");
				cyp450Mode = Integer.parseInt(mode.trim());
				System.out.println("Cyp450Mode = " + cyp450Mode);
			}
								
			if(commandLine.getOptionValue("k") != null){
				task = commandLine.getOptionValue("k").trim();
				if( !(task.contentEquals("pred") || task.contentEquals("cid")) ){
					throw new IllegalArgumentException("Invalid task(\"" +  task + "\") entered. Enter either 'pred' (prediction) or 'cid' (compound identification)");
				}				
			}
			else {
				throw new IllegalArgumentException("\n\tThe task type is missing. You must select either 'pred' (prediction) or 'cid' (compound identification)");
			}

			if(commandLine.getOptionValue("m") != null){
				masses = commandLine.getOptionValue("m").trim();
				if(masses.length() == 0){
					throw new MissingArgumentException("\n\tPlease enter a list of monoisotopic masses.");
				}			
//				System.out.println("MASS: " + masses);
			}
			
			if(commandLine.getOptionValue("f") != null){
				formulas = commandLine.getOptionValue("f").trim();
				if(formulas.length() == 0){
					throw new MissingArgumentException("\n\tPlease enter a list of chemical formulas.");
				}			
//				System.out.println("MASS: " + masses);
			}
			
//			System.out.println("ANNOTATE: " + annotate);
			if(task.contentEquals("cid")){
				if(commandLine.getOptionValue("t") != null){
					if(commandLine.getOptionValue("t").trim().length() == 0){
						throw new MissingArgumentException("\n\tThe option '-t' was used but the mass tolerance threshold is missing. Add a value or omit '-t' to use the default value (" + defaultMassToleranceThreshold +")");
					}
					else{
						massToleranceThreshold = Double.valueOf(commandLine.getOptionValue("t").trim());	
//						System.out.println("massTolerance: " + massTolerance);
					}
				}
					
				if(masses == null && formulas == null){
					throw new IllegalArgumentException("\n\tIdentification metadata are missing. Please add a list of masses (-m) or a list of formulas (-r)");
				}
				else if(formulas != null) {
					if(masses != null) {
						throw new IllegalArgumentException("\tList of masses and formulas provided. Please exclusively provide either a list of masses (-m) or a list of formulas (-r)");											
					}
					else if(massToleranceThreshold != null) {
						throw new IllegalArgumentException("\n\tA mass tolerance threshold is accepted only for mass-based, and not formula-based identification tasks. Please remove this argument.");											
						
					}
					opt = FinderOption.FORMULA;
					metadata_input = formulas;
//					System.out.println(opt + "\t" + metadata_input);
					/**
					 * Here, there user has not provided a mass threshold tolerance, but we set it to the default, just to avoid the NullPointerException
					 * The massToleranceThreshold argument will not be used for the FORMULA option.
					 */
					massToleranceThreshold = defaultMassToleranceThreshold;
				}
				else if(masses != null){
					opt = FinderOption.MASS;
					metadata_input = masses;
					if(massToleranceThreshold == null) {
						massToleranceThreshold = defaultMassToleranceThreshold;
					}
				}
			}

			final String biotransformerType = commandLine.getOptionValue("b");
			
			if(commandLine.getOptionValue("s") != null){
				nrOfSteps = Integer.valueOf(commandLine.getOptionValue("s"));
//				System.out.println("nrOfSteps: " + nrOfSteps);
			}	
			
			String bseq = commandLine.getOptionValue("q");
			if(bseq != null) {
				if(bseq.contentEquals("")) {
					throw new IllegalArgumentException("\n\tIllegalArgumentException. Please provide a valid biotranformer sequence.");				
				}
				else if(bseq.length()>1) {
					biotransformerSeqeuence = new BiotransformerSequence(bseq);
				}		
			}
			
			if(biotransformerType != null) {
				if(! optionsToBtTypes.containsKey(biotransformerType.toLowerCase())) {
					throw new IllegalArgumentException("IllegalArgumentException: The biotransformer type is invalid. Select between allHuman, cyp450, ecbased, env, hgut, phaseII, and superbio");			
				}
				if(commandLine.getOptionValue("q") != null) {
					throw new IllegalArgumentException("IllegalArgumentException: The parameters '-b' and '-seq' are mutually "
							+ "exclusive. While '-b' describes a specific biotransformer type, '-q' describes an order sequence of biotransformers to be applied.");
				}
			}
			
			else if(biotransformerType == null && commandLine.getOptionValue("q") == null){
				throw new MissingArgumentException("MissingArgumentException: Please select either '-b' for a bitransformer type, or '-seq' "
						+ "for a biotransformer sequence.");
			}
			if(commandLine.getOptionValue("s") != null && commandLine.getOptionValue("q") != null) {
					throw new IllegalArgumentException("IllegalArgumentException: The parameters '-s' and '-seq' are mutually "
							+ "exclusive. '-q' describes an order sequence of biotransformers to be applied, each for a specific number of steps");
			}					

			
			if(iFormat.contentEquals("smi")){
				String smi = commandLine.getOptionValue("ismi");
				singleInput = smiParser.parseSmiles(smi);			
			}
			else if(iFormat.contentEquals("mol")){
				inputFileName = commandLine.getOptionValue("imol");
				if(inputFileName == null){
					throw new MissingOptionException("\n\tPlease specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
//				if(outputF == null){
//					throw new MissingOptionException("A destination folder must be specified when your query molecules are provided in a file (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
//				}
			}
			
			else if(iFormat.contentEquals("sdf")){
				inputFileName = commandLine.getOptionValue("isdf");
				if(inputFileName == null){
					throw new MissingOptionException("\n\tPlease specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}
			else {
				throw new IllegalArgumentException("\n\tInvalid input format option(" + iFormat + ") entered. It must be one of 'ismi','imol', or 'isdf'. Type java -jar biotransformer-" + Version.current +".jar --help for help.");
			}

			if(oFormat.contentEquals("csv")){
				outputF = commandLine.getOptionValue("ocsv");
				if(outputF == null){
					throw new MissingOptionException("\n\tPlease specify an output file name (CSV). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}			
			else if(oFormat.contentEquals("sdf")){
				outputF = commandLine.getOptionValue("osdf");
				if(outputF == null){
					throw new MissingOptionException("\n\tPlease specify an output file name (SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}
			else {
				throw new IllegalArgumentException("\n\tInvalid output format option(" + oFormat + ") entered. It must be one of 'ocsv' or 'osdf'. Type java -jar biotransformer-" + Version.current +".jar --help for help.");
			}

			if(task.contentEquals("cid")){
				
				if(metadata_input != null){
					
					if(biotransformerType != null) {
						if(optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.ALLHUMAN){
							
							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
							
							
							if (singleInput !=null){
								MetaboliteFinder mtf = new MetaboliteFinder();
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								if(oFormat.contentEquals("csv")){
									mtf.findAllHumanMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt, cyp450Mode);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findAllHumanMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt, cyp450Mode);
								}
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findAllHumanMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, opt, cyp450Mode));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}
								
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
							}	
							
						}
						else if(optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.ENV){
							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
							
							
							if (singleInput !=null){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								MetaboliteFinder mtf = new MetaboliteFinder();
								
								if(oFormat.contentEquals("csv")){
									mtf.findAllEnvMicroMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findAllEnvMicroMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findAllEnvMicroMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, opt));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}							
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
							}					
						}
	
						else if(optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.SUPERBIO){

							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
								
							if (singleInput !=null){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								MetaboliteFinder mtf = new MetaboliteFinder();
								
								if(oFormat.contentEquals("csv")){
									mtf.findSuperbioMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, annotate, outputF, opt, cyp450Mode);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findSuperbioMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, annotate, outputF, opt, cyp450Mode);
								}
								
								
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findSuperbioMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, annotate, opt, cyp450Mode));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}							
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
	
							}					
						}
						else{
							throw new IllegalArgumentException("\n\tFor metabolite identification, the biotransformer type must be either allHuman, superbio, or env.");
						}
					}
					else if(biotransformerSeqeuence != null) {
						String[] mArr = metadata_input.trim().split(";");
						ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
						
						for(int k = 0; k < mArr.length; k++){
							try{
								dmassesOrFormulas.add(mArr[k].trim());
							}
							catch(Exception e){
								System.err.println(e.getMessage());
							}
						}		
						IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
						MetaboliteFinder mtf = new MetaboliteFinder();
						if (singleInput !=null){
							number_of_molecules++;
							
							metabolites.add(mtf.findMetabolitesFromSequence(singleInput, 
									biotransformerSeqeuence, dmassesOrFormulas, massToleranceThreshold, annotate, opt, cyp450Mode));

							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);					

								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										 metabolites.add(mtf.findMetabolitesFromSequence(atc, 
													biotransformerSeqeuence, dmassesOrFormulas, massToleranceThreshold, annotate, opt, cyp450Mode));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}	
								}
							}				
						}
						
						if(oFormat.contentEquals("csv")){
							FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
						}
						else if(oFormat.contentEquals("sdf")){							
							SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
							sdfWriter.write(metabolites);
							sdfWriter.close();
						}
					}
				}
				else{
					throw new IllegalArgumentException("\n\tFor metabolite identification, you must enter a list of masses, and/or molecular formulas");
				}
				
			}
			else if(task.contentEquals("pred")){
				if(biotransformerType != null) {
					
					if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.CYP450){
						Cyp450BTransformer cyp450bt = new Cyp450BTransformer(BioSystemName.HUMAN);
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//						int cyp450Mode = 1;
//						//Set the default mode as 1
//						if(mode!=null){
//							cyp450Mode = Integer.parseInt(mode);
//						}
						
						if (singleInput !=null){
							number_of_molecules++;
							//biotransformations = cyp450bt.predictCyp450BiotransformationChain(singleInput, true, true, nrOfSteps, scoreThreshold);
//							biotransformations = cyp450bt.predictCyp450BiotransformationsByMode(singleInput, cyp450Mode, true, true, scoreThreshold);
							biotransformations = cyp450bt.predictCyp450BiotransformationChainByMode(singleInput, true, true, nrOfSteps, scoreThreshold, cyp450Mode);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);					

								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										//biotransformations.addAll(cyp450bt.predictCyp450BiotransformationChain(atc, true, true, nrOfSteps, scoreThreshold));
										//biotransformations.addAll(cyp450bt.predictCyp450BiotransformationsByMode(singleInput, cyp450Mode, true, true, scoreThreshold));
										biotransformations.addAll(cyp450bt.predictCyp450BiotransformationChainByMode(atc, true, true, nrOfSteps, scoreThreshold, cyp450Mode));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
								}						
								
							}				
						}
										
						if(oFormat.contentEquals("csv")){
							cyp450bt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							cyp450bt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.ECBASED) {
						ECBasedBTransformer ecbt =  new ECBasedBTransformer(BioSystemName.HUMAN);
			
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = ecbt.simulateECBasedMetabolismChain(singleInput, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);

								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										biotransformations.addAll(ecbt.simulateECBasedMetabolismChain(atc, true, true, nrOfSteps, scoreThreshold));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
										continue;
									}										
								}	
							
							}				
						}
										
						if(oFormat.contentEquals("csv")){
							ecbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							ecbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.HGUT){
						HGutBTransformer hgut = new HGutBTransformer();
									
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hgut.simulateGutMicrobialMetabolism(singleInput, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {			
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hgut.inchiGenFactory);					
							if (containers.getAtomContainerCount()>0){
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										biotransformations.addAll(hgut.simulateGutMicrobialMetabolism(atc, true, true, nrOfSteps, scoreThreshold));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
										continue;
									}
									
								}	
							}				
						}	
						if(oFormat.contentEquals("csv")){
							hgut.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hgut.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}						
					}
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.PHASEII){
						Phase2BTransformer phase2b = new Phase2BTransformer(BioSystemName.HUMAN);
						
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(singleInput,
									true, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, phase2b.inchiGenFactory);					
							for(IAtomContainer atc : containers.atomContainers()){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								try {
									biotransformations.addAll(phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(atc, true, true, true, nrOfSteps, scoreThreshold));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}
						
						}
										
						if(oFormat.contentEquals("csv")){
							phase2b.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							phase2b.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.SUPERBIO){
						if(nrOfSteps!=12){
							System.out.println("\n\n=======>The configutration is set for this simulation. No need to set a number of steps for the super human transformer.\n\n");
						}
//						int cyp450Mode = 1;
//						if(mode!=null){
//							cyp450Mode = Integer.parseInt(mode);
//						}
						HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
						
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hsbt.simulateHumanSuperbioMetabolism(singleInput, scoreThreshold, cyp450Mode);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hsbt.getInChIGenFactory());					
							for(IAtomContainer atc : containers.atomContainers()){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								try {
									biotransformations.addAll(hsbt.simulateHumanSuperbioMetabolism(atc,scoreThreshold, cyp450Mode));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}									
							}							
						}			
		
						if(oFormat.contentEquals("csv")){
							hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
		
					}
					
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.ALLHUMAN){
						HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
//						int cyp450Mode = 1;
//						if(mode!=null){
//							cyp450Mode = Integer.parseInt(mode);
//						}
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hsbt.predictAllHumanBiotransformationChain(singleInput, nrOfSteps, scoreThreshold, cyp450Mode);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hsbt.getInChIGenFactory());					
							System.out.println("Nr. of molecules: " + containers.getAtomContainerCount() );
							for(IAtomContainer atc : containers.atomContainers()){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								try {
									biotransformations.addAll(hsbt.predictAllHumanBiotransformationChain(atc, nrOfSteps, scoreThreshold, cyp450Mode));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}									
							}							
						}			
		
						if(oFormat.contentEquals("csv")){
							hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					
					else if (optionsToBtTypes.get(biotransformerType.toLowerCase()) == bType.ENV){
						EnvMicroBTransformer ebt = new EnvMicroBTransformer();
						
						if (singleInput !=null){
							number_of_molecules++;
							if(oFormat.contentEquals("csv")){
		
								ebt.simulateEnvMicrobialDegradationAndSaveToCSV(singleInput, true, true, nrOfSteps, scoreThreshold, outputF, annotate);
								successful_predictions++;
							}
							else if(oFormat.contentEquals("sdf")){
								ebt.simulateEnvMicrobialDegradationAndSaveToSDF(singleInput, true, true, nrOfSteps, scoreThreshold, outputF, annotate);
								successful_predictions++;
							}
						}
						else {
							ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							for(IAtomContainer atc : containers.atomContainers()){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								try {
									biotransformations.addAll(ebt.applyEnvMicrobialTransformationsChain(atc, true, true, nrOfSteps, scoreThreshold));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}						
							
							if(oFormat.contentEquals("csv")){
								ebt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
							}
							else if(oFormat.contentEquals("sdf")){
								ebt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
							}					
						}						
					}
				}
				
				else if(biotransformerSeqeuence != null) {
					ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
					if (singleInput !=null){
						number_of_molecules++;
						biotransformations = biotransformerSeqeuence.runSequence(singleInput, scoreThreshold, cyp450Mode);
						successful_predictions++;
					}
					else {
						IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
						if (containers.getAtomContainerCount()>0){
							containers = FileUtilities.parseSdf(inputFileName);					

							for(IAtomContainer atc : containers.atomContainers()){
								number_of_molecules++;
								System.out.println("\nMolecule no. " + number_of_molecules + ": " + atc.getProperty(CDKConstants.TITLE));
								try {
									biotransformations.addAll(biotransformerSeqeuence.runSequence(atc, scoreThreshold, cyp450Mode));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
								}
								
							}						
							
						}				
					}
					UniversalBioTransformer hsbt = new UniversalBioTransformer();
					System.out.println("Format: " + oFormat);
					if(oFormat.contentEquals("csv")){
						System.out.println("saving to "+ outputF);
						hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
					}
					else if(oFormat.contentEquals("sdf")){
						System.out.println("saving to "+ outputF);
						hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
					}
					
				}
				System.out.println("Successfully completed metabolism prediction for " + successful_predictions + " out of " + number_of_molecules + " molecule(s).");
			}
			else {
				throw new IllegalArgumentException("You entered an invalid biotransformer option.\n"
						+ "Choose one of the following: ecbased, cyp450, hgut, phaseII, superbio, allHuman, or env.\nType java -jar biotransformer-" + Version.current +".jar --help for help.");
			}			
		}		
	}

}




