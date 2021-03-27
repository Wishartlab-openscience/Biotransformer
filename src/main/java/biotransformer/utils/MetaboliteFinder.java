/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
//import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.apache.commons.lang3.StringUtils;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
//import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
//import org.openscience.cdk.formula.MolecularFormulaChecker;
//import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import biotransformer.btransformers.Biotransformer.bType;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.transformation.Biotransformation;
import exception.BioTransformerException;

//

public class MetaboliteFinder{
	
	public HumanSuperBioTransformer hsbt;
	public EnvMicroBTransformer ebt;
	public UniversalBioTransformer ubt;
	public InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
	protected IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
	
	
	public MetaboliteFinder() throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException, CDKException{
		// TODO Auto-generated constructor stub
		hsbt 		= new HumanSuperBioTransformer();
		ebt 		= new EnvMicroBTransformer();
		ubt 		= new UniversalBioTransformer();
	
	}
	
	public enum FinderOption {
		MASS, FORMULA, MASSFORMULA
	}
		
	public void findAllHumanMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, String outputFileName, FinderOption opt, int cyp450Mode) throws Exception{
		
		IAtomContainerSet metabolites = findAllHumanMetabolites(startingCompound, mass_formulas, massTolerance, nrOfSteps, annotate, opt, cyp450Mode);
		
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(metabolites);
		sdfWriter.close();
	}
	
	public void findAllHumanMetabolitesToCSV(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, String outputFileName, FinderOption opt, int cyp450Mode) throws Exception{		
		IAtomContainerSet metabolites = findAllHumanMetabolites(startingCompound, mass_formulas, massTolerance, nrOfSteps, annotate, opt, cyp450Mode);
		FileUtilities.saveAtomContainerSetToCSV(metabolites, outputFileName);
	}
	
	public IAtomContainerSet findAllHumanMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, FinderOption opt, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> allBiotransformations = new ArrayList<Biotransformation>();
		
//		ArrayList<String> mass_formulas_cleaned = new ArrayList<String>();		
//		for(int s=0; s<mass_formulas.size(); s++){
//			
//		}
//		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas_cleaned.clone();
		
		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas.clone();
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		LinkedHashMap<String, ArrayList<Biotransformation>> substratesToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, ArrayList<Biotransformation>>productsToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, IAtomContainerSet> massFormulaToMolecules = new LinkedHashMap<String, IAtomContainerSet>();		
		LinkedHashMap<String, IAtomContainer> inchikeyToContainers = new LinkedHashMap<String, IAtomContainer>();
		LinkedHashMap<String, IAtomContainerSet> compoundToParents = new LinkedHashMap<String, IAtomContainerSet>();
		
		IAtomContainer startingCompoundStandardized = this.hsbt.standardizeMoleculeWithCopy(startingCompound);
			
		ChemStructureExplorer.addInChIandKey(startingCompoundStandardized);
		inchikeyToContainers.put((String) startingCompoundStandardized.getProperty("InChIKey"), startingCompoundStandardized);

		
		IAtomContainerSet pathways = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

		
		
		for(int k=0; k < mass_formulas.size(); k++){
			 massFormulaToMolecules.put(mass_formulas.get(k), DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class));
		}
		
		IAtomContainerSet currentCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentCompounds.addAtomContainer(startingCompound);
		int counter = 1;
		
		while(nrOfSteps >0 && remainingMassFormulas.size()>0){
			
			ArrayList<Biotransformation> currentBioT = this.hsbt.simulateOneStepAllHuman(currentCompounds, 0.5, cyp450Mode);
			IAtomContainerSet extractedCompounds = this.hsbt.extractProductsFromBiotransformationsWithTransformationData(currentBioT, annotate);
			allBiotransformations.addAll(currentBioT);
			currentCompounds.removeAllAtomContainers();
			currentCompounds.add(ChemStructureExplorer.uniquefy(extractedCompounds));
			
			for(IAtomContainer a : currentCompounds.atomContainers()){
				inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
			}
			
			for(IAtomContainer a : currentCompounds.atomContainers()){
//				inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
				String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//				System.out.println("precursorInChIKey: " + precursorInChIKey);
//				System.out.print((String) a.getProperty("InChIKey"));
//				System.out.println(" >> " + (String)a.getProperty("Precursor InChIKey"));
//				System.out.println(inchikeyToContainers.containsKey(precursorInChIKey));
				
//				if(inchikeyToContainers.containsKey(precursorInChIKey)){
					if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//						compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
						compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(a);
					}
					else {
						IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//						s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//						compoundToParents.put((String) a.getProperty("InChIKey"), s);	
						s.addAtomContainer(a);
						compoundToParents.put((String) a.getProperty("InChIKey"), s);
					}
//				}
//				System.out.println(compoundToParents.get((String) a.getProperty("InChIKey")) == null);
			}			
			
			
			for(int k=0; k < mass_formulas.size(); k++){				
				if(opt == FinderOption.FORMULA){
					IMolecularFormula formula1 = MolecularFormulaManipulator.getMolecularFormula(mass_formulas.get(k), this.builder);
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if(a.getProperty("Molecular formula") == null){
							a.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(a));
						}
						System.out.println("A: " + MolecularFormulaManipulator.getString(formula1));
						System.out.println("B: " + a.getProperty("Molecular formula").toString());
						if(MolecularFormulaManipulator.getString(formula1).contentEquals(a.getProperty("Molecular formula").toString().trim())){
							System.out.println(a.getProperty("Molecular formula").toString().trim());
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				else if(opt == FinderOption.MASS){
					Double mass = Double.valueOf(mass_formulas.get(k));
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if( Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass) <= massTolerance){
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				else if(opt == FinderOption.MASSFORMULA){
					String[] mf = mass_formulas.get(k).split(":");
					IMolecularFormula formula2 = MolecularFormulaManipulator.getMolecularFormula(mf[1], this.builder);
					Double mass2 = Double.valueOf(mf[0]);
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if( (MolecularFormulaManipulator.getString(formula2).contentEquals(mass_formulas.get(k))) && 
								(Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass2) <= massTolerance)){
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				if(massFormulaToMolecules.get(mass_formulas.get(k)).getAtomContainerCount() > 0){
					remainingMassFormulas.remove(mass_formulas.get(k));
				}
				
			}			
			nrOfSteps--;
			counter++;
//			System.err.println("Remaining masses: " + remainingMasses.size() + "\n");
//			System.err.println("\n\n\n\n\nNumber of possibilities: " + filteredCompounds.getAtomContainerCount() );
		}
		
		for(int i = 0; i < allBiotransformations.size(); i++){
			Biotransformation bt = allBiotransformations.get(i);
			for(IAtomContainer ac : bt.getSubstrates().atomContainers()){
				String ikey = ac.getProperty("InChIKey");
				if(ikey != null){
					if(substratesToBio.containsKey(ikey)){
						substratesToBio.get(ikey).add(bt);
					}else{
						substratesToBio.put(ikey,new ArrayList<Biotransformation>());
						substratesToBio.get(ikey).add(bt);
					}
				}
			}

			for(IAtomContainer at : bt.getProducts().atomContainers()){
				String ikey = at.getProperty("InChIKey");
				if(ikey != null){
					if(productsToBio.containsKey(ikey)){
						productsToBio.get(ikey).add(bt);
					}else{
						productsToBio.put(ikey,new ArrayList<Biotransformation>());
						productsToBio.get(ikey).add(bt);
					}
				}
			}
		}
		
		IAtomContainerSet uniqueFilteredCompounds = ChemStructureExplorer.uniquefy(filteredCompounds);
		for(IAtomContainer a : uniqueFilteredCompounds.atomContainers()){
//			System.out.println(this.ebt.smiGen.create(a) + " - " + a.getProperty("Major Isotope Mass"));
			pathways.addAtomContainer(findPathway(startingCompoundStandardized, a, compoundToParents, inchikeyToContainers, annotate));
		}		
				
		return pathways;
	}
	
	public void findSuperbioMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, boolean annotate, String outputFileName, FinderOption opt, int cyp450Mode) throws Exception{
		
		IAtomContainerSet metabolites = findSuperbioMetabolites(startingCompound, mass_formulas, massTolerance, annotate, opt, cyp450Mode);		
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(metabolites);
		sdfWriter.close();
	}
	
	
	public void findSuperbioMetabolitesToCSV(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, boolean annotate, String outputFileName, FinderOption opt, int cyp450Mode) throws Exception{		
		IAtomContainerSet metabolites = findSuperbioMetabolites(startingCompound, mass_formulas, massTolerance, annotate, opt, cyp450Mode);
		FileUtilities.saveAtomContainerSetToCSV(metabolites, outputFileName);
	}
	
	public IAtomContainerSet findSuperbioMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, boolean annotate, FinderOption opt, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> allBiotransformations = new ArrayList<Biotransformation>();
		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas.clone();
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet pathways = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		LinkedHashMap<String, ArrayList<Biotransformation>> substratesToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, ArrayList<Biotransformation>>productsToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, IAtomContainerSet> massFormulaToMolecules = new LinkedHashMap<String, IAtomContainerSet>();		
		LinkedHashMap<String, IAtomContainer> inchikeyToContainers = new LinkedHashMap<String, IAtomContainer>();
		LinkedHashMap<String, IAtomContainerSet> compoundToParents = new LinkedHashMap<String, IAtomContainerSet>();
		
		for(int k=0; k < mass_formulas.size(); k++){
			massFormulaToMolecules.put(mass_formulas.get(k), DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class));
		}
		
		IAtomContainer startingCompoundStandardized = this.hsbt.standardizeMoleculeWithCopy(startingCompound);
		ChemStructureExplorer.addInChIandKey(startingCompoundStandardized);
		inchikeyToContainers.put((String) startingCompoundStandardized.getProperty("InChIKey"), startingCompoundStandardized);
		
		
//		System.out.println("STARTING COMPOUND INCHIKEY: " + startingCompoundStandardized.getProperty("InChIKey"));
		
		allBiotransformations = this.hsbt.simulateHumanSuperbioMetabolism(startingCompoundStandardized, cyp450Mode);
		IAtomContainerSet extractedCompounds = this.hsbt.extractProductsFromBiotransformationsWithTransformationData(allBiotransformations, annotate);
		
		for(IAtomContainer a : extractedCompounds.atomContainers()){
			inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
//			String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//			System.out.println("precursorInChIKey: " + precursorInChIKey);
//			System.out.println(a.getProperties());
//			if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//				compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//			}
//			else {
//				IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//				s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//				compoundToParents.put((String) a.getProperty("InChIKey"), s);
//				
//			}
		}
		
		for(IAtomContainer a : extractedCompounds.atomContainers()){
//			inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
			String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//			System.out.println("precursorInChIKey: " + precursorInChIKey);
//			System.out.print((String) a.getProperty("InChIKey"));
//			System.out.println(" >> " + (String)a.getProperty("Precursor InChIKey"));
//			System.out.println(inchikeyToContainers.containsKey(precursorInChIKey));
			
//			if(inchikeyToContainers.containsKey(precursorInChIKey)){
				if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//					compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
					compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(a);
				}
				else {
					IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//					s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//					compoundToParents.put((String) a.getProperty("InChIKey"), s);	
					s.addAtomContainer(a);
					compoundToParents.put((String) a.getProperty("InChIKey"), s);
				}
//			}
//			System.out.println(compoundToParents.get((String) a.getProperty("InChIKey")) == null);
		}

		for(int k=0; k < mass_formulas.size(); k++){
			
			if(opt == FinderOption.FORMULA){
				IMolecularFormula formula1 = MolecularFormulaManipulator.getMolecularFormula(mass_formulas.get(k), this.builder);
				for(IAtomContainer a : extractedCompounds.atomContainers()){ 
					if(MolecularFormulaManipulator.getString(formula1).contentEquals(a.getProperty("Molecular formula").toString().trim())){
//						System.out.println(a.getProperty("Molecular formula").toString().trim());
						massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
						filteredCompounds.addAtomContainer(a);
//						break;
					}
				}					
			}
			else if(opt == FinderOption.MASS){
				Double mass = Double.valueOf(mass_formulas.get(k));
				for(IAtomContainer a : extractedCompounds.atomContainers()){ 
					if( Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass) <= massTolerance){
						massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
						filteredCompounds.addAtomContainer(a);
//						break;
					}
				}					
			}
			else if(opt == FinderOption.MASSFORMULA){
				String[] mf = mass_formulas.get(k).split(":");
				IMolecularFormula formula2 = MolecularFormulaManipulator.getMolecularFormula(mf[1], this.builder);
				Double mass2 = Double.valueOf(mf[0]);
				for(IAtomContainer a : extractedCompounds.atomContainers()){ 
					if( (MolecularFormulaManipulator.getString(formula2).contentEquals(mass_formulas.get(k))) && 
							(Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass2) <= massTolerance)){
						massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
						filteredCompounds.addAtomContainer(a);
//						break;
					}
				}					
			}
			
			if(massFormulaToMolecules.get(mass_formulas.get(k)).getAtomContainerCount() > 0){
				remainingMassFormulas.remove(mass_formulas.get(k));
			}
			
		}
		
	
		for(int i = 0; i < allBiotransformations.size(); i++){
			Biotransformation bt = allBiotransformations.get(i);
			for(IAtomContainer ac : bt.getSubstrates().atomContainers()){
				String ikey = ac.getProperty("InChIKey");
				if(ikey != null){
					if(substratesToBio.containsKey(ikey)){
						substratesToBio.get(ikey).add(bt);
					}else{
						substratesToBio.put(ikey,new ArrayList<Biotransformation>());
						substratesToBio.get(ikey).add(bt);
					}
				}
			}

			for(IAtomContainer at : bt.getProducts().atomContainers()){
				String ikey = at.getProperty("InChIKey");
				if(ikey != null){
					if(productsToBio.containsKey(ikey)){
						productsToBio.get(ikey).add(bt);
					}else{
						productsToBio.put(ikey,new ArrayList<Biotransformation>());
						productsToBio.get(ikey).add(bt);
					}
				}
			}
		}
		
		
//		System.out.println("Number of filtered compounds: " + filteredCompounds.getAtomContainerCount());
//		System.out.println("Number of explained masses: " + (mass_formulas.size() - remainingMassFormulas.size()) +  " out of " + mass_formulas.size() + "." );
//		for(IAtomContainer m : filteredCompounds.atomContainers()){
//			System.out.print(m.getProperty("InChIKey") + " ");
//			System.out.print(compoundToParents.get((String)m.getProperty("InChIKey")).getAtomContainerCount());
//			System.out.print(" parents\n");
//		}
		IAtomContainerSet uniqueFilteredCompounds = ChemStructureExplorer.uniquefy(filteredCompounds);
		System.out.println("Number of filtered compounds: " + uniqueFilteredCompounds.getAtomContainerCount());
		System.out.println("Number of explained masses: " + (mass_formulas.size() - remainingMassFormulas.size()) +  " out of " + mass_formulas.size() + "." );
		ArrayList<String> explained_masses_formulas = (ArrayList<String>) mass_formulas.clone();
		 explained_masses_formulas.removeAll(remainingMassFormulas);
		
		for(int x = 0; x <  explained_masses_formulas.size(); x++){
			System.out.println(explained_masses_formulas.get(x));
		}
		
		for(IAtomContainer a : uniqueFilteredCompounds.atomContainers()){
//			System.out.println(this.ebt.smiGen.create(a) + " - " + a.getProperty("Major Isotope Mass"));
			pathways.addAtomContainer(findPathway(startingCompoundStandardized, a, compoundToParents, inchikeyToContainers, annotate));
		}		
		
		return pathways;	
	}
	
	
	public void findAllEnvMicroMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, String outputFileName, FinderOption opt) throws Exception{
		
		IAtomContainerSet metabolites = findAllEnvMicroMetabolites(startingCompound, mass_formulas, massTolerance, nrOfSteps, annotate, opt);
		
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(metabolites);
		sdfWriter.close();
	}
	
	
	public void findAllEnvMicroMetabolitesToCSV(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, String outputFileName, FinderOption opt) throws Exception{		
		IAtomContainerSet metabolites = findAllEnvMicroMetabolites(startingCompound, mass_formulas, massTolerance, nrOfSteps, annotate, opt);
		FileUtilities.saveAtomContainerSetToCSV(metabolites, outputFileName);
	}
	
	public IAtomContainerSet findAllEnvMicroMetabolites(IAtomContainer startingCompound, ArrayList<String>mass_formulas, double massTolerance, int nrOfSteps, boolean annotate, FinderOption opt) throws Exception{
		ArrayList<Biotransformation> allBiotransformations = new ArrayList<Biotransformation>();
		
//		ArrayList<String> mass_formulas_cleaned = new ArrayList<String>();		
//		for(int s=0; s<mass_formulas.size(); s++){
//			
//		}
//		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas_cleaned.clone();
		
		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas.clone();
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		LinkedHashMap<String, ArrayList<Biotransformation>> substratesToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, ArrayList<Biotransformation>>productsToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, IAtomContainerSet> massFormulaToMolecules = new LinkedHashMap<String, IAtomContainerSet>();		
		LinkedHashMap<String, IAtomContainer> inchikeyToContainers = new LinkedHashMap<String, IAtomContainer>();
		LinkedHashMap<String, IAtomContainerSet> compoundToParents = new LinkedHashMap<String, IAtomContainerSet>();
		
		IAtomContainer startingCompoundStandardized = this.hsbt.standardizeMoleculeWithCopy(startingCompound);
			
		ChemStructureExplorer.addInChIandKey(startingCompoundStandardized);
		inchikeyToContainers.put((String) startingCompoundStandardized.getProperty("InChIKey"), startingCompoundStandardized);

		
		IAtomContainerSet pathways = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

		
		
		for(int k=0; k < mass_formulas.size(); k++){
			 massFormulaToMolecules.put(mass_formulas.get(k), DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class));
		}
		
		IAtomContainerSet currentCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentCompounds.addAtomContainer(startingCompound);
		int counter = 1;
		
		while(nrOfSteps >0 && remainingMassFormulas.size()>0){
			
//			System.out.println("nrOfSteps: " + nrOfSteps);
			ArrayList<Biotransformation> currentBioT = this.ebt.applyEnvMicrobialTransformations(currentCompounds, true, true, 0.5);
			
//			System.out.println("currentBioT: " + currentBioT.size());
			IAtomContainerSet extractedCompounds = this.ebt.extractProductsFromBiotransformationsWithTransformationData(currentBioT, annotate);
			allBiotransformations.addAll(currentBioT);
			currentCompounds.removeAllAtomContainers();
			currentCompounds.add(ChemStructureExplorer.uniquefy(extractedCompounds));
			
//			for(IAtomContainer a : currentCompounds.atomContainers()){
//				inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
//			}
	
			for(IAtomContainer a : currentCompounds.atomContainers()){
				inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
//				String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//				System.out.println("precursorInChIKey: " + precursorInChIKey);
//				System.out.println(a.getProperties());
//				if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//					compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//				}
//				else {
//					IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//					s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//					compoundToParents.put((String) a.getProperty("InChIKey"), s);
//					
//				}
			}
			
			for(IAtomContainer a : currentCompounds.atomContainers()){
//				inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
				String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//				System.out.println("precursorInChIKey: " + precursorInChIKey);
//				System.out.print((String) a.getProperty("InChIKey"));
//				System.out.println(" >> " + (String)a.getProperty("Precursor InChIKey"));
//				System.out.println(inchikeyToContainers.containsKey(precursorInChIKey));
				
//				if(inchikeyToContainers.containsKey(precursorInChIKey)){
					if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//						compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
						compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(a);
					}
					else {
						IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//						s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//						compoundToParents.put((String) a.getProperty("InChIKey"), s);	
						s.addAtomContainer(a);
						compoundToParents.put((String) a.getProperty("InChIKey"), s);
					}
//				}
//				System.out.println(compoundToParents.get((String) a.getProperty("InChIKey")) == null);
			}			
			
			
			for(int k=0; k < mass_formulas.size(); k++){				
				if(opt == FinderOption.FORMULA){
					IMolecularFormula formula1 = MolecularFormulaManipulator.getMolecularFormula(mass_formulas.get(k), this.builder);
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if(MolecularFormulaManipulator.getString(formula1).contentEquals(a.getProperty("Molecular formula").toString().trim())){
							System.out.println(a.getProperty("Molecular formula").toString().trim());
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				else if(opt == FinderOption.MASS){
					Double mass = Double.valueOf(mass_formulas.get(k));
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if( Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass) <= massTolerance){
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				else if(opt == FinderOption.MASSFORMULA){
					String[] mf = mass_formulas.get(k).split(":");
					IMolecularFormula formula2 = MolecularFormulaManipulator.getMolecularFormula(mf[1], this.builder);
					Double mass2 = Double.valueOf(mf[0]);
					for(IAtomContainer a : extractedCompounds.atomContainers()){ 
						if( (MolecularFormulaManipulator.getString(formula2).contentEquals(mass_formulas.get(k))) && 
								(Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass2) <= massTolerance)){
							massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
							filteredCompounds.addAtomContainer(a);
//							break;
						}
					}					
				}
				if(massFormulaToMolecules.get(mass_formulas.get(k)).getAtomContainerCount() > 0){
					remainingMassFormulas.remove(mass_formulas.get(k));
				}
				
			}			
			nrOfSteps--;
			counter++;
//			System.err.println("Remaining masses: " + remainingMasses.size() + "\n");
//			System.err.println("\n\n\n\n\nNumber of possibilities: " + filteredCompounds.getAtomContainerCount() );
		}
		
		for(int i = 0; i < allBiotransformations.size(); i++){
			Biotransformation bt = allBiotransformations.get(i);
			for(IAtomContainer ac : bt.getSubstrates().atomContainers()){
				String ikey = ac.getProperty("InChIKey");
				if(ikey != null){
					if(substratesToBio.containsKey(ikey)){
						substratesToBio.get(ikey).add(bt);
					}else{
						substratesToBio.put(ikey,new ArrayList<Biotransformation>());
						substratesToBio.get(ikey).add(bt);
					}
				}
			}

			for(IAtomContainer at : bt.getProducts().atomContainers()){
				String ikey = at.getProperty("InChIKey");
				if(ikey != null){
					if(productsToBio.containsKey(ikey)){
						productsToBio.get(ikey).add(bt);
					}else{
						productsToBio.put(ikey,new ArrayList<Biotransformation>());
						productsToBio.get(ikey).add(bt);
					}
				}
			}
		}
		
		IAtomContainerSet uniqueFilteredCompounds = ChemStructureExplorer.uniquefy(filteredCompounds);
		for(IAtomContainer a : uniqueFilteredCompounds.atomContainers()){
//			System.out.println(this.ebt.smiGen.create(a) + " - " + a.getProperty("Major Isotope Mass"));
			pathways.addAtomContainer(findPathway(startingCompoundStandardized, a, compoundToParents, inchikeyToContainers, annotate));
		}		
				
		return pathways;
	}



	public IAtomContainer findPathway(IAtomContainer startingCompound, IAtomContainer leafCompound,
			LinkedHashMap<String, IAtomContainerSet> compoundToParents, LinkedHashMap<String,IAtomContainer> inchikeyToContainers, boolean annotate
			
			
			) throws CloneNotSupportedException, CDKException{
		
		LinkedHashMap<Object, Object> props = new LinkedHashMap<Object, Object>();
		ArrayList<String> l = new ArrayList<String>();
		ArrayList<String> transversal = new ArrayList<String>();
	
//		IAtomContainerSet traversedNodes = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);;
		
		LinkedHashMap<String, IAtomContainer> traversedNodes = new LinkedHashMap<String, IAtomContainer>();
		
		props.put("InChI", (String) leafCompound.getProperty("InChI"));
		props.put("InChIKey", (String) leafCompound.getProperty("InChIKey"));
		
		props.put("Molecular formula", (String) leafCompound.getProperty("Molecular formula"));
		props.put("Major Isotope Mass", (String) leafCompound.getProperty("Major Isotope Mass"));
		
		if(annotate){
			if(leafCompound.getProperty("Synonyms") != null){
				String[] synonyms = ((String) leafCompound.getProperty("Synonyms")).split("\n");
//				for(int s = 0; s < synonyms.length; s++){
//					System.out.println(synonyms[s]);
//				}
				props.put(CDKConstants.TITLE, synonyms[0]);
			}
		}
		
		
		int step = 0;
		
		if(((String)leafCompound.getProperty("Precursor InChIKey")).contentEquals((String)startingCompound.getProperty("InChIKey"))){
			String annotation = "";
//			System.out.println("Reaction: " + leafCompound.getProperty("Reaction").toString());
//			System.out.println("Precursor InChI: " + leafCompound.getProperty("Precursor InChI").toString());
//			annotation += "Biotransformation_" + String.valueOf(step+1) + "\n";
			annotation += "Reaction Type: " + leafCompound.getProperty("Reaction") + 
			" (" + leafCompound.getProperty("Reaction ID") + ")\n";
			String[] enz = ((String) leafCompound.getProperty("Enzyme(s)")).split("\n");
			annotation += "Enzyme(s): " + StringUtils.join(enz , "; ") + "\n";
			annotation += "Precursor InChI: " + leafCompound.getProperty("Precursor InChI") + "\n";
			annotation += "Precursor InChIKey: " + leafCompound.getProperty("Precursor InChIKey") + "\n";					
			transversal.add(annotation);
			traversedNodes.put((String) leafCompound.getProperty("InChIKey"), leafCompound);			
		}
//		else if (ChemStructureExplorer.atomContainerInclusionHolds(compoundToParents.get((String)leafCompound.getProperty("InChIKey")), startingCompound)){
//					
//		
//		
//		}
		else{
			
//			System.out.println("STARTING COMPOUND: " + (String)startingCompound.getProperty("InChIKey"));
			
			IAtomContainer currentNode = leafCompound;
			while(! ((String)currentNode.getProperty("InChIKey")).contentEquals((String)startingCompound.getProperty("InChIKey"))){
				IAtomContainerSet allParentsOfCurrentNode = compoundToParents.get((String)currentNode.getProperty("InChIKey"));
				
				String annotation = "";
				for (IAtomContainer a : allParentsOfCurrentNode.atomContainers()){
					IAtomContainer atc = inchikeyToContainers.get(a.getProperty("Precursor InChIKey"));
//					boolean moveOn = false;
//					for(IAtomContainer atc : allParentsOfCurrentNode.get(i).getSubstrates().atomContainers() ){
						if(!traversedNodes.containsKey((String)atc.getProperty("InChIKey"))){
							
							
//							System.out.println((String)atc.getProperty("InChIKey"));
//							System.out.println("Reaction: " + a.getProperty("Reaction").toString());
//							System.out.println("Precursor InChI: " + a.getProperty("Precursor InChI").toString());
//							annotation += "Biotransformation_" + String.valueOf(step+1) + "\n";
							annotation += "Reaction Type: " + a.getProperty("Reaction") + 
							" (" + currentNode.getProperty("Reaction ID") + ")\n";
							String[] enz = ((String) a.getProperty("Enzyme(s)")).split("\n");
							annotation += "Enzyme(s): " + StringUtils.join(enz , "; ") + "\n";
							annotation += "Precursor InChI: " + a.getProperty("Precursor InChI") + "\n";
							annotation += "Precursor InChIKey: " + a.getProperty("Precursor InChIKey") + "\n";					
							transversal.add(annotation);
							traversedNodes.put((String) atc.getProperty("InChIKey"), atc);
							step++;
							currentNode = atc;
//							moveOn = true;
							break;
						}
//						if(moveOn == true){
//							break;
//						}
//					}
				}
//				System.out.println(((String)currentNode.getProperty("InChIKey")).contentEquals((String)startingCompound.getProperty("InChIKey")));
			}			
		}
		
		String p = "";
		for(int z = 0; z < transversal.size(); z++){
//			p += "Biotransformation_" + String.valueOf(transversal.size() - z) + "\n";
//			p += transversal.get(z) + "\n";
			p += "Biotransformation_" + String.valueOf(z + 1) + "\n";
			p += transversal.get(transversal.size() - z - 1) + "\n";
			
//			p += "Biotransformation_" + String.valueOf(z+1) + "\n";
//			p += transversal.get(z) + "\n";
		}

		
		props.put("Pathway", p);
		
		
		
//		currentNode = inchikeyToContainers.get((String) currentParent.getProperty("Precursor InChIKey")); 
					
//		System.out.println("currentParent: " + (String)currentParent.getProperty("InChIKey") + "\n"); 
		
		
//		System.out.println("(String)currentParent.getProperty(InChIKey): " + (String)currentParent.getProperty("InChIKey"));
//		System.out.println("(String)startingCompound.getProperty(InChIKey): " + (String)startingCompound.getProperty("InChIKey"));
//		traversedNodes.put(currentParent.getProperty("InChIKey"), currentParent);
		
		
		
//		while(currentParent != null && (String)currentParent.getProperty("InChIKey") != (String)startingCompound.getProperty("InChIKey")){
//			System.out.println((String)currentParent.getProperty("InChIKey") != (String)startingCompound.getProperty("InChIKey"));
//			step++;
//			String annotation = "";
//			IAtomContainerSet allParentsOfCurrentCompound = (IAtomContainerSet) productsToBio.get((String)currentParent.getProperty("InChIKey"));
//
//			for(IAtomContainer atc : allParents)
//			
//			
//			annotation += "Reaction Type: " +  (String) currentParent.getProperty("Reaction") + 
//					" (" + currentParent.getProperty("Reaction ID") + ")\n";
//			String[] enz = ((String) currentParent.getProperty("Enzyme(s)")).split("\n");
//			annotation += "Enzyme(s): " + StringUtils.join(enz , "; ") + "\n";
//			annotation += "Precursor InChI: " + currentParent.getProperty("Precursor InChI") + "\n";
//			annotation += "Precursor InChIKey: " + currentParent.getProperty("Precursor InChIKey") + "\n";
//			
////			props.put("Biotransformation_" + String.valueOf(step), annotation);
//			l.add(annotation);
//			currentParent = inchikeyToContainers.get((String) currentParent.getProperty("Precursor InChIKey")); 
//			System.out.println(annotation);			
////			System.out.println("currentParent: " + (String)currentParent.getProperty("InChIKey") + "\n"); 
//		}
		
//		for(int a = 0; a < l.size(); a++){
//			props.put("Biotransformation_" + String.valueOf(a+1), l.get(l.size()-a-1));
//		}
		
		IAtomContainer path = leafCompound.clone();
		path.setProperties(props);
		
		return path;
	}
	

	public IAtomContainerSet findMetabolitesFromSequence(IAtomContainer startingCompound, BiotransformerSequence bseq, ArrayList<String>mass_formulas, double massTolerance, boolean annotate, FinderOption opt, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> allBiotransformations = new ArrayList<Biotransformation>();
		
		
		ArrayList<String> remainingMassFormulas = (ArrayList<String>) mass_formulas.clone();
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		LinkedHashMap<String, ArrayList<Biotransformation>> substratesToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, ArrayList<Biotransformation>>productsToBio = new LinkedHashMap<String, ArrayList<Biotransformation>>();
		LinkedHashMap<String, IAtomContainerSet> massFormulaToMolecules = new LinkedHashMap<String, IAtomContainerSet>();		
		LinkedHashMap<String, IAtomContainer> inchikeyToContainers = new LinkedHashMap<String, IAtomContainer>();
		LinkedHashMap<String, IAtomContainerSet> compoundToParents = new LinkedHashMap<String, IAtomContainerSet>();
		
		
		
		
		IAtomContainer startingCompoundStandardized = this.hsbt.standardizeMoleculeWithCopy(startingCompound);
			
		ChemStructureExplorer.addInChIandKey(startingCompoundStandardized);
		inchikeyToContainers.put((String) startingCompoundStandardized.getProperty("InChIKey"), startingCompoundStandardized);
		
		IAtomContainerSet pathways = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

		
		for(int k=0; k < mass_formulas.size(); k++){
			 massFormulaToMolecules.put(mass_formulas.get(k), DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class));
		}
		
		IAtomContainerSet currentCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentCompounds.addAtomContainer(startingCompound);
	

		for (BiotransformerSequenceStep btStep :bseq.sequence) {
			
			if (remainingMassFormulas.size() == 0) {
				break;
			}
			else {

//				int counter = 0;
				int nrOfSteps = btStep.nOfIterations;
				while(nrOfSteps >0 && remainingMassFormulas.size()>0){	
//					for(String x : remainingMassFormulas) {
//						System.out.print(x + "\t");
//					}
//					System.out.print("\n");
					ArrayList<BiotransformerSequenceStep> currentsteps = new ArrayList<BiotransformerSequenceStep>();
					currentsteps.add(new BiotransformerSequenceStep(btStep.btype, 1, bseq.scoreThreshold));
					
					BiotransformerSequence currentSeq = new BiotransformerSequence(currentsteps , bseq.scoreThreshold);				
//					System.out.println("Current number of compounds: " + currentCompounds.getAtomContainerCount());
					ArrayList<Biotransformation> currentBioT = currentSeq.runSequence(currentCompounds, bseq.scoreThreshold, cyp450Mode);				
					System.out.println("Current biotransformations: " + currentBioT.size());
					IAtomContainerSet extractedCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//					extractedCompounds = this.ubt.extractProductsFromBiotransformationsWithTransformationData(currentBioT, annotate);
					if (btStep.btype != bType.ENV) {
						extractedCompounds = this.hsbt.extractProductsFromBiotransformationsWithTransformationData(currentBioT, annotate);
					}
					else {
						extractedCompounds = this.ebt.extractProductsFromBiotransformationsWithTransformationData(currentBioT, annotate);
					}
					
					allBiotransformations.addAll(currentBioT);
					currentCompounds.removeAllAtomContainers();
					currentCompounds.add(ChemStructureExplorer.uniquefy(extractedCompounds));
					
					for(IAtomContainer a : currentCompounds.atomContainers()){
						inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
					}
					
					for(IAtomContainer a : currentCompounds.atomContainers()){
//						inchikeyToContainers.put((String) a.getProperty("InChIKey"), a);
						String precursorInChIKey = (String) a.getProperty("Precursor InChIKey");
//						System.out.println("precursorInChIKey: " + precursorInChIKey);
//						System.out.print((String) a.getProperty("InChIKey"));
//						System.out.println(" >> " + (String)a.getProperty("Precursor InChIKey"));
//						System.out.println(inchikeyToContainers.containsKey(precursorInChIKey));
						
//						if(inchikeyToContainers.containsKey(precursorInChIKey)){
							if(compoundToParents.containsKey((String) a.getProperty("InChIKey"))){
//								compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
								compoundToParents.get((String) a.getProperty("InChIKey")).addAtomContainer(a);
							}
							else {
								IAtomContainerSet s = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//								s.addAtomContainer(inchikeyToContainers.get(precursorInChIKey));
//								compoundToParents.put((String) a.getProperty("InChIKey"), s);	
								s.addAtomContainer(a);
								compoundToParents.put((String) a.getProperty("InChIKey"), s);
							}
//						}
//						System.out.println(compoundToParents.get((String) a.getProperty("InChIKey")) == null);
					}			
					
					
					for(int k=0; k < mass_formulas.size(); k++){				
						if(opt == FinderOption.FORMULA){
							IMolecularFormula formula1 = MolecularFormulaManipulator.getMolecularFormula(mass_formulas.get(k), this.builder);
							for(IAtomContainer a : extractedCompounds.atomContainers()){ 
								if(a.getProperty("Molecular formula") == null){
									a.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(a));
								}
								System.out.println("A: " + MolecularFormulaManipulator.getString(formula1));
								System.out.println("B: " + a.getProperty("Molecular formula").toString());
								if(MolecularFormulaManipulator.getString(formula1).contentEquals(a.getProperty("Molecular formula").toString().trim())){
									System.out.println(a.getProperty("Molecular formula").toString().trim());
									massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
									filteredCompounds.addAtomContainer(a);
//									break;
								}
							}					
						}
						else if(opt == FinderOption.MASS){
							Double mass = Double.valueOf(mass_formulas.get(k));
							for(IAtomContainer a : extractedCompounds.atomContainers()){ 
								if( Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass) <= massTolerance){
									massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
									filteredCompounds.addAtomContainer(a);
//									break;
								}
							}					
						}
						else if(opt == FinderOption.MASSFORMULA){
							String[] mf = mass_formulas.get(k).split(":");
							IMolecularFormula formula2 = MolecularFormulaManipulator.getMolecularFormula(mf[1], this.builder);
							Double mass2 = Double.valueOf(mf[0]);
							for(IAtomContainer a : extractedCompounds.atomContainers()){ 
								if( (MolecularFormulaManipulator.getString(formula2).contentEquals(mass_formulas.get(k))) && 
										(Math.abs(Double.valueOf((String) a.getProperty("Major Isotope Mass")) - mass2) <= massTolerance)){
									massFormulaToMolecules.get(mass_formulas.get(k)).addAtomContainer(a);
									filteredCompounds.addAtomContainer(a);
//									break;
								}
							}					
						}
						if(massFormulaToMolecules.get(mass_formulas.get(k)).getAtomContainerCount() > 0){
							remainingMassFormulas.remove(mass_formulas.get(k));
						}
						
					}			
					nrOfSteps--;
//					counter++;
//					System.err.println("Remaining masses: " + remainingMasses.size() + "\n");
//					System.err.println("\n\n\n\n\nNumber of possibilities: " + filteredCompounds.getAtomContainerCount() );
				}				
			}

		}
		
		for(int i = 0; i < allBiotransformations.size(); i++){
			Biotransformation bt = allBiotransformations.get(i);
			for(IAtomContainer ac : bt.getSubstrates().atomContainers()){
				String ikey = ac.getProperty("InChIKey");
				if(ikey != null){
					if(substratesToBio.containsKey(ikey)){
						substratesToBio.get(ikey).add(bt);
					}else{
						substratesToBio.put(ikey,new ArrayList<Biotransformation>());
						substratesToBio.get(ikey).add(bt);
					}
				}
			}

			for(IAtomContainer at : bt.getProducts().atomContainers()){
				String ikey = at.getProperty("InChIKey");
				if(ikey != null){
					if(productsToBio.containsKey(ikey)){
						productsToBio.get(ikey).add(bt);
					}else{
						productsToBio.put(ikey,new ArrayList<Biotransformation>());
						productsToBio.get(ikey).add(bt);
					}
				}
			}
		}
		
		IAtomContainerSet uniqueFilteredCompounds = ChemStructureExplorer.uniquefy(filteredCompounds);
		for(IAtomContainer a : uniqueFilteredCompounds.atomContainers()){
//			System.out.println(this.ebt.smiGen.create(a) + " - " + a.getProperty("Major Isotope Mass"));
			pathways.addAtomContainer(findPathway(startingCompoundStandardized, a, compoundToParents, inchikeyToContainers, annotate));
		}		
				
		return pathways;
	}
	

	
	
	
	
	public boolean compoundsHasToleratedMass(IAtomContainer atc, ArrayList<Double> masses, double tolerance){
		boolean pass = false;
		double mass = atc.getProperty("Major Isotope Mass");
//		if(mass == null){
//			// FILL THIS IN. calculate the mass and update the properties
//		}
		
		for(int i = 0; i < masses.size(); i++){
			if(Math.abs(mass - masses.get(i)) <= tolerance){
				return true;
			}
		}		
		return pass;
	}
}
