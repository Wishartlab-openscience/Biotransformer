/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.BufferedReader;
//import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
//import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
//import java.io.OutputStreamWriter;
//import java.io.Writer;
//import java.nio.file.Path;
//import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.apache.commons.io.FilenameUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
//import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

//import biotransformer.transformation.Biotransformation;
//import biotransformer.transformation.MetabolicReaction;
 


public class FileUtilities {

	public FileUtilities() {
		// TODO Auto-generated constructor stub
//		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
	
	}
	
	
	public static IAtomContainerSet parseSdf(String sdfFileName) throws IOException {
		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(sdfFileName), bldr);
		
		while (sdfr.hasNext()){
			IAtomContainer mol = sdfr.next();
			containers.addAtomContainer(mol);	
		}
		
		sdfr.close();
		return containers;
		
	}
	
	
	public static IAtomContainerSet parseSdfAndAddTitles(String sdfFileName, InChIGeneratorFactory igf) throws CDKException, IOException {
		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(sdfFileName), bldr);
		
		while (sdfr.hasNext()){
			IAtomContainer molecule = sdfr.next();
			String identifier = molecule.getProperty(CDKConstants.TITLE);
			if(identifier == null){
				identifier = molecule.getProperty("Name");
				if(identifier == null){
					identifier = molecule.getProperty("$MolName"); 
					if(identifier == null){
						identifier = molecule.getProperty("InChIKey");
						if(identifier == null){
							identifier = igf.getInChIGenerator(molecule).getInchiKey();
						}
					}

				}
				molecule.setProperty(CDKConstants.TITLE, identifier);
			}
			
			containers.addAtomContainer(molecule);	
		}
		
		sdfr.close();
		return containers;
		
	}
	
	
	
	public static int countUniqueCompounds(String sdfFileName) throws CDKException, IOException{
		int count = 0;
		
		IAtomContainerSet containers = parseSdf(sdfFileName);
		
		LinkedHashMap<String, IAtomContainer> lmh  = new LinkedHashMap<String, IAtomContainer>();
		
		for (IAtomContainer atc : containers.atomContainers()){
			String ikey = atc.getProperty("InChIKey");
			if(ikey == null){
				InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
				ikey = gen.getInchiKey();			
			}
			
			if(! lmh.containsKey(ikey)){
				lmh.put(ikey, atc);
				count++;
			}
		}
		
		
		return count;
	}

	
	public static void divideSdfFile(String sdfFileName, int limit) throws CDKException, IOException{
		IAtomContainerSet containers = parseSdf(sdfFileName);
		int part_nr = 1;
		int at_nr=0;
		
		if(containers.getAtomContainerCount()>limit){
			IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			System.out.println(sdfFileName);
			System.out.println(FilenameUtils.getBaseName(sdfFileName));
			System.out.println(FilenameUtils.getFullPathNoEndSeparator(sdfFileName));
			
			for(IAtomContainer atc : containers.atomContainers()){
				at_nr++;
				System.out.println(at_nr);
				molecules.addAtomContainer(atc);
				if(at_nr % limit == 0){
					SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(sdfFileName) + 
							"/" + FilenameUtils.getBaseName(sdfFileName)  + "_part_"+ part_nr + ".sdf" ));
					for(IAtomContainer a : molecules.atomContainers()){
						sdfWriter.write(a);						
					}
					sdfWriter.close();
					part_nr++;
					molecules.removeAllAtomContainers();
				}
			}
		
		
		}
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
		
		while((line = bRead.readLine()) !=null && !line.contains("SMILES")){
			counter++;
			System.out.println(counter);
			String[] sline = line.split("\t");
			
			if(!sline[1].contentEquals("NULL")){
				
				System.out.println(counter + " " + sline[0]);
				IAtomContainer atc= smiParser.parseSmiles(sline[1]);
				AtomContainerManipulator.suppressHydrogens(atc);
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
				atc.setProperty(CDKConstants.TITLE, sline[0]);
				atc.setProperty("Name", sline[0]);
				
//				if(gen.getInchiKey().trim().contentEquals(sline[3].trim())){
					atc.setProperty("InChIKey", gen.getInchiKey());
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
				
				atc.setProperty("SoMs", "");
				atc.setProperty("References","");
				sdfWriter.write( ChemStructureManipulator.preprocessContainer(atc) );
			}
		}
		
		sdfWriter.close();
		bRead.close();
	}

	
	
	public static void saveAtomContainerSetToCSV(IAtomContainerSet products, String outputFileName) throws Exception{

		LinkedHashMap<Object, Object> properties = new LinkedHashMap<Object, Object>();
		
		if(products.getAtomContainerCount() > 0){
			
		}
		
		try{
			
			if(products.getAtomContainerCount() > 0){				
				
				ArrayList<String> header = new ArrayList<String>();	
				for( Object prop : products.getAtomContainer(0).getProperties().keySet()){
					header.add(String.valueOf(prop));
				}
				FileWriter fWriter = new FileWriter(outputFileName);
				CSVPrinter csvPrinter 	= new CSVPrinter(fWriter, CSVFormat.DEFAULT);
				csvPrinter.printRecord(header);
				for(IAtomContainer atc : products.atomContainers()){
					ArrayList<String> props = new ArrayList<String>();
					for (int k = 0; k < header.size(); k++){
						props.add((String) atc.getProperty(header.get(k)));
					}
					csvPrinter.printRecord(props);
				}				
//				csvPrinter.flush();
				csvPrinter.close();
				System.out.println("The results were saved to the following file: " + outputFileName + "\n");
				
			}
			else{
				System.out.println("The number of metabolites to save is " + (products.getAtomContainerCount()));
			}
	 
		}catch (Exception e) {
			
			e.printStackTrace();

		}	
		
	}
	
	public static void saveAtomContainerSetToSDF(IAtomContainerSet containers, String outputFileName) throws CDKException, IOException{
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(containers);
		sdfWriter.close();
	}

	public static String saveAtomContainersToString(IAtomContainerSet containers) throws CDKException, IOException{
		ByteArrayOutputStream containers_to_s = new ByteArrayOutputStream();
		SDFWriter sdfWriter = new SDFWriter(containers_to_s);		
		sdfWriter.write(containers);
		sdfWriter.close();
		
		return containers_to_s.toString();
		
	}
	

	public static String saveAtomContainerToString(IAtomContainer container) throws CDKException, IOException{
		ByteArrayOutputStream containers_to_s = new ByteArrayOutputStream();
		SDFWriter sdfWriter = new SDFWriter(containers_to_s);		
		sdfWriter.write(container);
		sdfWriter.close();
		
		return containers_to_s.toString();
		
	}	
}
