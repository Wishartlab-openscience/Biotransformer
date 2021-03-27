/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.fingerprint;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.isomorphism.matchers.smarts.SmartsMatchers;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.hp.hpl.jena.util.FileUtils;

import ambit2.smarts.query.SmartsPatternCDK;
import biotransformer.utils.ChemStructureExplorer;
import phase2filter.utils.ChemStructureManipulator;

public class ChemStructureFingerprinter {

//	public ChemStructureFingerprinter() {
//		// TODO Auto-generated constructor stub
//		
//	}

	public List<List<Integer>> getMatchingAtoms(IAtomContainer molecule,
			SMARTSQueryTool smartsPattern) throws Exception {
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		boolean occurrences = smartsPattern.matches(molecule);
		List<List<Integer>> indexes = smartsPattern.getMatchingAtoms();
//		for (int i = 0; i < indexes.size(); i++) {
//			List<Integer> currentAtomIndexes = indexes.get(i);
//			System.out.println("\nMatch nr." + (i + 1) + ": " + (currentAtomIndexes));
//			for (int j = 0; j < currentAtomIndexes.size(); j++) {
//				System.out.println(currentAtomIndexes.get(j) + ": "
//						+ molecule.getAtom((int) currentAtomIndexes.get(j)));
//			}
//		}
		return indexes;
	}

	
	public static LinkedHashMap<String, String> getMiniFingerprintPatterns(){
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl_aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("ketone","[#6]-[#6;X3](-[#6])=[O;X1]");

		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("alcohol","[OX2H][CX4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("aromatic_alcohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("phenol",
				"[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");		
		
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]");
		queries.put("carboxylic_ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		
		queries.put("amine", "[NX3+0,NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]");		
		queries.put("quaternary_ammonium",
				"[NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("n_acyl_aromatic_alpha_amino_acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put(
				"n_acyl_aliphatic_alpha_amino_acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put("nitrile", "C#N");
		
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
				
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");

		queries.put(
				"organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");

		queries.put(
				"steroid",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put(
				"coenzyme_a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty_acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put(
				"_2_N_linked_ribose_deriv",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");

		queries.put("glycerol_3_phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");
		queries.put("glycerolipid", 
				"[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]");
		queries.put("glycerophospholipid", 
				"[#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("glycerophosphonolipid", 
				"[#6;A]P([#8;A;X2H1,X1-])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("sphingolipid", 
				"[$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8]),$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X4]=[#6;A;H1X4]-[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8])]");
		
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("flavonoid","[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]1=,:[#6;R1][#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8]1),"
				+ "$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]1[#6;R1]=,:[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8]1),"
				+ "$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8+]=,:1),"
				+ "$([#6;R1]-,=1-[#6;R1]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#8;R1]-[#6;R1]-,=1-[#6;R1]-1=[#6;R1]-[#6]=[#6]-[#6]=[#6;R1]-1)]");
		
		queries.put("isoflavonoid","[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]=,:1[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1][#6;R1]=,:1),"
				+ "$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]=,:1[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1][#6;R1]=,:1)]");
		
		queries.put("neoflavonoid","[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]1[#6;R1]=,:[#6;R1][#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2),"
				+ "$([#6;R1]1[#6;R1][#6;R1]([#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1]1)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)]");		
		
		queries.put("glucuronic_acid","[#8]!@-[#6;A;H1X4R1]1[#8;R1][#6;A;H1X4R1]([#6;A;H1X4R1](!@-[#8])[#6;A;H1X4R1](!@-[#8])[#6;A;H1X4R1]1!@-[#8])[#6;X3]([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("sulfuric_acid",
				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("glutathione", "[#7;H2X3]-[#6;H1X3]([#6;A;H2X4][#6;A;H2X4][#6;X3](=O)[#7;A;H1X3][#6;H1X3]([#6;A;H2X4][#16;A;X1-,X2])-[#6;X3](=O)[#7;A;H1X3][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=O)-[#6;X3]([#8;A;X2H1,X1-])=O");
		queries.put("methoxy_group","[#6]!@-[#8;X2]!@-[#6;A;H3X4]");
		queries.put("glycine","[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]");
		queries.put("n_hydroxylarylamine_o_acetate", "[#6;H3X4]-[#6;X3](=[O;X1])-[#8;X2][#7;A;H1X3]!@-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("n_hydroxylamine_n_acetate","[H]C([H])([H])[#6](=O)-[#8][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6]");
		
		return queries;
	}
	
	
	/**
	 * 
	 * @return : A HashMap with the SMARTS expressions for functional groups and
	 *         structural patterns
	 * @throws Exception
	 */
	public static LinkedHashMap<String, String> getFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic", "[*;a]");
		queries.put("sulfuric acid",
				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide", "[F,Cl,Br,I]");
		queries.put("acyl chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put(
				"organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[c;R1]1[c;R1][c;R1][o;R1][c;R1]1");
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put(
				"hydroxylamine",
				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		// queries.put("hydroxylamine","[NX3H2,$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][OX2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
		queries.put("phenol",
				"[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");
		queries.put("primary carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("tertiary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl halide", "[F,Cl,Br,I][c]");
		queries.put("C ONS bond", "[#6]~[#7,#8,#16]");
		queries.put(
				"13-dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary aliph ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary arom ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]-1=[O;X1]");
		queries.put("p-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1](=O)-[#6;R1]=[#6;R1]-1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																											// this
																											// class
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																						// this
																						// class
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put(
				"thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put(
				"amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put(
				"n-hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7H1]-[#7H1]-[#6;a]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1-(acyloxy)ethyl ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put(
				"(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine_monocyclic", "[#6;R1]1=,:[#6;R1][#7;R1]=,:[#6;R1][#7;R1]=,:[#6;R1]1");
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[SX2][CX2]#[NX1]");
		queries.put("alpha,beta-unsaturated system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
		queries.put(
				"CH-acidic",
				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put(
				"CH-acidic strong",
				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N-aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary aliph amine",
				"[NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary aliph amine",
				"[NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary aliph amine",
				"[NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("primary arom amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary arom amine",
				"[NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary arom amine",
				"[NX3H0+0,NX4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("peroxo group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put(
				"thioamide S-oxide derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put(
				"thioamide SS-dioxide derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n-nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c-nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic acid ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put(
				"sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put(
				"steroid",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene_monocyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#16;R1][#6;R1]=,:1");
		queries.put("thiophene", "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1");		
		queries.put("thiophene s-oxide", "[#8;X1-][S+]-,:1-,:[#6]=,:[#6][#6]=,:[#6]-,:1");
		queries.put("1,3-thiazole", "[#6]1=,:[#6][#16][#6]=,:[#7]1"); // SMARTCyp - a 2D-method for Prediction of Cytochrome P450 Mediated Drug Metabolism_Supplement
		queries.put("1,2-thiazole", "[#6]=,:1[#6]=,:[#7][#16][#6]=,:1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "[#6]1=,:[#6][#6]=,:[#7][#7]=,:[#6]1");
		queries.put("1,3-dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put(
				"biphenyl",
				"[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put(
				"1,2-aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur compound", "[#6]~[#16]");
		queries.put("secondary aliphatic/aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put(
				"propargyl-type 13-dipolar organic compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put(
				"diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put(
				"phenol ether",
				"[#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("heteroaromatic", "[a;!c]"); // 107 bits
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put(
				"hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put(
				"hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put(
				"steroid",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("phosphate monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3-acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n-acyl-aromatic alpha-amino-acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put(
				"n-acyl-aliphatic alpha-amino-acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[$([#6]-[#7;X3](-[#1,#6])-[#7;X3+](-[#8;X1-])=[O;X1]),$([#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3])]");
		queries.put("1,2-oxazole", "[#8]-1-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("1,3-oxazole", "[#8]-1-[#6]=[#6]-[#7]=[#6]-1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("glycerol-3-phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");

		queries.put(
				"coenzyme a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put(
				"2-N-linked ribose deriv.",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic aocohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("alpha-amino acid or derivatives", 
				"[#7;A][#6;X4]-[#6;X3]([!#1!#6])=[O;X1]");  //modified from ClasyFire's alpha-amino-acid-erivative-1 by adding X4 on the C2 carbon.
		queries.put("alpha-amino acid", 
				"[#7;A;X3,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X4]-[#6;X3]([#8;A;X2H1,X1-])=[O;X1]"); //modified from ClasyFire's alpha-amino-acid-1 by adding X4 on the C2 carbon.
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");

		queries.put("tryptamine", 
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6;R1]-1=[#6;R1]-[#7;R1]-[#6;R2]-2=[#6;R2]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
				
		queries.put("flavonol", 
				"[H][#8]-[c;R1]1[#6;R1](=O)c2-,:cc-,:cc-,:c2[#8;R1][c;R1]1-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1]([#8;A;H1X2])[c;R1]-,:1");
		queries.put("flavone", 
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:1)-[#6;R1]=,:1[#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1](=O)[#6;R1]=,:1-[!#8]"); // PMID: 15914008
		queries.put("coumarin", 
				"[$(O=[#6]-1-[#8]-c2ccccc2-[#6]=[#6]-1),$(O=[#6]-1-[#8]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]=[#6]-1)]");

		return queries;
	}
	
	public static LinkedHashMap<String, String> getRINFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic", "[*;a]");
		queries.put("sulfuric acid",
				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide", "[F,Cl,Br,I]");
		queries.put("acyl chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put(
				"organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[c;R1]1[c;R1][c;R1][o;R1][c;R1]1");
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put(
				"hydroxylamine", 
				"[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
//				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		// queries.put("hydroxylamine","[NX3H2,$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][OX2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
		queries.put("phenol",
				"[#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("primary carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("tertiary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl halide", "[F,Cl,Br,I][c]");
		queries.put("C ONS bond", "[#6]~[#7,#8,#16]");
		queries.put(
				"13-dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary aliph ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary arom ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o-quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]1=[O;X1]");
		queries.put("p-quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1](=O)[#6;R1]=,:[#6;R1]1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																											// this
																											// class
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																						// this
																						// class
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put(
				"thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put(
				"amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put(
				"n-hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7;X3](-[#6,#1])-[#7;X3](-[#6,#1])-[#6,#1]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1-(acyloxy)ethyl ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put(
				"(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine_monocyclic", "[#6;R1]1=,:[#6;R1][#7;R1]=,:[#6;R1][#7;R1]=,:[#6;R1]1");
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[SX2][CX2]#[NX1]");
		queries.put("alpha_beta-unsaturated_system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
		queries.put(
				"CH-acidic",
				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put(
				"CH-acidic strong",
				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N-aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary aliph amine",
				"[NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary aliph amine",
				"[NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary aliph amine",
				"[NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("primary arom amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary arom amine",
				"[NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary arom amine",
				"[NX3H0+0,NX4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("peroxo group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put(
				"thioamide S-oxide derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put(
				"thioamide SS-dioxide derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n-nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c-nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic acid ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put(
				"sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put(
				"steroid",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene_monocyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#16;R1][#6;R1]=,:1");
		queries.put("thiophene", "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1");		
		queries.put("thiophene s-oxide", "[#8;X1-][S+]-,:1-,:[#6]=,:[#6][#6]=,:[#6]-,:1");
		queries.put("1,3-thiazole", "[#6]1=,:[#6][#16][#6]=,:[#7]1"); // SMARTCyp - a 2D-method for Prediction of Cytochrome P450 Mediated Drug Metabolism_Supplement
		queries.put("1,2-thiazole", "[#6]=,:1[#6]=,:[#7][#16][#6]=,:1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "[#6]1=,:[#6][#6]=,:[#7][#7]=,:[#6]1");
		queries.put("13-dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put(
				"biphenyl",
				"[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put(
				"12-aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur compound", "[#6]~[#16]");
		queries.put("secondary aliphatic/aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put(
				"propargyl-type 13-dipolar organic compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put(
				"diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put(
				"phenol ether",
				"[$([#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");
		queries.put("heteroaromatic", "[a;!c]"); // 107 bits
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put(
				"hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put(
				"hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put(
				"steroid",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("phosphate monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3-acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n-acyl-aromatic alpha-amino-acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put(
				"n-acyl-aliphatic alpha-amino-acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[$([#6]-[#7;X3](-[#1,#6])-[#7;X3+](-[#8;X1-])=[O;X1]),$([#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3])]");
		queries.put("12-oxazole", "[#6]=,:1[#6]=,:[#7][#8][#6]=,:1");
		queries.put("13-oxazole", "[#6]1=,:[#6][#8][#6]=,:[#7]1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put(
				"coenzyme a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put(
				"2-N-linked ribose deriv.",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic aocohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("glycerol-3-phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");

//		Newly added ones
		queries.put("biguanide", 
				"[$([#7]!@-[#6](=[#7])!@-[#7]!@-[#6](!@-[#7])=[#7]),$([#7]!@-[#6](!@-[#7])=[#7]!@-[#6](!@-[#7])=[#7])]");
		queries.put("glycerolipid", 
				"[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]");
		queries.put("glycerophospholipid", 
				"[#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("glycerophosphonolipid", 
				"[#6;A]P([#8;A;X2H1,X1-])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("sphingolipid", 
				"[$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8]),$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X4]=[#6;A;H1X4]-[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8])]");
		queries.put("glycero-3-dithiophosphocholine", 
				"[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#8;X2]!@-P(!@-[#16])(!@=[S;X1])!@-[#8]!@-[#6;A;H2X4]!@-[#6;A;H1X4](!@-[#8])!@-[#6;A;H2X4]!@-[#8]");
		queries.put("glycero-3-thiophosphocholine", 
				"[#6;A;H3X4][N;R0+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#8]P([#8;A;X2H1,X1-])(=[S;X1])[#8;R0][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("bisphosphatidic acid", 
				"[#6]!@-[#6]!@-[#6]!@-[#6](!@=O)!@-[#8]!@-[#6]!@-[#6](!@-[#6]!@-[#8]!@-P(!@-[#8])(!@=O)!@-[#8]!@-[#6]!@-[#6](!@-[#6]!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6])!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6])!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6]");
		queries.put("lysobisphosphatidic acid", 
				"[$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0]),$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#8]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0]),$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#8]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])]");
		queries.put("phosphatidylethanol", 
				"[#6;H3X4R0]-[#6;H2X4R0]-[#8;R0]P([!#1!#6;$([OX2H1]),$([O-])])(=O)[#8;R0]-[#6;H2X4R0]-[#6;H1X4R0](-[#6;H2X4R0]-[#8;R0]-[#6;R0](-[#6,#1])=[O;R0])-[#8;R0]-[#6;R0](-[#6,#1])=[O;R0]");		
//		queries.put("phosphoglycosphingolipid", 
//				"");
		queries.put("semilysobisphosphatidic acid", 
				"[H][#8;R0]-[#6;R0](-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])-[#6;R0]-[#8]P([#8])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])-[#8]-[#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]");
		queries.put("3-benzazepine", 
				"N1C=CC2=CC=CC=C2C=C1");
		queries.put("dibenzo_p_dioxin", 
				"[#8]-1-c2ccccc2-[#8]-c2ccccc-12");		
		queries.put("13-oxazin-2-one",
				"[$(O=C1OC=CC=N1),$(O=C1NC=CCO1)]"); // dx.doi.org/10.1021/ml500297n | ACS Med. Chem. Lett. 2014, 5, 1156−1161
		queries.put("13-oxazine", 
				"[$(C1NC=C=CO1),$([#6]-1-[#8]-[#6]=[#7]-[#6]=[#6]-1),$([#6]-1-[#6]=[#6]-[#8]-[#6]=[#7]-1),$([#6]-1-[#8]-[#6]=[#6]-[#6]=[#7]-1),$([#6]-1-[#8]-[#6]-[#7]=C=[#6]-1)]"); // dx.doi.org/10.1021/ml500297n | ACS Med. Chem. Lett. 2014, 5, 1156−1161
		queries.put("benzo-14-diazepine", 
				"[$([#6]-,=1-,=[#6]-[#7]=[#6]-c2ccccc2-[#7]-,=1),$([#6]-,=1-,=[#6]-[#7]=[#6]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#7]-,=1)]");		
		queries.put("1-benzofuran", 
				"[$(c1cc2ccccc2o1),$(O1C=CC2=CC=CC=C12)]");
		queries.put("benzimidazole", 
				"[$(c1nc2ccccc2n1),$(N1C=NC2=CC=CC=C12)]");		
		queries.put("sulfonylurea", 
				"[$([#7;X3]!@-[#6](=O)!@-[#7]=S(=O)=O),$([#6]S(=O)(=O)[#7]-[#6](-[#7;X3])=O),$([#6]S(=O)(=O)-[#7]=[#6](/[#7;X3])-[#8]),$([#6]S(=O)(=O)[#7]-[#6](-[#8])=[#7;X2])]");	
		
		queries.put("phosphate diester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2]-[#6]");
		queries.put("alpha-amino acid or derivatives", 
				"[#7;A][#6;X4]-[#6;X3]([!#1!#6])=[O;X1]");  //modified from ClasyFire's alpha-amino-acid-erivative-1 by adding X4 on the C2 carbon.
		queries.put("alpha-amino acid", 
				"[#7;A;X3,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X4]-[#6;X3]([#8;A;X2H1,X1-])=[O;X1]"); //modified from ClasyFire's alpha-amino-acid-1 by adding X4 on the C2 carbon.
		queries.put("tryptamine", 
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6;R1]-1=[#6;R1]-[#7;R1]-[#6;R2]-2=[#6;R2]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("flavonol", 
				"[H][#8]-[c;R1]1[#6;R1](=O)c2-,:cc-,:cc-,:c2[#8;R1][c;R1]1-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1]([#8;A;H1X2])[c;R1]-,:1");
		queries.put("flavone", 
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:1)-[#6;R1]=,:1[#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1](=O)[#6;R1]=,:1-[!#8]"); // PMID: 15914008
		queries.put("coumarin", 
				"[$(O=[#6]-1-[#8]-c2ccccc2-[#6]=[#6]-1),$(O=[#6]-1-[#8]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]=[#6]-1)]");
		queries.put("alkoxy group", 
				"[#6;A;X4][#8;X2][#6;A;H2X4][#7;X3]");
		queries.put("n_isopropyl group", 
				"[#6;A;H3X4][#6;H2X4]([#6;A;H3X4])-[#7;X3]");
		queries.put("thiono group", 
				"[*]=[SX1]");
		queries.put("s_2_hydroxycarboxylate", 
				"[#6]-[#6@H](-[#8;H1X2])-[#6;R0](-[#8;X1-])=[O;X1]"); 
		queries.put("hexoside", 
				"[#6]-[#8;H0X2]-[#6;R1]-1-[#8]-[#6;R1](-[#6])-[#6;R1](-[#8])-[#6;R1](-[#8])-[#6;R1]-1-[#8]");
		queries.put("phenothiazine", 
				"[$([#7]-1-c2ccccc2-[#16]-c2ccccc-12),$(N1C2=C(SC3=C1C=CC=C3)C=CC=C2)]"); // doi:10.1016/j.pharep.2015.04.005
		queries.put("phenoxazine", 
				"[$([#7]-1-c2ccccc2-[#8]-c2ccccc-12),$(N1C2=C(OC3=C1C=CC=C3)C=CC=C2)]");			
		queries.put("pyrrole", 
				"[$(c1ccnc1),$(N1C=CC=C1)]");			
		queries.put("pyrroline", 
				"C1CC=NC1");			
		queries.put("thioxanthene", 
				"[$([#6]-1-c2ccccc2-[#16]-c2ccccc-12),$(C1C2=CC=CC=C2SC2=CC=CC=C12)]");		
		queries.put("phthalazine", 
				"[#6]-1[#7]-[#7][#6]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-1-2");			
		queries.put("isoflavan", 
				"[H][C;R1]-,:1([H])-,:[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1][#6;R1]=,:[#6;R1]-,:1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");	

		// Added April 26th, 2016
		queries.put("pyridine", 
				"[$([#6;R1]-1=[#6;R1]-[#6;R1]=[#7;R1]-[#6;R1]=[#6;R1]-1),$([c;R1]1[c;R1][c;R1][n;R1][c;R1][c;R1]1)]");		
		
		queries.put("xanthone", 
				"O=[#6]1[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");		
		
		queries.put("fluorene", 
				"[#6]1[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");			
		
		queries.put("naphthoquinone", 
				"O=[#6]1[#6]=,:[#6][#6](=O)[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");
		
		queries.put("urea", 
				"[#7;X3;!$([#7][!#6])][#6;X3]([#7;X3;!$([#7][!#6])])=[O;X1]");
		
		queries.put("thiourea", 
				"[#7;X3;!$([#7][!#6])][#6;X3]([#7;X3;!$([#7][!#6])])=[S;X1]");
		
		queries.put("guanidine", 
				"[#7;A;v3X3,v4X4+][#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		
		queries.put("nitrosothiol", 
				"[#6]-[#16;X2]-[#7;X2]=[O;X1]");
		
		queries.put("sugar_pattern_1", 
				"[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
		
		queries.put("sugar_pattern_2", 
				"[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");		
		
		queries.put("sugar_pattern_combi", 
				"[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");		
		
		queries.put("sugar_pattern_2_reducing", 
				"[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
		
		queries.put("sugar_pattern_2_alpha", 
				"[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");		
		
		queries.put("sugar_pattern_2_beta", 
				"[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
		
		queries.put("thiolactam", 
				"[#6;R][#6;X3R]([#7;X3;$([H1][#6;!$(C=[O,N,S])]),$([H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[S;X1]");
		
		queries.put("carbothiolactone", 
				"[#6][#6X3R](=[SX1])[#16X2][#6;!$(C=[O,N,S])]");
		
		queries.put("imidothiolactone", 
				"[#6R][#6X3R](=,:[#7X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])[SX2][#6;!$(C=[O,N,S])]");
		
		queries.put("barbiturate", 
				"[$([#8]-,=[#6]-1-,=[#7]-[#6](-,=[#8])-[#6]-[#6](-,=[#8])-[#7]-1),$([#8]-,=[#6]-1-[#6]-[#6](-,=[#8])-[#7]=[#6](-[#8])-[#7]-1),$([#8]-,=[#6]-1-[#6]-[#6](-,=[#8])-[#7]-[#6](=O)-[#7]-1),$([$([O-]C1=NC(=O)NC(=O)C1),$([O-]C1=NC([O-])=NC(=O)C1),$([O-]C1=NC(=O)N=C([O-])C1)])]");		

		queries.put("benzodioxole", "[$(C1Oc2ccccc2O1),$(C1OOc2ccccc12)]");
		
		queries.put("glucosinolate","[#6][#6](-[#16]-[#6;R1]-1-[#8]-[#6](-[#6]-[#8])-[#6;R1]([!#1!#6;O,$([O-])])-[#6;R1]([!#1!#6;O,$([O-])])-[#6;R1]-1[!#1!#6;O,$([O-])])=[#7]/[#8]S([!#1!#6;O,$([O-])])(=O)=O"); // Slightly modified from the ClassyFire SMARTS
		
		queries.put("psoralen", "O=C1Oc2cc3occc3cc2C=C1"); // Guo LQ et al / Acta Pharmacol Sin 2004 Feb; 25 (2): 129-136
	
		queries.put("carbonate_salt", "[$([+1,+2,+3,+4,+5,+6,+7].[#8;A;X2H1,X1-]!@-[#6](!@-[#8;X1-])=O),$([#8;A;X2H1,X1-]!@-[#6](=O)!@-[#8;X2]-[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86])]");
		
		queries.put("carbonic_acid_derivative","[O,X]-[#6](-[O,X])=O");		
		
		queries.put("dibenzocyclooctadiene","[#6;R1]-1-[#6;R1]-[#6;R1]-c2cccc[c;R2]2-[c;R2]2ccccc2-[#6;R1]-1"); // 10.1124/dmd.104.000646 DMD December 2004 vol. 32 no. 12 1351-1358
		
		queries.put("benzenesulfonamide","[#7]S(=O)(=O)!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1"); // 10.1124/dmd.104.000646 DMD December 2004 vol. 32 no. 12 1351-1358
		
		queries.put("3_phenylpyrazole","c1cc(nn1)-!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		
		queries.put("hetero_n_basic_h","[nX3H1+0]"); // doi: 10.3389/fphar.2015.00123 Modeling of interactions between xenobiotics and cytochrome P450 (CYP) enzymes
		
		queries.put("hetero_n_basic_no_h","[nX3H0+0]");	// doi: 10.3389/fphar.2015.00123 Modeling of interactions between xenobiotics and cytochrome P450 (CYP) enzymes

		queries.put("hetero_non_basic","[nX2,nX3+]");
		
		queries.put("piperidine","[#6;R1]-1-[#6;R1]-[#6;R1]-[#7;R1]-[#6;R1]-[#6;R1]-1");
		queries.put("4-aminopiperidine_monocyclic","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6;R1]-1-[#6;R1]-[#6;R1]-[#7;R1]-[#6;R1]-[#6;R1]-1");
		queries.put("4-aminopiperidine","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6]-1-[#6]-[#6]-[#7]-[#6]-[#6]-1");
		queries.put("halobenzene","[F,Cl,Br,I]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("benzene","[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");

		// added on May 8, 2016	
		
		queries.put("14_dihydropyridine","C1C=CNC=C1");
		queries.put("dihydropyrimidinethione","[O;X1R0]=[#6;R1]-1-[#7]-[#6]-[#6]=[#6]-[#7]-1");
		queries.put("quinazoline","[#6]=,:1[#6]=,:[#6][#6]=,:2[#7][#6]=,:[#7][#6][#6]=,:2[#6]=,:1");
		queries.put("2_phenylpyrimidine_unfused_benzene","[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)!@-[#6]-1=[#7]-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("3_phenylpyrimidine_unfused_benzene","[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6]-1=[#7]-[#6]=[#7]-[#6]=[#6]-1");

		queries.put("13_thiazolidine","C1CSCN1");
		queries.put("carboxylic-acid-imide-N-substituted","[#6]-[#7;X3](-[#6;X3](-[#6,#1])=[O;X1])-[#6;X3](-[#6,#1])=[O;X1]");
		queries.put("carboxylic-acid-imide-N-unsubstituted","[H][#7;X3](-[#6;X3](-[#6,#1])=[O;X1])-[#6;X3](-[#6,#1])=[O;X1]");
		queries.put("vinylogous-halide","[#9,#17,#35,#53;A;X1]-,:[#6;X3]=,:[#6;X3][#6;X3]=[O;X1]"); // http://ccj.springeropen.com/articles/10.1186/1752-153X-7-49
		queries.put("vinylogous-amide","[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;X3]=,:[#6;X3][#6;X3]=[O;X1]");
		queries.put("tetralin","[#6]1=,:[#6][#6]=,:[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]2=,:[#6]1");
		queries.put("n_phenylurea","[#7]!@-[#6](=O)!@-[#7]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("piperazine","C1CNCCN1");
		queries.put("azacycle","[#6]@[#7]");
		queries.put("oxacycle","[#6]@[#8]");
		queries.put("organoheterocycle","[!C!c!R0]");
		queries.put("filter_bridged_rings_66_heterocyclic","*~1~*~*~*(~*~*~1)!@*~1~*~*~*~*~*~1");
		queries.put("filter_bridged_rings_65_heterocyclic","*~1~*~*~*(~*~1)!@*~1~*~*~*~*~*~1");
		queries.put("filter_bridged_rings_55_heterocyclic","*~1~*~*~*(~*~1)!@*~1~*~*~*~*~1");

		queries.put("filter_fused_rings_67_heterocyclic","*~1~*~*~*~2~*~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_66_heterocyclic","*~1~*~*~*~2~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_65_heterocyclic","*~1~*~*~2~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_55_heterocyclic","*~1~*~*~2~*~*~*~*~2~*~1");
		
		
		
		queries.put("aliphatic_chain_13","[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]");
		queries.put("aliphatic_chain_4","[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]");
		queries.put("ortho_substituted_ring","*-!:aa-!:*");
		queries.put("meta_substituted_ring","*-!:aaa-!:*");
		queries.put("para_substituted_ring","*-!:aaaa-!:*");
		queries.put("aryl_fluoride","F[$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("1-methyl-23-dihydro-1H-imidazole-2-thione_monocyclic","[#6;A;H3X4][#7;A;R1]1[#6;A;R1]=[#6;A;R1][#7;A;R1][#6;A;R1]1=[S;X1]");
		queries.put("23-dihydro-1H-imidazole-2-thione_monocyclic","[#6;A;R1]1[#7;A;R1][#6;A;R1]=[#6;A;R1][#7;A;R1]1");
		queries.put("terminal_alkyne","[H]C!@#C*");
		queries.put("w_1_acetylene","[#6;A;H3X4]C!@#C*");			
		queries.put("dithiocarboxylic_acid_amide","[#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])])][#16;X2]-[#6;X3]([#6,#1;A])=[S;v2]");
		queries.put("thiocarboxylic_acid_amide","[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put("terminal_cc_bond","[#6]-[#6]=[#6;A;H2X3]");	
		queries.put("thiophene_sulfoxide","O=S1[#6]=[#6]-[#6]=[#6]1");
		queries.put("trifluoroalkyl_group","[#6;A;X4][C;X4](F)(F)F");
		queries.put("12_dithiole_3_thione","S=C1SSC=C1");	
		queries.put("piperidine_monocyclic","[#6;A;R1]1[#6;A;R1][#6;A;R1][#7;A;R1][#6;A;R1][#6;A;R1]1");
		queries.put("sulfenyl_group","[#6]-[#16;v2X2R0]-*");
		queries.put("sulfenyl_halide","[#6]-[#16;v2X2R0]-[F,Cl,Br,I]");	
		queries.put("thiophosphate_diester","[#6]-[#8;X2][P;X4]([#8;A;H1X2])(=[S;X1])[#8;X2]-[#6]");
		queries.put("45_dihydro_123_oxadiazol_5_one","O=C1CN=NO1");
		queries.put("alkyl_sulfanyl_carbothioyl_amine","[#7;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;X3R0](-[#16;X2])=[S;X1R0]");	
		queries.put("123_triazolidine","[#6]-1-[#6]-[#7]=[#7]-[#7]-1");
		queries.put("1_amino_123_triazolidine","NN1CCN=N1");
		queries.put("n_cyclobutylbenzylamine","[#6;R0](-[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("n-cyclopropylbenzylamine","[#6;R0](-[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("cyclopropylamine","[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-1");
		queries.put("cyclobutylamine","[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-1");	
		queries.put("123_thiadiazole","c1csnn1");	
		queries.put("124_thiadiazole","c1ncsn1");	
		queries.put("125_thiadiazole","c1cnsn1");	
		queries.put("134_thiadiazole","c1nncs1");	
		queries.put("xanthate","[$([!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86;+].[#6;A;X4][#8;X2]-[#6;X3](!@-[#16-])=[SX1]),$([#6;A;X4][#8;X2]-[#6;X3](=[SX1])!@-[#16;X2]-[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86])]");	
		queries.put("pyrazine","c1cnccn1");	 // pp440,445 in Robert M. Rydzewski. Real World Drug Discovery: A Chemist's Guide to Biotech and Pharmaceutical Research
		queries.put("11_disubstituted_hydrazine","[#6]-[#7;X3](-[#6])!@-[#7;A;H2X3]");	// Inhibition of Cytochrome P450 Enzymes p.266
		queries.put("monosubstituted_hydrazine","[#6]-[#7;H1X3]!@-[#7;A;H2X3]");	// Inhibition of Cytochrome P450 Enzymes p.266
		queries.put("carboxylic_acid_hydrazide","[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)]-[#7](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#7](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#6](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])=O");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("12-dihydroquinoline","[#6;A;X4]1[#7]-[#6]-2=[#6](-[#6]=[#6]1)-[#6]=[#6]-[#6]=[#6]-2");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano 
		queries.put("alkylhydrazine","[#6;A;X4][#7;X3](-[#6,#1])-[#7;X3](-[#6,#1])-[#6,#1]");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("hydrazone","[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#7;X2]=[#6]");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("lactone","[#6;!$(C=[O,N,S])][#8;X2][#6;X3R]([#6])=[O;X1]");	
		queries.put("lactam","[#6;R][#6;X3R]([#7;X3;$([H1][#6;!$(C=[O,N,S])]),$([H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[O;X1]");	
		queries.put("linear_furanocoumarin_pattern_1","O=[#6]-1-[#8]-[#6]-2=[#6](-[#6]=[#6]-1)-[#6;R1]=[#6]-1-[#6]=[#6;R1]-[#8]-[#6]-1=[#6;R1]-2");	// J Pharm Pharm Sci. 2001 Sep-Dec;4(3):217-27.Inhibition of human CYP3A4 activity by grapefruit flavonoids, furanocoumarins and related compounds.
		queries.put("stilbene","[$(*=[#6]-1-[#6]=[#6]-[#6]=[#6]-[#6;R1]-1!@-[#6]!@=[#6]~[c;R1]1ccccc1),$([#6](!@=[#6]~/[c;R1]1ccccc1)!@-[#6;R1]-1=[#6]-[#6]=[#6]-[#6]=[#6]1)]");	// Mol Nutr Food Res. 2008 Jun;52 Suppl 1:S77-83. doi: 10.1002/mnfr.200700202.
		queries.put("anthocyanidin","[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1+]=,:1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("o_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6]-1-[#6;R0]=[#6;R0]-[#6;R0](-*)=O");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("m_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6;R1]-1=[#6]-[#6](-[#6;R0]=[#6;R0]-[#6;R0](-*)=O)=[#6;R1]-[#6;R1]=[#6;R1]-1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("p_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6](-[#6;R0]=[#6;R0]-[#6;R0](-*)=O)-[#6]=[#6;R1]-1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("o_hydroxybenzoic_acid_der","[#8;A;X2R0][#6]-1=[#6](-[#6]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("m_hydroxybenzoic_acid_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6]-[#6](=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("p_hydroxybenzoic_acid_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6](-[#6]=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("1_benzopyran","[#6]-1-,=[#6]-c2ccccc2-[#8]-[#6]-1");	
		queries.put("2_benzopyran","[#6]-1-,=[#6]-c2ccccc2-[#6]-[#8]-1");	
		queries.put("chalcone_or_der","[#8;R0]-,=[#6;R0](-,=[#6;R0]-,=[#6;R0]!@-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1)!@-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1");	
		queries.put("phenoxy_compound","*-[#8;X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("thienopyridine","[$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#6]=[#7]-[#6]=[#6]-2),$([#6]~1~[#7]~[#6]-[#6]-2=[#6](-[#6]~1)-[#16]-[#6]=[#6]-2),$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#6]=[#6]-[#7]=[#6]-2),$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#7]=[#6]-[#6]=[#6]-2)]");  //http://www.sciencedirect.com/science/article/pii/S0968089614007147
		queries.put("1_benzothiophene","[#6]1=,:[#6][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#16]1"); //http://www.sciencedirect.com/science/article/pii/S0968089614007147
		queries.put("2_benzothiophene","[#6]=,:1[#16][#6]=,:[#6]2[#6]=,:[#6][#6]=,:[#6][#6]=,:12"); //http://www.sciencedirect.com/science/article/pii/S0968089614007147
//		queries.put("kavalactone","");	 // Robert S. FOTI and Jan L.WAHLSTROM; The role of dietary supplements in cytochrome P450-mediated drug interactions 
//		queries.put("ginkgolide_or_bilabolide","");		 // Robert S. FOTI and Jan L.WAHLSTROM; The role of dietary supplements in cytochrome P450-mediated drug interactions 
		queries.put("cyclic_carboxylic_ester","[#6;R;!$(C=[O,N,S])]-[#8;X2R][#6;A;X3R;$([R0][#6]),$([H1R0])]=[O;X1R0]");	// for macrolides
//		queries.put("macrocycle","[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11]");	// for macrolides	
		queries.put("macrocycle_r12_r40","[r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r40]");	// for macrolides
		queries.put("chloroethylene_deriv","[Cl]-[#6](-[#6,#1])=[#6](/Cl)[#17,#1;A]");
		queries.put("conjug_bond","[*,#1]-[#6]=[#6]-[#6]=[#6]-[*,#1]");
		
		queries.put("curcumin","[$([#8]~[#6;A](~[#6;A;R0]~[#6;A]~[#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A]~[#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A]~[#6;A;R0][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)]");	
		queries.put("unfused_benzene","[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");	
		queries.put("s_in_aromatic_5_ring_with_lone_pair","[sX2r5]");	
		queries.put("o_in_aromatic_5_ring_with_lone_pair","[oX2r5]");			
		queries.put("n_in_aromatic_5_ring","[nX2r5]");	
		//	 1 However, no N-oxides of nitrogen atoms in five-membered aromatic rings have been found. Instead, such nitrogen atoms tend to bind directly to the heme iron atom, leading to CYP inhibition. http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=3412406&fileOId=3412409
	
		queries.put("s_in_aromatic_6_ring_with_lone_pair","[sX2r6]");	
		queries.put("o_in_aromatic_6_ring_with_lone_pair","[oX2r6]");			
		queries.put("n_in_aromatic_6_ring","[nX2r6]");
		queries.put("pattern_c1", 
				"[#6][#8,#16;A][#6;H2X4R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // Tioconazole Miconazole sulconazole Miconazole
		queries.put("pattern_c2","[#6;A]-,=1-,=[#6;A][#6;A;R2]2[#6;A]-,=[#6;A]-,=[#6;A]-,=[#6;A][#6;A;R2]2[#6;A]-,=1");
		queries.put("pattern_c3","*-[#6;X3](=O)!@-[#6;X3](-[F,Cl,Br,I])-[F,Cl,Br,I]");	
		queries.put("pattern_c4","[#6]-[#16;X2]!@-[#6;X3](-[#6])=O");	// from BioTransformer fingerprint-based clustering
		queries.put("pattern_c5","[#6]-[#7;R0]=[#6;R0]-[#7;R0]-[#8;H1X2R0]");	
		queries.put("cyclic","[*;R]");		
		queries.put("134_oxadiazole","c1nnco1");	
		queries.put("124_oxadiazole","c1ncno1");		
		queries.put("cc_double_bond","[#6]=[#6]");	
		queries.put("ring_of_size3","*~1~*~*~1");
		queries.put("ring_of_size5","*~1~*~*~*~*~1");	
		queries.put("ring_of_size6","*~1~*~*~*~*~*~1");
		queries.put("ring_of_size7","*~1~*~*~*~*~*~*~1");
		queries.put("heteroring_of_size5","*~1~*~*~[!#1!#6]~*~1");	
		queries.put("heteroring_of_size6","*~1~*~*~[!#1!#6]~*~*~1");	
		queries.put("heteroring_of_size7","*~1~*~*~*~[!#1!#6]~*~*~1");	
		queries.put("oxolane","C1CCOC1");
		queries.put("imidazoline","C1CN=CN1");
		queries.put("benzo_13_dioxolane","[#6]-1-[#8]-c2ccccc2-[#8]-1");

		
		// J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280
		queries.put("pattern_c6","*-*(-*)-[*;X2]-*(-*)-*");		//atom has two neighbors, each with three or more neighbors (including the central atom). (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("pattern_c7","[!#1!#6]*[!#1!#6]"); // atom with at least two heteroatom neighbors (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("pattern_c8","*!@-*!@-*"); //* atom with more than one chain bond (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("n_halide_bond","[#7]-[F,Cl,Br,I]"); // N-X bond (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("nn_double_bond","[#7]=[#7]");
		queries.put("cn_double_bond","[#6]=[#7]");
		queries.put("cs_double_bond","[#6]=[#16]");
		queries.put("cn_single_bond","[#6]-[#7]");
		queries.put("co_single_bond","[#6]-[#8]"); // C-O single bond
		
		
		// Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition (PMID: 26108525). J Chem Inf Model. 2015 Jul 27;55(7):1426-34. doi: 10.1021/acs.jcim.5b00130. Epub 2015 Jul 2.
		queries.put("r5_1ht_sat","[#6]-1-[#6]-[#6][!#1!#6][#6]-1"); // saturated 5-membered ring with 1 heteroatom. Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition 
		queries.put("r5_2ht_sat","[$([#6]-1-[#6][!#1!#6][!#1!#6][#6]-1),$([#6]1-[#6][!#1!#6][#6][!#1!#6]1)]"); // saturated 5-membered ring with 2 heteroatoms
		queries.put("r5_3ht_sat","[$([#6]1-[#6][!#1!#6][!#1!#6][!#1!#6]1),$([#6]1[!#1!#6][#6][!#1!#6][!#1!#6]1),$([#6]1[!#1!#6][#6][!#1!#6][!#1!#6]1)]");  // saturated 5-membered ring with 3 heteroatoms
		queries.put("r5_4ht_sat","[#6]1[!#1!#6][!#1!#6][!#1!#6][!#1!#6]1");  // saturated 5-membered ring with 4 heteroatoms
		queries.put("r5_1ht_aro","c1cc[!#1!#6]c1"); // aromatic 5-membered ring with 1 heteroatom. Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition 
		queries.put("r5_2ht_aro","[$(c1c[!#1!#6][!#1!#6]c1),$(c1c[!#1!#6]c[!#1!#6]1)]"); // aromatic 5-membered ring with 2 heteroatoms
		queries.put("r5_3ht_aro","[$(c1c[!#1!#6][!#1!#6][!#1!#6]1),$(c1[!#1!#6]c[!#1!#6][!#1!#6]1),$(c1[!#1!#6]c[!#1!#6][!#1!#6]1)]");  // aromatic 5-membered ring with 3 heteroatoms
		queries.put("r5_4ht_aro","c1[!#1!#6][!#1!#6][!#1!#6][!#1!#6]1");  // aromatic 5-membered ring with 4 heteroatoms
		queries.put("carbon_atom","[#6]");
		queries.put("oxygen_atom","[#8]");
		queries.put("nitrogen_atom","[#7]");
		queries.put("sulfur_atom","[#16]");
		queries.put("phos_atom","[#15]");		

		queries.put("het_atom","[!#1!#6]");
		queries.put("organometallic","[#6]~[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86]");		

		queries.put("aliph_carboxylic","[#6;X4]-[#6]([#8;A;X2H1,X1-])=O");
		queries.put("sat_carb_r6","C1CCCCC1"); //  saturated carbon-only ring size 6

		queries.put("1_sulfenylpiperidine","O=[S][#7]-1-[#6]-[#6]-[#6]-[#6]-[#6]-1");
		queries.put("pattern_c10","CCC(C)N"); // *
		queries.put("n_o_single_bond","[#7]-[#8]");
		queries.put("ketone","[#6]-[#6;X3](-[#6])=[O;X1]");

		// Generate using SOCO Chepelevet al.
		
		queries.put("pattern_c9","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;X2;!$([#7][!#6])])[#8;X2;!$([#8][!#6]),OX1-]");
//		queries.put("pattern_c9","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;X2;!$([#7][!#6])])[#8;X2;!$([#8][!#6]),OX1-]");
		queries.put("pattern_c53","[#6;R][#6;R][#6;R][#6;R][#6;R]");
		queries.put("pattern_c11","[#6][#6][#6]([#6])[#6]");
		queries.put("pattern_c12","[#6]~[#6](~[#6])[#6]");
		queries.put("pattern_c13","[#6]-[#9,#17,#35,#53,#85]");
		queries.put("pattern_c14","[#6X3](=[OX1])[#6X3]=,:[#6X3][#7,#8,#16,F,Cl,Br,I]");
		queries.put("pattern_c15","[$([*@@](~*)(~*)(*)*),$([*@@H](*)(*)*),$([*@@](~*)(*)*),$([*@@H](~*)~*)]");
		queries.put("pattern_c16","C=CCCCC");
		queries.put("pattern_c17","C=CC=CCC");
		queries.put("pattern_c18","[#6X4H3]CN");
		queries.put("pattern_c19","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("pattern_c20","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("pattern_c21","[!#6][#6X3](=[!#6])[!#6]");
		queries.put("pattern_c22","[#7X3;!$([#7][!#6])][#6X3](=[OX1])[#7X3;!$([#7][!#6])]");
		queries.put("pattern_c23","[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put("pattern_c24","[R]~[*;!R]~[*;!R]~[*;!R]~[R]");
		queries.put("pattern_c25","[R]~[*;!R]~[*;!R]~[*;!R]~[*;!R]~[R]");
		queries.put("pattern_c26","[#6][#7][#6X4H3]");
		queries.put("pattern_c27","[#6X4H3][#6][#7]");
		queries.put("pattern_c28","[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("pattern_c29","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]");
		queries.put("pattern_c30","[#6X3](=[OX1])[#6X3]=,:[#6X3][#6;!$(C=[O,N,S])]");
		queries.put("pattern_c31","[$([#7X2,OX1,SX1]=,:**=,:*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=**=*),$([#7X3,OX2,SX2;!H0]*=,:**:n)]");
		queries.put("pattern_c32","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*)]");
		queries.put("pattern_c33","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*),$([*@@](~*)(~*)(*)*),$([*@@H](*)(*)*),$([*@@](~*)(*)*),$([*@@H](~*)~*)]");
		queries.put("pattern_c34","OCC");
		queries.put("pattern_c35","N(C)C");
		queries.put("pattern_c36","[*](=O)[!#6][*;!H]");
		queries.put("pattern_c37","[$([#6X3H0][#6]),$([#6X3H])](=[!#6])[!#6]");
		queries.put("pattern_c38","[$([#7X2,OX1,SX1]=*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=*),$([#7X3,OX2,SX2;!H0]*:n)]");
		queries.put("pattern_c39","[#6]~[#6]~[#7]");
		queries.put("pattern_c40","[#6]#[#7]");			
		queries.put("pattern_c41","[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("pattern_c42","[FX1][CX4;!$([H0][Cl,Br,I]);!$([F][C]([F])([F])[F])]([FX1])([FX1])");	
		queries.put("pattern_c42","[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("pattern_c43","[FX1][CX4;!$([H0][Cl,Br,I]);!$([F][C]([F])([F])[F])]([FX1])([FX1])");
		queries.put("pattern_c44","OCCCC");
		queries.put("pattern_c45","[#6][#6][#8][#6]");
		queries.put("pattern_c46","[#6]~[#8]~[#6]");
		queries.put("pattern_c47","[#6]~[#6]~[#6]~[#8]~[#6]");
		queries.put("pattern_c48","[#6][#6][#8][#6][#6]");
		queries.put("pattern_c49","[#6;R][#6;R][#6;R][#8;R]");
		queries.put("pattern_c50","[#6;R][#6;R][#6;R][#6;R][#8;R]");
		
// these should also be added to the fingerprint for SOM prediction
//		queries.put("acetanilide","[#6;X4]-[#6;X3](=[O;X1])-[#7;X3]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("aniline_der","[#7;X3;H1,H2]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("o_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1-[#8;H1X2]");   // Kalgutkar AS. et al (PMID: 15975040; Current Drug Metabolism, 2005, 6, 161-225)
		queries.put("p_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1](-[#8;H1X2])[c;R1][c;R1]1"); // Kalgutkar AS. et al (PMID: 15975040; Current Drug Metabolism, 2005, 6, 161-225)
		queries.put("aryl_urea","[#6;a]-[#7;X3]-[#6;X3](-[#7;X3;!$([#7][!#6])])=[O;X1]");
		queries.put("aryl_carbamate","[#6;a;R1]!@-[#7]!@-[#6](=O)-[#8]-[#6;!$(C=[O,N,S])]");
		queries.put("hydroquinone","[H][#8]!@-[c;R1]1[c;R1](-[*,#1;!$([OH])])[c;R1](-[*,#1;!$([OH])])[c;R1](!@-[#8][H])[c;R1](-[*,#1;!$([OH])])[c;R1]1-[*,#1;!$([OH])]");
		queries.put("n_nitrosourea","[#7]!@-[#6;X3](!@=[O;X1])!@-[#7;X3]!@-[#7;X2]!@=O");
		queries.put("n_mustard","[$(Cl[#6]-[#6]-[#7]-[#6]-[#6]Cl),$(F[#6]-[#6]-[#7]-[#6]-[#6]F),$(Br[#6]-[#6]-[#7]-[#6]-[#6]Br),$(I[#6]-[#6]-[#7]-[#6]-[#6]I)]");
		
// Contributions of Human Enzymes in Carcinogen Metabolism : Chem Res Toxicol. 2012 July 16; 25(7): 1316–1383. doi:10.1021/tx300132k.
		queries.put("aza_aromatic","[#7;a;R]");
		queries.put("oxazaphosphorin_2_amine_2_oxide","[#7][P;R1]1(=[O;X1])[#7;X3R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#8;R1]1");
		queries.put("nitrosamine","[#7;!$(N*=O)]-[#7;X2]=[O;X1]");
		queries.put("beta_lactam","O=[#6]-1-[#6]-[#6]-[#7]-1");
//		queries.put("michael_acceptor","");	 

		
		queries.put("pattern_c51","[#6;H3X4]-[#6]-[#7;X3](-[#6;H3X4])-[#7;X2]=[O;X1]");
		queries.put("pattern_c52","[#6;A;H3X4][#6]-[#7;X3]([#6;A;H3X4])-[#6][#6;A;H3X4]");	
		

		// Electron donating group / Electron withdrawing group :  http://www.chem.ucalgary.ca/courses/350/Carey5th/Ch12/ch12-8b.html
		// Add those
		// Read this link: http://www.chem.ucalgary.ca/courses/350/Carey5th/Ch12/ch12-8d.html

//		Strongly activating EDG
		queries.put("phenoxide_ortho_subs","[#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*");
		queries.put("phenoxide_para_subs","[c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1");
		queries.put("prim_amin_ortho_subs","[#7;A;H2X3][c;R1]:[*;R1]!@-*");
		queries.put("prim_amin_para_subs","[#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*");
		queries.put("ter_amin_ortho_subs","*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]");
		queries.put("ter_amin_para_subs","*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]");	
		queries.put("phenol_ortho_subs","[#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*");
		queries.put("phenol_para_subs","[#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1");
		queries.put("phenolether_ortho_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenolether_para_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1");		
		queries.put("phenolether_ortho_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenolether_para_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1");		

//		 Moderately activating EDG		
		queries.put("anilide_n_ortho_subs","[#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("anilide_n_para_subs","[#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1");		
		queries.put("phenolester_ortho_subs","[#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenolester_para_subs","[#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1");	

//		 Weakly activating EDG	
		queries.put("phenylalkyl_ortho_subs","[#7;X3]!@-[#6]([#6,#1;A])!@=O.[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenylalkyl_para_subs","[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1");
		queries.put("biphenyl_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("biphenyl_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("phenylvinyl_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]");
		queries.put("phenylvinyl_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1");
		
		
//		Weakly deactivating EWG
		queries.put("phenylhalide_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]");		
		queries.put("phenylhalide_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1");
		
		// Modelately deactivating EWG	
		queries.put("benzaldehyde_meta_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1");
		queries.put("phenylketone_meta_subs","[#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");		
		queries.put("benzoate_meta_subs","[#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzoic_acid_meta_subs","[#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzoyl_chloride_meta_subs","[Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");		

		// Strongly deactivating EWG		
		
		queries.put("benzotrifluoride_meta_subs","[F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzotrichloride_meta_subs","[Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");			
		queries.put("benzonitrile_meta_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]");	
		queries.put("benzenesulfate_meta_subs","[#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("anilinium_meta_subs","[#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("phenyl_1qam_meta_subs","[#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	// benzene C1-substituted with a quaternary ammoniuam salt and C3- substituted any atom.
		queries.put("nitrobenzene_meta_subs","[#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	

		
		
		queries.put("edg","[$([#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1),$([#7;A;H2X3][c;R1]:[*;R1]!@-*),$([#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*),$(*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$(*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#7;X3]!@-[#6]([#6,#1;A])!@=O.[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1)]");	
		queries.put("ewg","[$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1),$([#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]),$([#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1)]");	
		
		
		
		// CYP isoform specificity toward drug metabolism: analysis using common feature hypothesis
		// J Mol Model (2012) 18:709–720: DOI 10.1007/s00894-011-1105-5
		
		queries.put("isopropyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-*");
		queries.put("isobutyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-[#6;A;X4;H2,H3]");
//		queries.put("edg","");
//		queries.put("ewg","");		
		// http://www.ifm.liu.se/compchem/msi/doc/life/catalyst46/tutorials/11_excludeOrTool.doc.html
		// prim. amine OR sec. amine OR ter. amine OR amidine OR guadidino OR amidineH OR (positively charged center except dipole)
		queries.put("posionazable","[$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]),$([#7;AH1;v3X3,v4X4+][#6;X3]([#7;AH1;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]),$([#7;A;H1X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6]);!$([C][N])]=[#7;A;X2;!$(NC=[O,S])]);!$([+1,+2,+3,+4,+5,+6,+7]~[+0,+1,+2,+3,+4,+5,+6,+7])]");		
//		queries.put("","");
//		queries.put("","");		
		
		
		queries.put("acyl_carnitine","[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6](!@-[#6;A;H2X4]!@-[#6](!@-[#8;A;X1-,X2H1])!@=O)!@-[#8]!@-[#6](!@-[#6,#1;A])!@=O");
		queries.put("acyl_choline","[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#8]!@-[#6]([#6,#1;A])!@=O");
		queries.put("thiohemiacetal","[SX2]([#6;!$(C=[O,S,N])])[CX4;!$(C(S)(S)[!#6])][OX2H]");		
		
		queries.put("thioacetal","[#6]-[#16;X2][C;X4]([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#16;X2]-[#6]");
		queries.put("hemiaminal","[NX3v3;!$(NC=[#7,#8,#15,#16])]([#6])[CX4;!$(C(N)(N)[!#6])][OX2H]");
		queries.put("aminal","[*,#1;CX4,#1]-[#7](-[*,#1;CX4,#1])C([*,#1;CX4,#1])([*,#1;CX4,#1])[#7](-[*,#1;CX4,#1])-[*,#1;CX4,#1]");
		queries.put("hemiacetal","[H][#8]C([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[#6]");		
		queries.put("acetal","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]C([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("imidothioester","[#6;!$(C=[O,N,S])]-[#16;X2][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");	
		queries.put("imidothioacid","[$([SX2H]),$([SX1-])][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");
		queries.put("Imidothioacid_cyclic","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])[$([SX2H]),$([SX1-])]");		
		queries.put("imidothiolactone","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])-[#16;X2]-[#6;!$(C=[O,N,S])]");
		queries.put("imidolactam","[#6][#6;X3R;$([H0](=[NX2;!$(N(=[#6X3][#7X3])C=[O,S])])[#7X3;!$(N([#6X3]=[#7X2])C=[O,S])]),$([H0](-[NX3;!$(N([#6X3]=[#7X2])C=[O,S])])=,:[#7X2;!$(N(=[#6X3][#7X3])C=[O,S])])]");	
		queries.put("imidoyl_halide","[#9,#17,#35,#53;A;X1][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");
		queries.put("imidoyl_halide_cyclic","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])-,:[#9,#17,#35,#53;A;X1]");
		queries.put("carbxoxylic_acid_amidrazone","[$([$([#6X3][#6]),$([#6X3H])](=[#7X2v3])[#7X3v3][#7X3v3]),$([$([#6X3][#6]),$([#6X3H])]([#7X3v3])=[#7X2v3][#7X3v3])]");	
		queries.put("alpha_hydroxyacid","[H][#8;v2]-[#6](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#6](=O)-[#8][H]");
		queries.put("alpha_hydroxyaldehyde","[#6]-,=[#6;R0](-[#8;R0][H])-[#6;R0]([H])=[O;R0]");
		queries.put("carboxylic_acid_orthoester","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]C([#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("carbonic_acid_ester_halide","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](=[O;X1])-[#8;X2][#9,#17,#35,#53;A;X1]");
		queries.put("carbonic_acid_monoester","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](-[$([OX2H]),$([OX1-])])=[O;X1]");	
		queries.put("carbonic_acid_diester","[#6;!$(C=[O,N,S])][#8;X2][#6;X3](=[O;X1])[#8;X2][#6;!$(C=[O,N,S])]");
		queries.put("thiocarbonic_acid_derivative","[$([*,#1;!$(C=[O,N,S])][#8;X2][#6;X3](=O)[#16;X2][*,#1;!$(C=[O,N,S])]),$([*,#1;!$(C=[O,N,S])][#8;X2][#6;X3](=[O;X1])[#16;X2][*,#1;!$(C=[O,N,S])])]");	
		queries.put("thiocarbonic_acid_ester_halide","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](=[S;X1])-[#8;X2][#9,#17,#35,#53;A;X1]");
		queries.put("thiocarbonic_acid_monoester","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](-[$([OX2H]),$([OX1-])])=[S;X1]");	
		queries.put("thiocarbonic_acid_diester","[#6;!$(C=[O,N,S])][#8;X2][#6;X3](=[S;X1])[#8;X2][#6;!$(C=[O,N,S])]");
		queries.put("isourea","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;A;X2;!$([#7][!#6])])[#8X2!$([#8][!#6]),OX1-]");	
		queries.put("aziridinone","[O;X1]=[#6;X3]-1-[#6]-[#7]-1");
			
		queries.put("heteroaromatic_s","[sX2]");
		queries.put("heteroaromatic_o","[o]");	

		queries.put("nitrosamide","O=[#6]-[#7]-[#7;X2]=[O;X1]");	
		queries.put("sulfuric_acid_monoester","[#6]-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		queries.put("sulfuric_acid_diester","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]S(=O)(=O)[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("sulfuric_acid_diamide","[#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][S;X4]([#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])(=[O;X1])=[O;X1]");
		queries.put("sulfuric_acid_amide_ester","[#15;A;X4;$([H3]=[OX1]),$([H2](=[OX1])[#6]),$([H1](=[OX1])([#6])[#6]),$([H0](=[OX1])([#6])([#6])[#6])]");
		queries.put("sulfuric_acid_derivative","[#7,#8;A]S([#7,#8;A])(=O)=O");	
		queries.put("sulfonic_acid","[C,$([cX3](:*):*),$([cX2+](:*):*)][S;v6]([!#1!#6;!S])(=O)=O");
		queries.put("hydroxylamine_o_sulfonic_acid","[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		queries.put("sulfonic_acid_ester","[#6;!$(C=[O,N,S])]-[#8][S;v6]([#6,#1;A])(=[O;X1])=[O;X1]");	
		queries.put("sulfonic_acid_derivative","[C,$([cX3](:*):*),$([cX2+](:*):*)][S;v6]([!#1!#6;!S])(=O)=O");
		queries.put("phosphine","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#15;v3X3](-[C,$([cX3](:*):*),$([cX2+](:*):*)])-[C,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("phosphine_imide","[$([#6]P([#6,#1;A])([#6,#1;A])=[#7;X2][#6,#1;A]),$([#6][P;X4+]([#6,#1;A])([#6,#1;A])[#7-][#6,#1;A])]");	
		queries.put("phosphine_sulfide","[$([#6]P([#6,#1;A])([#6,#1;A])=[S;X1]),$([#6][P;X4+]([#16;X1-])([#6,#1;A])[#6,#1;A])]");
		queries.put("phosphinoxide","[#15;A;X4;$([H3]=[OX1]),$([H2](=[OX1])[#6]),$([H1](=[OX1])([#6])[#6]),$([H0](=[OX1])([#6])([#6])[#6])]");	
		queries.put("phosphonium","[P+;!$([P]~[!#6]);!$([P]*~[#7,#8,#15,#16])]");
		queries.put("phosphonic_acid_derivative","[O,N,X]P([O,N,X])([C,$([cX3](:*):*),$([cX2+](:*):*)])=O");	
		queries.put("phosphonic_acid_ester","[#6;!$(C=[O,N,S])]-[#8;X2]P([#6;!$(C=[O,N,S])])([O,N,X])=[O;X1]");
		queries.put("phosphonic_acid_diester","[#6;!$(C=[O,N,S])]-[#8;X2][#15;A;X4;$([H1]),$([H0][#6])](=[O;X1])[#8;X2]-[#6;!$(C=[O,N,S])]");	
		queries.put("phosphonic_amide_ester","[#6;!$(C=[O,N,S])]-[#8;X2][#15;A;X4;$([H1]),$([H0][#6])]([#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[O;X1]");
		queries.put("phosphonic_dihydrazide","[#6]P(=O)([#7]([#6,#1;A])-[#7]([#6,#1;A])[#6,#1;A])[#7]([#6,#1;A])-[#7]([#6,#1;A])[#6,#1;A]");	
		queries.put("trialkylborane","[#6;X4]-[#5](-[#6;X4])-[#6;X4]");
//	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
//
		//  From OpenBabel's FP4
		
		//hits chloromethylenethers and other reactive alkylating agents
		queries.put("Halogen_acetal_like","[NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1]");
		
		// # includes all of the above and other combinations (S-C-N, hydrates, ...), but still no aminomethylenesters and similar
		queries.put("acetal_like","[NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1,NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		// # also reactive alkylating agents. Acid does not have to be carboxylic acid, also S- and P-based acids allowed
		queries.put("halogenmethylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1]");
		
		// # Same as above, but N,O or S instead of halogen. Ester/amide allowed only on one side 

		queries.put("nos_methylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		// # Combination of the last two patterns
		queries.put("hetero_methylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1,NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		queries.put("cyanhydrine","[NX1]#[CX2][CX4;$([CH2]),$([CH]([CX2])[#6]),$(C([CX2])([#6])[#6])][OX2H]");
		
		// # includes aminals, silylacetals, ketenesters, etc. C=C DB is not aromatic, everything else may be
		queries.put("Ketenacetal","[#7X2,#8X3,#16X2;$(*[#6,#14])][#6X3]([#7X2,#8X3,#16X2;$(*[#6,#14])])=[#6X3]");
				
		queries.put("spiro_compound","[D4R;$(*(@*)(@*)(@*)@*)]");
		
		// # two different rings sharing exactly two atoms
		queries.put("annelated_rings","[R;$(*(@*)(@*)@*);!$([R2;$(*(@*)(@*)(@*)@*)])]@[R;$(*(@*)(@*)@*);!$([R2;$(*(@*)(@*)(@*)@*)])]");
		
		
		//	# 5 or 6-membered ring containing one O and at least one (r5) or two (r6) oxygen-substituents.		
		queries.put("sugar_pattern_1","[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
 
		// # 5 or 6-membered ring containing one O and an acetal-like bond at postion 2.
		queries.put("sugar_pattern_2","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				 
		// # combination of the two above
		queries.put("sugar_pattern_combi","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_reducing","[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_alpha","[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_beta","[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				
		// # pattern1 occours more than once (in same molecule, but moieties don't have to be adjacent!)
		queries.put("poly_sugar_1","([OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)].[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)])");
				
		// # pattern2 occours more than once (in same molecule, but moieties don't have to be adjacent!)
		queries.put("poly_sugar_2","([OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)].[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)])");
				
		queries.put("conjugated_double_bond","*=*[*]=,#,:[*]");
		queries.put("conjugated_tripple_bond","*#*[*]=,#,:[*]");
		
		//		# only one single-bonded substituent on each DB-atom. no aromats. 
		//		# only found when character of DB is explicitely stated.
//		queries.put("cis_double_bond","*/[D2]=[D2]\*");
		//		# only one single-bonded substituent on each DB-atom. no aromats. 
		//		# only found when character of DB is explicitely stated.		
		queries.put("trans_double_bond:","*/[D2]=[D2]/*");
		
		// # should hits all combinations of two acids
		queries.put("mixed_anhydride","[$(*=O),$([#16,#14,#5]),$([#7]([#6]=[OX1]))][#8X2][$(*=O),$([#16,#14,#5]),$([#7]([#6]=[OX1]))]");
		queries.put("halogen_on_hetero","[FX1,ClX1,BrX1,IX1][!#6]");
		queries.put("halogen_multi_subst","[F,Cl,Br,I;!$([X1]);!$([X0-])]");
		
		// # 1,3 migration of H allowed. Includes keto/enol and amide/enamide. 
		// # Aromatic rings must stay aromatic - no keto form of phenol 
		queries.put("_13-tautomerizable","[$([#7X2,OX1,SX1]=*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=*),$([#7X3,OX2,SX2;!H0]*:n)]");
		queries.put("_15-tautomerizable","[$([#7X2,OX1,SX1]=,:**=,:*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=**=*),$([#7X3,OX2,SX2;!H0]*=,:**:n)]");
		
		// # Hits atoms with tetrahedral chirality, if chiral center is specified in the SMILES string
		// # depictmach does not find oxonium, sulfonium, or sulfoxides!
		queries.put("chiral_center_specified","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*)]");
		
		// # Hits atoms with tetrahedral chirality, if chiral center is not specified in the SMILES string
		// # "@?" (unspecified chirality) is not yet supported in Open Babel Version 2.0 
		queries.put("chiral_center_unspecified","[$([*@?](~*)(~*)(*)*),$([*@?H](*)(*)*),$([*@?](~*)(*)*),$([*@?H](~*)~*)]");
//		queries.put("","");
//		queries.put("","");
//		queries.put("","");
//		queries.put("","");

		queries.put("13_benzodiazole", "c1:[nH]:c(:c2:n:1):c:c:c:c:2");
		queries.put("13_benzothiazole", "c1:s:c(:c2:n:1):c:c:c:c:2");
		queries.put("23-dihydro_14-benzodioxine", "c1:s:c(:c2:n:1):c:c:c:c:2");
		queries.put("arene_epoxide", "[#6]=,:1[#6]=,:[#6][#6]2-[#8]-[#6]2[#6]=,:1");
		queries.put("1_methoxy2_unsubstituted_arene", "[H][#6]1=,:[#6][#6]=,:[#6][#6]=,:[#6]1-[#8;X2][C;X4]([H])([H])[H]");
		
		queries.put("neonicotinoid_pattern1", "[H]C([H])([#7;R1]1-[#6;R1]-[#6;R1]-[#7;R1]-[#6;R1]-1=[#7;X2]/[#7;X3+](-[#8;X1-])=[O;X1])[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](Cl)-[#7]=[#6;R1]-1");
		queries.put("neonicotinoid_pattern2", "[H][#7;X3](-[#6;R1]-1=[#7;R1+](-[#6;R1]-[#6;R1]-[#7;R1]-1)C([H])([H])[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](Cl)-[#7]=[#6;R1]-1)-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("neonicotinoid_pattern3", "[H]C([H])([#7;X3R0]-[#6;R0](-[#7;X3R0]-[#6;R0])=[#7;X2]-[#7;X3+](-[#8;X1-])=[O;X1])[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](Cl)-[#7]=[#6;R1]-1");
		queries.put("neonicotinoid_pattern4", "[H]C([H])([#7;R1]-1-[#6;R1]-[#6;R1]-[#7;R1]([#1,#6;A])-[#6;R1]-1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](Cl)-[#7]=[#6;R1]-1");

		
		
		
		
		
		queries.put("pattern_c51", "ccc(c)O");
		queries.put("pattern_c52", "cccccO");	
		queries.put("pattern_c53", "cc(c)O");	
		queries.put("pattern_c54", "ccccO");	
		queries.put("pattern_c55", "ccccccO");	
		queries.put("pattern_c56", "cccc(O)cc");	
		queries.put("pattern_c57", "ccccc(c)O");	
		queries.put("pattern_c58", "Oc1ccccc1");	
		queries.put("pattern_c59", "cccc(c)O");	
		queries.put("pattern_c60", "ccnc");	
		queries.put("pattern_c61", "cnc");	
		queries.put("pattern_c62", "ccn");	
		queries.put("pattern_c63", "cccnc");	
		queries.put("pattern_c64", "cccn");	
		queries.put("pattern_c65", "ccN");	
		queries.put("pattern_c66", "cccN");	
		queries.put("pattern_c67", "ccC");	
		queries.put("pattern_c68", "Cc1ccccc1");	
		queries.put("pattern_c69", "cccC");
		queries.put("pattern_c70", "ccc(C)cc");	
		queries.put("pattern_c71", "cccc");	
		queries.put("pattern_c72", "ccccC");	
		queries.put("pattern_c73", "ccccc");	
		queries.put("pattern_c74", "cc(c)C");
		queries.put("pattern_c75", "ccc");	
		queries.put("pattern_c76", "c1ccccc1");	
		queries.put("pattern_c77", "cccccc");
		queries.put("pattern_c78", "cccc(C)cc");	
		queries.put("pattern_c79", "ccc(c)C");	
		queries.put("pattern_c80", "cccc(c)C");
		queries.put("pattern_c81", "cccccC");	
		queries.put("pattern_c82", "ccccc(c)C");	
		queries.put("pattern_c83", "ccccccC");
		queries.put("pattern_c84", "CCN");	
		queries.put("pattern_c85", "CNC");	
		queries.put("pattern_c86", "CCC");
		
		// CYP2A6
		queries.put("pattern_c87", "CCNC");	
		queries.put("pattern_c88", "CCCC");	
		queries.put("pattern_c89", "CCO");
		queries.put("pattern_c90", "cccccn");	
		queries.put("pattern_c91", "ccncc");	
		queries.put("pattern_c92", "ccccnc");
		queries.put("pattern_c93", "ccccn");	
		queries.put("pattern_c94", "c1ccncc1");	
		queries.put("pattern_c95", "cccncc");
		
		// CYP2B6
		queries.put("pattern_c96", "CCc1ccccc1");	
		queries.put("pattern_c97", "cccc(cc)CC");	
		queries.put("pattern_c98", "ccccccCC");
		queries.put("pattern_c99", "ccccc(c)CC");	
		queries.put("pattern_c100", "ccc(cc)CC");	
		queries.put("pattern_c101", "cccc(c)CC");
		queries.put("pattern_c102", "cccccCC");	
		queries.put("pattern_c103", "ccc(c)CC");	
		queries.put("pattern_c104", "ccccCC");
		queries.put("pattern_c105", "cc(c)CC");	
		queries.put("pattern_c106", "ccc(O)cc");	
		queries.put("pattern_c107", "NC=O");
		queries.put("pattern_c108", "cCN");	
		queries.put("pattern_c109", "ccCC");	
		queries.put("pattern_c110", "ccCN");

		
		// Veith, H. et al. (2009); Nature Biotechnology; Volume 27, Issue 11, November 2009, Pages 1050-1055
		// Comprehensive characterization of cytochrome P450 isozyme selectivity across chemical libraries
		// removed patterns no: 19, 69, 75, 78, 131, 146, 148, 250, 253, 255, 262, 264, 271 (2_phenylpyrimidine),
		// 272 (3_phenylpyrimidine)
		
		queries.put("pattern_c111", "CCCCCCCCC");
		queries.put("pattern_c112", "CCCCCCCCCC=O");
		queries.put("pattern_c113", "CCCCCCCCC=CCCC");
		queries.put("pattern_c114", "CNCN");
		queries.put("pattern_c115", "CCOC(=O)C");
		queries.put("pattern_c116", "CCCNC(=O)CCC");
		queries.put("pattern_c117", "CNC=CC=O");
		queries.put("pattern_c118", "CCC(N)C(=O)O");
		queries.put("pattern_c119", "CC=C(NC(=O)C)C=O");
		queries.put("pattern_c120", "CCNc1:c:c:c:c:c:1");
		queries.put("pattern_c121", "Cc1:c:c:c:c(N):c:1");
		queries.put("pattern_c122", "Cc1:c:c:c(N):c:c:1");
		queries.put("pattern_c123", "NCCNc1:c:c:c:c:c:1");
		queries.put("pattern_c124", "CCNCNc1:c:c:c:c:c:1");
		queries.put("pattern_c125", "CN(C)CNc1:c:c:c:c:c:1");
		queries.put("pattern_c126", "CCN(C)CNc1:c:c:c:c:c:1");
		queries.put("pattern_c127", "CCCCNCNc1:c:c:c:c:c:1");
		queries.put("pattern_c128", "NCCCNCNc1:c:c:c:c:c:1");
		queries.put("pattern_c129", "NC(=O)CNc1:c:c:c:c:c:1");
		queries.put("pattern_c130", "NC(=S)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c131", "CCNC(=S)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c132", "c1:c:c(N(=O)=O):c:c:c:1");
		queries.put("pattern_c133", "Nc1:c:c:c:c:c:1N(=O)=O");
		queries.put("pattern_c134", "C=Nc1:c:c:c:c:c:1");
		queries.put("pattern_c135", "N=CNc1:c:c:c:c:c:1");
		queries.put("pattern_c136", "NNc1:c:c:c:c:c:1");
		queries.put("pattern_c137", "NCc1:c:c:c:c:c:1");
		queries.put("pattern_c138", "c1:c:c:c(CNC):c:c:1");
		queries.put("pattern_c139", "CCNCc1:c:c:c:c:c:1");
		queries.put("pattern_c140", "CN(C)Cc1:c:c:c:c:c:1");
		queries.put("pattern_c141", "CCCNCc1:c:c:c:c:c:1");
		queries.put("pattern_c142", "O=CNCc1:c:c:c:c:c:1");
		queries.put("pattern_c143", "CCCCCNCc1:c:c:c:c:c:1");
		queries.put("pattern_c144", "CC(=O)NCc1:c:c:c:c:c:1");
		queries.put("pattern_c145", "CC(N)c1:c:c:c:c:c:1");
		queries.put("pattern_c146", "COc1:c:c:c:c(CN):c:1");
		queries.put("pattern_c147", "NCc1:c:c:c:c:c:1O");
		queries.put("pattern_c148", "CN(C)Cc1:c:c:c(O):c:c:1");
		queries.put("pattern_c149", "NCc1:c:c:c:c:c:1Cl");
		queries.put("pattern_c150", "CN(C)Cc1:c:c:c:c:c:1O");
		queries.put("pattern_c151", "NCc1:c:c:n:c:c:1");
		queries.put("pattern_c152", "CNCc1:c:c:n:c:c:1");
		queries.put("pattern_c153", "NCc1:c:c:c:c:n:1");
		queries.put("pattern_c154", "N=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c155", "CN=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c156", "NN=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c157", "CNN=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c158", "O=CNN=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c159", "CC(=O)NN=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c160", "C(c1:c:c:c:c:c:1)=NO");
		queries.put("pattern_c161", "COc1:c:c:c:c:c:1");
		queries.put("pattern_c162", "Cc1:c:c:c(O):c:c:1");
		queries.put("pattern_c163", "COc1:c:c:c(C):c:c:1");
		queries.put("pattern_c164", "CCOc1:c:c:c:c:c:1");
		queries.put("pattern_c165", "COc1:c:c:c:c(C):c:1");
		queries.put("pattern_c166", "Cc1:c(OC):c:c:c:c:1");
		queries.put("pattern_c167", "CCc1:c:c:c(OC):c:c:1");
		queries.put("pattern_c168", "CCCc1:c:c:c(O):c:c:1");
		queries.put("pattern_c169", "NCCc1:c:c:c(O):c:c:1");
		queries.put("pattern_c170", "COc1:c:c:c(C):c:c:1OC");
		queries.put("pattern_c171", "OCCOc1:c:c:c:c:c:1");
		queries.put("pattern_c172", "C(CC)Cc1:c:c:c(O):c:c:1");
		queries.put("pattern_c173", "O=CCOc1:c:c:c:c:c:1");
		queries.put("pattern_c174", "COc1:c:c:c:c(OC):c:1");
		queries.put("pattern_c175", "CC=Cc1:c:c:c(O):c:c:1");
		queries.put("pattern_c176", "CC(O)COc1:c:c:c:c:c:1");
		queries.put("pattern_c177", "CCOc1:c:c:c:c:c:1OC");
		queries.put("pattern_c178", "Oc1:c:c:c:c:c:1C=O");
		queries.put("pattern_c179", "COc1:c:c:c:c:c:1C=O");
		queries.put("pattern_c180", "COc1:c:c:c(C=O):c:c:1OC");
		queries.put("pattern_c181", "Oc1:c:c:c:c(C=C):c:1");
		queries.put("pattern_c182", "COc1:c:c:c(Cl):c:c:1");
		queries.put("pattern_c183", "NCCCOc1:c:c:c:c:c:1");
		queries.put("pattern_c184", "O=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c185", "NC(c1:c:c:c:c:c:1)=O");
		queries.put("pattern_c186", "CNC(c1:c:c:c:c:c:1)=O");
		queries.put("pattern_c187", "CCNC(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c188", "c1:c:c:c(C(C)=O):c:c:1");
		queries.put("pattern_c189", "OC(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c190", "COC(c1:c:c:c:c:c:1)=O");
		queries.put("pattern_c191", "CCC(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c192", "CC=CC(c1:c:c:c:c:c:1)=O");
		queries.put("pattern_c193", "O=Cc1:o:c:c:c:1");
		queries.put("pattern_c194", "O=Cc1:c:c:c:[nH]:1");
		queries.put("pattern_c195", "NC(=O)c1:o:c:c:c:1");
		queries.put("pattern_c196", "Cn1:c:c:c:c:1C=O");
		queries.put("pattern_c197", "O=Cc1:c:c:c:n:c:1");
		queries.put("pattern_c198", "CCNC(=O)c1:o:c:c:c:1");
		queries.put("pattern_c199", "O=CNc1:c:c:c:c:c:1");
		queries.put("pattern_c200", "CC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c201", "NC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c202", "CCCC(Nc1:c:c:c:c:c:1)=O");
		queries.put("pattern_c203", "NCC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c204", "CN(C)C(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c205", "CC(=O)Nc1:c:c:c:c(C):c:1");
		queries.put("pattern_c206", "CC(=O)Nc1:c:c:c:c:c:1C");
		queries.put("pattern_c207", "CC(=O)Nc1:c:c:c(C):c:c:1");
		queries.put("pattern_c208", "CC(C)C(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c209", "CSCC(=O)N");
		queries.put("pattern_c210", "CSCC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c211", "CC=CC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c212", "CCN(C)CC(=O)Nc1:c:c:c:c:c:1");
		queries.put("pattern_c213", "O=CNc1:n:c:c:s:1");
		queries.put("pattern_c214", "O=CNc1:n:n:c:s:1");
		queries.put("pattern_c215", "n1c(N)ccnc1");
		queries.put("pattern_c216", "CNc1:c:c:n:c:n:1");
		queries.put("pattern_c217", "CCNc1:c:c:n:c:n:1");
		queries.put("pattern_c218", "CNc1:n:c:c:c(N):n:1");
		queries.put("pattern_c219", "Cc1:c:c:n:c:c:1");
		queries.put("pattern_c220", "Cc1:c:c:n:c:n:1");
		queries.put("pattern_c221", "Cc1:n:c:c:c(N):n:1");
		queries.put("pattern_c222", "CCn1:c:c:c:c:1");
		queries.put("pattern_c223", "c1:c:c(Nc(:c:c:c:c2):c:2):n:c:n:1");
		queries.put("pattern_c224", "Cc1:c:c:c(C):n:1C");
		queries.put("pattern_c225", "Cc1:c:c:n:c(C):c:1");
		queries.put("pattern_c226", "Cc1:c:c(N):n:c:n:1");
		queries.put("pattern_c227", "Cc1:c:c(C):n:c:n:1");
		queries.put("pattern_c228", "Nc1:c(:c2:n:c:n:1):c:c:c:c:2");
		queries.put("pattern_c229", "CNc1:n:c:n:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c230", "CNc1:n:c(C):n:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c231", "c1:c:c:c(:c:c:1)c2:n:c:c3:c:c:c:c:c:3:n:2");
		queries.put("pattern_c232", "CNc1:n:c(:n:c2:c:c:c:c:c:1:2)c3:c:c:c:c:c:3");
		queries.put("pattern_c233", "CNc1:n:c:n:c2:c:c:c(:c:c:1:2)c3:c:c:c:c:c:3");
		queries.put("pattern_c234", "c1:c:c(:c2:c:c:1):c(N(C)C):n:c:n:2");
		queries.put("pattern_c235", "Cc1:c:c:c2:c:c:c:c:c:2:n:1");
		queries.put("pattern_c236", "CCCc1:c:c:c:c:c:1");
		queries.put("pattern_c237", "C=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c238", "CC(C)Cc1:c:c:c:c:c:1");
		queries.put("pattern_c239", "c1:c:c:c(C=CC):c:c:1");
		queries.put("pattern_c240", "CCc1:c:c:c:c:n:1");
		queries.put("pattern_c241", "CC(C)(C)c1:c:c:c:c:c:1");
		queries.put("pattern_c242", "CCCn1:c:n:c:c:1");
		queries.put("pattern_c243", "CCCCn1:c:c:n:c:1");
		queries.put("pattern_c244", "NCCc1:c:c:c:c:c:1");
		queries.put("pattern_c245", "CNCCc1:c:c:c:c:c:1");
		queries.put("pattern_c246", "CN(CCc1:c:c:c:c:c:1)C");
		queries.put("pattern_c247", "C(Cc1:c:c:c:c:c:1)=O");
		queries.put("pattern_c248", "NC(=O)Cc1:c:c:c:c:c:1");
		queries.put("pattern_c249", "CNC(Cc1:c:c:c:c:c:1)=O");
		queries.put("pattern_c250", "CCNC(=O)Cc1:c:c:c:c:c:1");
		queries.put("pattern_c251", "OCCCc1:c:c:c:c:c:1");
		queries.put("pattern_c252", "O=CC=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c253", "CCOCc1:c:c:c:c:c:1");
		queries.put("pattern_c254", "OC(=O)CCc1:c:c:c:c:c:1");
		queries.put("pattern_c255", "CC(O)c1:c:c:c:c:c:1");
		queries.put("pattern_c256", "CN(C)CCCc1:c:c:c:c:c:1");
		queries.put("pattern_c257", "NCC=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c258", "NC(=O)C=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c259", "NC(CCc1:c:c:c:c:c:1)=O");
		queries.put("pattern_c260", "CCC(O)c1:c:c:c:c:c:1");
		queries.put("pattern_c261", "OC(=O)Cc1:c:c:c:c:c:1");
		queries.put("pattern_c262", "OC(=O)C=Cc1:c:c:c:c:c:1");
		queries.put("pattern_c263", "CCN1CCCCC1");
		queries.put("pattern_c264", "O=CN1CCCCC1");
		queries.put("pattern_c265", "CC(=O)N1CCCCC1");
		queries.put("pattern_c266", "CCN1CCNCC1");
		queries.put("pattern_c267", "CCN1CCN(C)CC1");
		queries.put("pattern_c268", "CCN1CCOCC1");
		queries.put("pattern_c269", "CCCN1CCCCC1");
		queries.put("pattern_c270", "CN1CCN(CC1)C(=O)C");
		queries.put("pattern_c271", "CCCN1CCN(C)CC1");
		queries.put("pattern_c272", "CCC1CCCCN1");
		queries.put("pattern_c273", "CCCN1CCC(C)CC1");
		queries.put("pattern_c274", "CNCCN1CCOCC1");
		queries.put("pattern_c275", "CC1CCCC(C)N1");
		queries.put("pattern_c276", "OCC1OCC=CC1O");
		queries.put("pattern_c277", "O=CC1CCCCC1");
		queries.put("pattern_c278", "NC(=O)C1CCCCC1");
		queries.put("pattern_c279", "CNC1CCCCC1");
		queries.put("pattern_c280", "CC1CCC(CC1)C(=O)N");
		queries.put("pattern_c281", "CC(=O)NC1CCCCC1");
		queries.put("pattern_c282", "O=C1NC(=O)C2CCCCC12");
		queries.put("pattern_c283", "OC1CC(=N)C2CCC3C(C2C1O)C(=O)NC3=O");
		queries.put("pattern_c284", "CN1C(=O)C2CCC3C(C(O)C(O)CC3=N)C2C1=O");
		queries.put("pattern_c285", "CC12CCCCC1CCC2");
		queries.put("pattern_c286", "O=S(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c287", "NS(=O)(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c288", "CNS(=O)(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c289", "CCNS(=O)(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c290", "CCN(CC)S(=O)(=O)c1:c:c:c:c:c:1");
		queries.put("pattern_c291", "Cc1:c:c:c(:c:c:1)S(=O)=O");
		queries.put("pattern_c292", "S(N1CCCCC1)(=O)=O");
		queries.put("pattern_c293", "Nc1:c:c:c(:c:c:1)S(=O)=O");
		queries.put("pattern_c294", "Nc1:c:c:c(:c:c:1)S(=O)(=O)N");
		queries.put("pattern_c295", "O=S(=O)(N1CCCCC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c296", "O=S(=O)(Nc1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c297", "Sc1:c:c:c:c:c:1");
		queries.put("pattern_c298", "c1:c:c:c(SC):c:c:1");
		queries.put("pattern_c299", "SCc1:c:c:c:c:c:1");
		queries.put("pattern_c300", "CSCc1:c:c:c:c:c:1");
		queries.put("pattern_c301", "CCSc1:c:c:c:c:c:1");
		queries.put("pattern_c302", "CCSCc1:c:c:c:c:c:1");
		queries.put("pattern_c303", "CCSc1:n:n:c(C):n:1C");
		queries.put("pattern_c304", "Cc1:n:n:c(SCC=O):n:1C");
		queries.put("pattern_c305", "Clc1:c:c:c:c:c:1");
		queries.put("pattern_c306", "Cc1:c:c:c:c:c:1Cl");
		queries.put("pattern_c307", "Cc1:c:c:c(Cl):c:c:1");
		queries.put("pattern_c308", "Cc1:c:c:c:c(Cl):c:1");
		queries.put("pattern_c309", "FC(F)(F)c1:c:c:c:c:c:1");
		queries.put("pattern_c310", "c1:c:c:c2:[nH]:c:c:c:2:c:1/c1:c:c2:c:c:c:c:c:2:n:1");
		queries.put("pattern_c311", "Cc1:c:[nH]:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c312", "CCc1:c:[nH]:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c313", "Cn1:c:c:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c314", "Cc1:[nH]:c(:c2:c:1):c:c:c:c:2");
		queries.put("pattern_c315", "Cc1:n:c2:c:c:c:c:c:2:[nH]:1");
		queries.put("pattern_c316", "Cn1:c(:c2:n:c:1):c:c:c:c:2");
		queries.put("pattern_c317", "Cn1:c:n:c2:c:n:c:n:c:1:2");
		queries.put("pattern_c318", "c1:c:n:c2:c:c:c:c:c:2:c:1");
		queries.put("pattern_c319", "c1:c:c:c2:c:c:c:c:c:2:c:1");
		queries.put("pattern_c320", "Cc1:c:c:n:c2:c:c:c:c:c:1:2");
		queries.put("pattern_c321", "C1Oc2:c:c:c:c:c:2O1");
		queries.put("pattern_c322", "Cc1:c:c:c2OCOc:2:c:1");
		queries.put("pattern_c323", "CN1CCc2:c:c:c:c:c:2C1");
		queries.put("pattern_c324", "C(c1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c325", "c1:c:c:c(:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c326", "C(Nc1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c327", "c1:c:c:c(:c:c:1)c2:n:c:c:c:n:2");
		queries.put("pattern_c328", "c1:c:c:c(:c:c:1)c2:c:c:n:c:n:2");
		queries.put("pattern_c329", "C(Cc1:c:c:c:c:c:1)NCc2:c:c:c:c:c:2");
		queries.put("pattern_c330", "O=C(Nc1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c331", "N(c1:c:c:c:c:c:1)c2:c:c:n:c:n:2");
		queries.put("pattern_c332", "c1:c(CNCc2:c:c:c:c:c:2):c:c:c:c:1");
		queries.put("pattern_c333", "C(Oc1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c334", "CC(c1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c335", "c1:c:c:c(:c:c:1)n2:c:c:c:n:2");
		queries.put("pattern_c336", "C(Nc1:c:c:n:c:n:1)c2:c:c:c:s:2");
		queries.put("pattern_c337", "Nc1:c:c(:n:c:n:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c338", "COc1:c:c:c:c:c:1CNc2:c:c:n:c:n:2");
		queries.put("pattern_c339", "CNc1:c:c(:n:c:n:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c340", "O=C(Cc1:c:c:c:c:c:1)Nc2:c:c:c:c:c:2");
		queries.put("pattern_c341", "C(Cc1:c:c:c:c:c:1)c2:c:c:c:c:c:2");
		queries.put("pattern_c342", "c1:c:c(Nc2:c:c:c:c:c:2):c:c:c:1");
		queries.put("pattern_c343", "CN(Cc1:c:c:c:c:c:1)Cc2:c:c:c:c:c:2");
		queries.put("pattern_c344", "C(N1CCCCC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c345", "O=C(N1CCCCC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c346", "C1CN(CCN1)c2:c:c:c:c:c:2");
		queries.put("pattern_c347", "CN1CCN(CC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c348", "CCN1CCN(CC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c349", "NC(C1CC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c350", "CN1CCN(Cc2:c:c:c:c:c:2)CC1");
		queries.put("pattern_c351", "O=C(Nc1:c:c:c:c:c:1)N2CCCCC2");
		queries.put("pattern_c352", "C1CCN(CC1)c2:c:c:c:c:c:2");
		queries.put("pattern_c353", "O=C1CCCN1c2:c:c:c:c:c:2");
		
// // added on January 17, 2017
//		Structural Alerts of Mutagens and Carcinogens: Current Computer-Aided Drug Design, 2006, Vol. 2, No. 2
//				
//		Alkyl phosphonate
//		Alkyl sulfonates
//		Alkyl hydrazines
//		Alkyl aldehydes
//		N-methylol derivatives
//		2-chloroethyl mustards
//		N-chloroamines
//		propiolactones
//		propiosultones
//		aromatic aziridinyl derivatives
//		aliphatic aziridinyl derivatives
//		alkyl N-nitrosoamines
//		cyclic N-nitrosamines
		
		
		queries.put("aldoxime", "[H][#8;X2]-[#7;X2]=[#6;X3](/[H])-[#6]");
		queries.put("alkylaldoxime", "[H][#8;X2]-[#7;X2]=[#6;X3]([H])-[#6;X4]"); // Boucher JL et al. PMID: 8011645, Biochemistry. 1994 Jun 28;33(25):7811-8
		queries.put("arylaldoxime", "[H][#8;X2]-[#7;X2]=[#6;X3]([H])-[#6;a]"); // Boucher JL et al. PMID: 8011645, Biochemistry. 1994 Jun 28;33(25):7811-8
		queries.put("ketoxime", "[H][#8;X2]-[#7;X2]=[#6;X3](/[#6])-[#6]");

		// Kalgutkar AS. et al (PMID: 15975040; Current Drug Metabolism, 2005, 6, 161-225)
		queries.put("n_alkylformamide", "[H][#7;X3](-[#6;X4])-[#6;X3]([H])=[O;X1]");
		queries.put("n_phenylformamide", "[H][#7;X3](-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#6;X3]([H])=[O;X1]");
		queries.put("dialkylthiourea", "[H][#7;X3](-[#6;X4])-[#6;X3](=[S;X1])-[#7;X3]([H])-[#6;X4]");
		queries.put("alkylarylthiourea", "[H][#7;X3](-[#6;a])-[#6;X3](=[S;X1])-[#7;X3]([H])-[#6;X4]");
		queries.put("diarylthiourea", "[H][#7;X3](-[#6;a])-[#6;X3](=[S;X1])-[#7;X3]([H])-[#6;a]");
		queries.put("masqued_o_quinone", "[#6]-[#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1-[#8;X2]-[#6]");
		queries.put("masqued_p_quinone", "[#6]-[#8;X2]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8;X2]-[#6]");
		queries.put("catechol_monoether", "[H][#8;X2]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1-[#8;X2]-[#6]");
		queries.put("p_quinone_monoether", "[H][#8;X2]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#8;X2]-[#6]");
		queries.put("p_cresol", "[H][#6][C;R1]=,:1([H])-,:[#6;R1]=,:[#6;R1]-,:[C;R1]([H])(=,:[#6;R1][#6;R1]=,:1)[#8][H]");
		queries.put("13_benzodioxole_bicyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2-[#8;R1]-[#6;R1]-[#8;R1]-[#6;R2]=,:2[#6;R1]=,:1");
		queries.put("13_benzodioxole", "[#6]=,:1[#6]=,:[#6][#6]=,:2-[#8]-[#6]-[#8]-[#6]=,:2[#6]=,:1");

		
		queries.put("5_alkoxyindole", "[#1,$([CX4][H])]-[#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]-2=,:[#6;R2]([#6;R1]=,:1)-[#6;R]=[#6;R1]-[#7;R1]-2-[#1,$([CX4][H])]");
		queries.put("3_Methylindole", "[H]C([H])([H])[#6;R1]-1=[#6;R1]-[#7;R1]-[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]-1=,:2");
		queries.put("a_aroylthiophene", "[#6]-[#6;X3](=O)-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#16;R1]1");
		queries.put("dibenzazepine", "[#6]=,:1[#6]=,:[#6][#6]=,:2[#7][#6]=,:3[#6]=,:[#6][#6]=,:[#6][#6]=,:3[#6]=,:[#6][#6]=,:2[#6]=,:1");
		queries.put("acyl_halide", "[#1,#6]-[#6;X3](-[F,Cl,Br,I])=[O;X1]");

		// F. Peter Guengerich.Guengerich Chem. Res. Toxicol., Vol. 14, No. 6, 2001
		queries.put("pyrrolizidine_alkaloid", "[H]C([H])([#8;X2]-[#1,#6])[#6;R1]-1=[#6;R1]-[#6;R1]-[#7;X3R2]-2-[#6;R1]-[#6;R1]-[#6;R1](-[#8;X2]-[#1,#6])-[#6;R2]-1-2");
		queries.put("14_dihydropyridine", "[H]C-,:1([#6;X4])-,:[#6;R1]=,:[#6;R1][#7][#6;R1]=,:[#6;R1]-,:1");	
		queries.put("cyclohexene", "*-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("alkyl_aldehyde_pattern", "[H][#6;X3](=[O;X1])[C;X2]([H])([H])C([H])([H])[#1,#6]");
		queries.put("23_dihydro_1H_indol_3_one", "[H]C-,:1([H])-,:[#7][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6]-,:1=O");
		queries.put("fatty_acid_hydroperoxide", "[H]-[#6;X3](-[#6])=[#6;X3](/[H])-[#6;X3](-[H])=[#6;X3](/[H])[C;X4]([H])([#6])[#6]([#8;A;X2H1,X1-])=[O;X1]");
		
		// Andrew Parkinson: BIOTRANSFORMATION OF XENOBIOTICS; http://www.farmasi.unud.ac.id/ind/wp-content/uploads/Bio-Transformation-of-Xenobiotics.pdf
		queries.put("thiophosphoric_acid_ester", "[#6]-[#8;X2][P;X4](=[S;X1])([#8;X2]-[#1,#6])[#8;X2]-[#1,#6]");
		queries.put("ester_deriv", "[#6;!$(C=[O,N,S])]-[#8;X2][#6;A;X3]([#7X3,#6])=[O;X1]");
		queries.put("4_hydroxyanilide", "[H][#8]!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6](-[#6,#1])=[O;X1]");
		queries.put(" 4_monosubstituted_14_dihydropyridine", "[H][#7;R1]1[#6;R1]=,:[#6;R1]-,:[C;R1]([H])(*)-,:[#6;R1]=,:[#6;R1]1");
		queries.put("n_aminopyrrolidine_derivative", "[H][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6]-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#7;R1]-1");	
		queries.put("n_organo_substituted_pyrrolidine", "[#6]-[#7;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-1");	
		queries.put("s_alkyl_group", "[#6;X4]-[#16;X2]");
		queries.put("112_trisubstituted_alkene", "[H]-[#6;X3](-[#6])=[#6;X3](/[#6])-[#6]");
		queries.put("12_disubstituted_alkene", "[H]-[#6;X3](-[#6])=[#6;X3](/[H])-[#6]");
		queries.put("desaturation_precursor", "[H][#6;X3](-[#6])-[#6;X3]([H])-[#6]");
		queries.put("aromatic_methyl", "[H]C=,:1([#6;A;H3X4])-,:[#6]=,:[#6][#6]=,:[#6][#6]=,:1");

		// B. Clement; Oxidation and reduction of nitrogen via CYP450: Importance of the reduction of genotoxic N-hydroxylated functional groups
		// From the book Drug Metabolism: Towards the Next Millennium by N. J. Gooderham
		queries.put("primary_amino_azaheterocycle", "[$([*;R][#7;R]=,:[#6;R1][#7;A]([H])[H]),$([*;R]=,:[#7;R][#6;R1][#7;A]([H])[H])]");
		queries.put("amino_azaheterocycle", "[$([*;R][#7;R]=,:[#6;R1][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]),$([*;R]=,:[#7;R][#6;R1][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])])]");
		queries.put("tertiary_heteroaromatic", "[*;R][#7;X2R]=,:[*;R]");
		queries.put("n_arylamide", "[H][#7;X3](-[*;a])[#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("imine", "[H][#7;A;X2]=[#6;A;X3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])]");
		queries.put("amidinohydrazone", "[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])](=[#7;A;X2;!$(NC=[O,S])])[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#7;X2]=[#6]");

		// Clement, B. (2002) Reduction of N-hydroxylated compounds: Amidoximes (N-hydroxyamidines) as pro-drugs of amidines
		queries.put("benzamidoxime", "[H][#8;X2]-[#7;X2]=[#6;X3](-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)-[#7;X3]([H])[H]");

		//TERRENCE R. BURKE, JR.(1980) MECHANISM OF DEFLUORINATION OF ENFLURANE; DRUG METABOLISM AND DIsPosITIoN; Vol.9, No.1	
		// KENNETH E. THUMMEL (1993); HUMAN  LIVER MICROSOMAL ENFLURANE DEFLUORINATION CATALYZED BYnCYTOCHROME P-450 2E1
		
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
		
		// HAO CHEN et al.(2003) METABOLISM OF (S)-5,6-DIFLUORO-4-CYCLOPROPYLETHYNYL-4-TRIFLUOROMETHYL- 3, 4-DIHYDRO-2(1H)-QUINAZOLINONE, A NON-NUCLEOSIDE REVERSE TRANSCRIPTASE INHIBITOR, IN HUMAN LIVER MICROSOMES. METABOLIC ACTIVATION AND ENZYME KINETICS
		// DRUG METABOLISM AND DISPOSITION; Vol. 31, No. 1
		queries.put("organofluoride_pattern1", "F[#6][#6]F");
		
		
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
		
		// Kharasch  et al (1993). Identification of cytochrome P450 2E1 as the predominant enzyme catalyzing human liver microsomal defluorination of sevoflurane, isoflurane, and methoxyflurane 
		queries.put("organofluoride_pattern2", "[H][#6;X4](F)-[#6;X4](-[F,Cl,Br,I])-[F,Cl,Br,I]");  // to [F,Cl,Br,I]-[#6;X4](-[F,Cl,Br,I])-[#6;X3]=O
		queries.put("12_dihalo_arene", "[*,#1]-[#6]1=,:[#6][#6]=,:[#6](-[F,Cl,Br,I])[#6](-[F,Cl,Br,I])=,:[#6]1");	
		queries.put("2_chloromethoxy_111_trifluoroethane_derivative", "[H][#6;X4](-[#8;X2]-[#6;X4]-[F,Cl,Br,I])[C;X4](F)(F)F");
		
// 		Clinical Isoflurane Metabolism by Cytochrome P450 2E1; Anesthesiology 3 1999, Vol.90, 766-771.
//		isoflurane to trifluoroacetic acid
		
		
		// Kharasch E.D et al.(2000) Human halothane metabolism, lipid peroxidation, and cytochromes P4502A6 and P4503A4; E J Clin Pharmacol (2000) 55: 853. doi:10.1007/s002280050707		
		// halothane oxidation to trifluoroacetic acid (CYP2E1)	
 		// halothane reduction to chlorotrifluoroethane (CYP2A6)
		queries.put("hydrochlorofluorocarbon_pattern1", "[H][C;X4]([F,Cl,Br,I])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("hydrochlorofluorocarbon_pattern2", "[H][C;X4]([H])([F,Cl,Br,I])[C;X4](F)(F)[F,Cl,Br,I]");
		queries.put("hydrochlorofluorocarbon_pattern3", "[H][#6;X2][C;X4]([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]");
		
		// Jayaprakasam Bolleddula, Kevin DeMent, James P. Driscoll, Philip Worboys, Patrick J. Brassil & David L. Bourdet (2014)
		// Biotransformation and bioactivation reactions of alicyclic amines in drug molecules, Drug Metabolism Reviews, 46:3, 379-419,
		// DOI: 10.3109/03602532.2014.924962
		
		queries.put("azepane", "C1CCCNCC1");
		queries.put("morpholine", "C1COCCN1");
		queries.put("pyrrolidine", "C1CCNC1");
		queries.put("azetidine", "C1CNC1");
		queries.put("aziridine", "C1CN1");
		queries.put("dibenzo_14_diazepine", "[$([#6]=,:1[#6]=,:[#6][#6]=,:2[#7][#6]3=,:[#6][#6]=,:[#6][#6]=,:[#6]3[#6]=,:[#7][#6]=,:2[#6]=,:1),$([#6]=,:1[#6]=,:[#6][#6]=,:2[#7][#6]3=,:[#6][#6]=,:[#6][#6]=,:[#6]3=,:[#6][#7][#6]=,:2[#6]=,:1)]");
		queries.put("nicotine", "[#6]-[#7]-1-[#6]-[#6]-[#6]-[#6]-1!@[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("benzo_14_diazepine", "[$([#6]1=,:[#6][#6]=,:[#6]2[#7][#6]=,:[#6][#7]=,:[#6][#6]2=,:[#6]1),$(c1ccc2nc3ccccc3ncc2c1)]");
		queries.put("5_phenyl_23_dihydro_1h_14-benzodiazepine_derivative", "[*,#1]-[#7;R1]-1-[#6;R2]-2=[#6;R2](-[#6;R1]=[#6;R1](-[*,#1])-[#6;R1]=[#6;R1]-2)-[#6;R1](=[#7;R1]-[#6;R1]-[#6;R1]-1=*)-[#6;R1]-1=[#6;R1](-[*,#1])-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");				
		queries.put("34_dihydro_2H_pyrrole", "C1CC=NC1");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
//		queries.put("", "");
		
		
		// MACCS 322
		
		queries.put("maccs_322_001", "*~*~*~*");
		queries.put("maccs_322_002", "[!#1!#6]");	
//		queries.put("maccs_322_003", "");	
		queries.put("maccs_322_004", "*~*~*~*~*");	
		queries.put("maccs_322_005", "[$(*~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6])]");
		queries.put("maccs_322_006", "[$(*~[!#1!#6]~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6]~[!#1!#6])]");	
		queries.put("maccs_322_007", "[!#1!#6][H]");	
		queries.put("maccs_322_008", "[$(*~*~[#6]([H])[H]),$(*~C(~*)([H])[H])]");		
		queries.put("maccs_322_009", "*~C([H])([H])[H]");
//		queries.put("maccs_322_010", "[F,Cl,Br,I]");	//halogen (already here)
		queries.put("maccs_322_011", "[$(*-*(-*)-*)]");	 // [$(),$(*-*(-*)-*)]
//		queries.put("maccs_322_012", "");	
		queries.put("maccs_322_013", "[$(*@-*(@-*)@-*)]"); // [$(),$(*@-*(@-*)@-*)]
		queries.put("maccs_322_014", "*@-*!@-*@-*");	
		queries.put("maccs_322_015", "*!:*:*!:*");	
		queries.put("maccs_322_016", "*!@-*!@-*");
		queries.put("maccs_322_017", "*!@-*@-*!@-*");
		queries.put("maccs_322_018", "*:*!:*:*");	
//		queries.put("maccs_322_019", "");	// Same as organoheterocycle
		queries.put("maccs_322_020", "[$(*-*(-*)(-*)(-*)-*),$(*@-*(@-*)(@-*)(@-*)@-*),$([!#1!#6!#7!#8!#16!#9!#17!#35!#53])]");	
		queries.put("maccs_322_021", "[!+0]");
//		queries.put("maccs_322_022", ""); // Same as nitrogen_atom	
//		queries.put("maccs_322_023", ""); // Same as sulfur_atom	
//		queries.put("maccs_322_024", ""); // Same as oxygen_atom
//		queries.put("maccs_322_025", "");
//		queries.put("maccs_322_026", "");	
		
		queries.put("maccs_322_027", "C(CC)");	
		queries.put("maccs_322_028", "C(CCC)");	
		queries.put("maccs_322_029", "C(CN)");
		queries.put("maccs_322_030", "C(CCN)");	
		queries.put("maccs_322_031", "C(NN)");	
		queries.put("maccs_322_032", "C(NNC)");
		queries.put("maccs_322_033", "C(NNN)");	
		queries.put("maccs_322_034", "C(CO)");	
		queries.put("maccs_322_035", "C(CCO)");
		queries.put("maccs_322_036", "C(NO)");
		queries.put("maccs_322_037", "C(NCO)");	
		queries.put("maccs_322_038", "C(NNO)");	
		queries.put("maccs_322_039", "C(OO)");	
		queries.put("maccs_322_040", "C(COO)");
		queries.put("maccs_322_041", "C(NOO)");	
		queries.put("maccs_322_042", "C(OOO)");	
		queries.put("maccs_322_043", "[!#1!#6](CC)");
		queries.put("maccs_322_044", "[!#1!#6](CCC)");	
		queries.put("maccs_322_045", "[!#1!#6](CN)");	
		queries.put("maccs_322_046", "[!#1!#6](CCN)");
		queries.put("maccs_322_047", "[!#1!#6](NN)");
		queries.put("maccs_322_048", "[!#1!#6](CNN)");	
		queries.put("maccs_322_049", "[!#1!#6](NNN)");	
		queries.put("maccs_322_050", "[!#1!#6](CO)");	
		queries.put("maccs_322_051", "[!#1!#6](CCO)");
		queries.put("maccs_322_052", "[!#1!#6](NO)");	
		queries.put("maccs_322_053", "[!#1!#6](CNO)");	
		queries.put("maccs_322_054", "[!#1!#6](NNO)");
		queries.put("maccs_322_055", "[!#1!#6](OO)");	
		queries.put("maccs_322_056", "[!#1!#6](COO)");	
		queries.put("maccs_322_057", "[!#1!#6](NOO)");
		queries.put("maccs_322_058", "[!#1!#6](OOO)");
		queries.put("maccs_322_059", "C-C");	
		queries.put("maccs_322_060", "C-N");	
		queries.put("maccs_322_061", "C-O");	
		queries.put("maccs_322_062", "C-S");
		queries.put("maccs_322_063", "C-[#17]");	
		queries.put("maccs_322_064", "C-P");	
		queries.put("maccs_322_065", "C-F");
		queries.put("maccs_322_066", "C-[Br]");
		queries.put("maccs_322_067", "C-[Si]");
		queries.put("maccs_322_068", "C-I");
		queries.put("maccs_322_069", "C-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_070", "N-N");
		queries.put("maccs_322_071", "N-O");
		queries.put("maccs_322_072", "N-S");
		queries.put("maccs_322_073", "N-[#17]");
		queries.put("maccs_322_074", "N-P");
		queries.put("maccs_322_075", "N-F");
		queries.put("maccs_322_076", "N-[Br]");
		queries.put("maccs_322_077", "N-[Si]");
		queries.put("maccs_322_078", "N-I");
		queries.put("maccs_322_079", "N-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_080", "O-O");
		queries.put("maccs_322_081", "O-S");
		queries.put("maccs_322_082", "O-[#17]");
		queries.put("maccs_322_083", "O-P");
		queries.put("maccs_322_084", "O-F");
		queries.put("maccs_322_085", "O-[Br]");
		queries.put("maccs_322_086", "O-[Si]");
		queries.put("maccs_322_087", "O-I");
		queries.put("maccs_322_088", "O-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_089", "S-S");
		queries.put("maccs_322_090", "S-[#17]");
		queries.put("maccs_322_091", "S-P");
		queries.put("maccs_322_092", "S-F");
		queries.put("maccs_322_093", "S-[Br]");
		queries.put("maccs_322_094", "S-[Si]");
		queries.put("maccs_322_095", "S-I");
		queries.put("maccs_322_096", "S-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_097", "[#17]-[#17]");
		queries.put("maccs_322_098", "[#17]-P");
		queries.put("maccs_322_099", "[#17]-F");
		queries.put("maccs_322_100", "[#17]-[Br]");
		queries.put("maccs_322_101", "[#17]-[Si]");
		queries.put("maccs_322_102", "[#17]-I");
		queries.put("maccs_322_103", "[#17]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_104", "P-P");
		queries.put("maccs_322_105", "P-F");
		queries.put("maccs_322_106", "P-[Br]");
		queries.put("maccs_322_107", "P-[Si]");
		queries.put("maccs_322_108", "P-I");
		queries.put("maccs_322_109", "P-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_110", "F-F");
		queries.put("maccs_322_111", "F-[Br]");
		queries.put("maccs_322_112", "F-[Si]");
		queries.put("maccs_322_113", "F-I");
		queries.put("maccs_322_114", "F-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_115", "[Br]-[Br]");
		queries.put("maccs_322_116", "[Br]-[Si]");
		queries.put("maccs_322_117", "[Br]-I");
		queries.put("maccs_322_118", "[Br]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_119", "[Si]-[Si]");
		queries.put("maccs_322_120", "[Si]-I");
		queries.put("maccs_322_121", "[Si]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_122", "I-I");
		queries.put("maccs_322_123", "I-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_124", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_125", "C=C");
		queries.put("maccs_322_126", "C=N");
		queries.put("maccs_322_127", "C=O");
		queries.put("maccs_322_128", "C=S");
		queries.put("maccs_322_129", "C=[#17]");
		queries.put("maccs_322_130", "C=P");
		queries.put("maccs_322_131", "C=F");
		queries.put("maccs_322_132", "C=[Br]");
		queries.put("maccs_322_133", "C=[Si]");
		queries.put("maccs_322_134", "C=I");
		queries.put("maccs_322_135", "C=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_136", "N=N");
		queries.put("maccs_322_137", "N=O");
		queries.put("maccs_322_138", "N=S");
		queries.put("maccs_322_139", "N=O");
		queries.put("maccs_322_140", "N=[#17]");
		queries.put("maccs_322_141", "N=P");
		queries.put("maccs_322_142", "N=F");
		queries.put("maccs_322_143", "N=[Br]");
		queries.put("maccs_322_144", "N=[Si]");
		queries.put("maccs_322_145", "N=I");
		queries.put("maccs_322_146", "N=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_147", "O=O");
		queries.put("maccs_322_148", "O=S");
		queries.put("maccs_322_149", "O=[#17]");
		queries.put("maccs_322_150", "O=P");
		queries.put("maccs_322_151", "O=[Br]");
		queries.put("maccs_322_152", "O=[Si]");
		queries.put("maccs_322_153", "O=I");
		queries.put("maccs_322_154", "O=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_155", "S=S");
		queries.put("maccs_322_156", "S=[#17]");
		queries.put("maccs_322_157", "S=P");
		queries.put("maccs_322_158", "S=F");
		queries.put("maccs_322_159", "S=[Br]");
		queries.put("maccs_322_160", "S=[Si]");
		queries.put("maccs_322_161", "S=I");
		queries.put("maccs_322_162", "S=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_163", "[#17]=[#17]");
		queries.put("maccs_322_164", "[#17]=P");
		queries.put("maccs_322_165", "[#17]=[#9]");
		queries.put("maccs_322_166", "[#17]=[Br]");
		queries.put("maccs_322_167", "[#17]=[Si]");
		queries.put("maccs_322_168", "[#17]=I");
		queries.put("maccs_322_169", "[#17]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_170", "P=P");
		queries.put("maccs_322_171", "P=[#9]");
		queries.put("maccs_322_172", "P=[Br]");
		queries.put("maccs_322_173", "P=[Si]");
		queries.put("maccs_322_174", "P=I");
		queries.put("maccs_322_175", "P=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_176", "[#9]=[#9]");
		queries.put("maccs_322_177", "[#9]=[Br]");
		queries.put("maccs_322_178", "[#9]=[Si]");
		queries.put("maccs_322_179", "[#9]=I");
		queries.put("maccs_322_180", "[#9]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_181", "[Br]=[Br]");
		queries.put("maccs_322_182", "[Br]=[Si]");
		queries.put("maccs_322_183", "[Br]=I");
		queries.put("maccs_322_184", "[Br]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_185", "[Si]=[Si]");
		queries.put("maccs_322_186", "[Si]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_187", "I-I");
		queries.put("maccs_322_188", "I=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_189", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_190", "C#C");
		queries.put("maccs_322_191", "C#N");
		queries.put("maccs_322_192", "C#O");
		queries.put("maccs_322_193", "C#S");
		queries.put("maccs_322_194", "C#S");
		queries.put("maccs_322_195", "C#[#17]");		
		queries.put("maccs_322_196", "C#P");
		queries.put("maccs_322_197", "C#F");
		queries.put("maccs_322_198", "C#[Br]");
		queries.put("maccs_322_199", "C#[Si]");			
		queries.put("maccs_322_200", "C#I");
		queries.put("maccs_322_201", "C#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_202", "N#N");
		queries.put("maccs_322_203", "N#O");	
		queries.put("maccs_322_204", "N#S");	
		queries.put("maccs_322_205", "N#[#17]");	
		queries.put("maccs_322_206", "N#P");	
		queries.put("maccs_322_207", "N#F");	
		queries.put("maccs_322_208", "N#[Br]");	
		queries.put("maccs_322_209", "N#[Si]");	
		queries.put("maccs_322_210", "N#I");
		queries.put("maccs_322_211", "N#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");	
		queries.put("maccs_322_212", "O#O");	
		queries.put("maccs_322_213", "O#S");	
		queries.put("maccs_322_214", "O#[#17]");	
		queries.put("maccs_322_215", "O#P");	
		queries.put("maccs_322_216", "O#F");	
		queries.put("maccs_322_217", "O#[Br]");	
		queries.put("maccs_322_218", "O#[Si]");	
		queries.put("maccs_322_219", "O#I");	
		queries.put("maccs_322_220", "O#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");	
		queries.put("maccs_322_221", "S#S");		
		queries.put("maccs_322_222", "S#[#17]");		
		queries.put("maccs_322_223", "S#P");
		queries.put("maccs_322_224", "S#[#9]");
		queries.put("maccs_322_225", "S#[Br]");
		queries.put("maccs_322_226", "S#[Si]");
		queries.put("maccs_322_227", "S#I");
		queries.put("maccs_322_228", "S#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_229", "[#17]#[#17]");
		queries.put("maccs_322_230", "[#17]#[#15]");
		queries.put("maccs_322_231", "[#17]#[#9]");
		queries.put("maccs_322_232", "[#17]#[Br]");
		queries.put("maccs_322_233", "[#17]#[Si]");
		queries.put("maccs_322_234", "[#17]#I");
		queries.put("maccs_322_235", "[#17]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_236", "[#15]#[#15]");
		queries.put("maccs_322_237", "[#15]#[#9]");
		queries.put("maccs_322_238", "[#15]#[Br]");
		queries.put("maccs_322_239", "[#15]#[Si]");
		queries.put("maccs_322_240", "[#15]#[#53]");
		queries.put("maccs_322_241", "[#15]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_242", "[#9]#[#9]");
		queries.put("maccs_322_243", "[#9]#[Br]");
		queries.put("maccs_322_244", "[#9]#[Si]");
		queries.put("maccs_322_245", "[#9]#[#53]");
		queries.put("maccs_322_246", "[#9]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_247", "[#35]#[Br]");
		queries.put("maccs_322_248", "[#35]#[Si]");
		queries.put("maccs_322_249", "[#35]#[#53]");
		queries.put("maccs_322_250", "[#35]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_251", "[Si]#[Si]");		
		queries.put("maccs_322_252", "[Si]#[#53]");
		queries.put("maccs_322_253", "[Si]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_254", "[#53]#[#53]");
		queries.put("maccs_322_255", "[#53]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_256", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_257", "[#6]@-[#6]");		
		queries.put("maccs_322_258", "[#6]@-[#7]");
		queries.put("maccs_322_259", "[#6]@-[#8]");		
		queries.put("maccs_322_260", "[#6]@-[#16]");
		queries.put("maccs_322_261", "[#6]@-[#17]");		
		queries.put("maccs_322_262", "[#6]@-[#15]");		
		queries.put("maccs_322_263", "[#6]@-[#9]");
		queries.put("maccs_322_264", "[#6]@-[#35]");		
		queries.put("maccs_322_265", "[#6]@-[#14]");		
		queries.put("maccs_322_266", "[#6]@-[#53]");
		queries.put("maccs_322_267", "[#6]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_268", "[#7]@-[#7]");		
		queries.put("maccs_322_269", "[#7]@-[#8]");
		queries.put("maccs_322_270", "[#7]@-[#16]");
		queries.put("maccs_322_271", "[#7]@-[#17]");
		queries.put("maccs_322_272", "[#7]@-[#15]");		
		queries.put("maccs_322_273", "[#7]@-[#9]");		
		queries.put("maccs_322_274", "[#7]@-[#35]");
		queries.put("maccs_322_275", "[#7]@-[#14]");		
		queries.put("maccs_322_276", "[#7]@-[#53]");
		queries.put("maccs_322_277", "[#7]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_278", "[#8]@-[#8]");		
		queries.put("maccs_322_279", "[#8]@-[#16]");
		queries.put("maccs_322_280", "[#8]@-[#17]");		
		queries.put("maccs_322_281", "[#8]@-[#15]");
		queries.put("maccs_322_282", "[#8]@-[#9]");		
		queries.put("maccs_322_283", "[#8]@-[#35]");		
		queries.put("maccs_322_284", "[#8]@-[#14]");
		queries.put("maccs_322_285", "[#8]@-[#53]");		
		queries.put("maccs_322_286", "[#8]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_287", "[#16]@-[#16]");
		queries.put("maccs_322_288", "[#16]@-[#17]");		
		queries.put("maccs_322_289", "[#16]@-[#15]");
		queries.put("maccs_322_290", "[#16]@-[#9]");
		queries.put("maccs_322_291", "[#16]@-[#35]");
		queries.put("maccs_322_292", "[#16]@-[#14]");		
		queries.put("maccs_322_293", "[#16]@-[#53]");		
		queries.put("maccs_322_294", "[#16]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_295", "[#17]@-[#17]");		
		queries.put("maccs_322_296", "[#17]@-[#15]");		
		queries.put("maccs_322_297", "[#17]@-[#9]");
		queries.put("maccs_322_298", "[#17]@-[#35]");		
		queries.put("maccs_322_299", "[#17]@-[#14]");
		queries.put("maccs_322_300", "[#17]@-[#53]");		
		queries.put("maccs_322_301", "[#17]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_302", "[#15]@-[#15]");
		queries.put("maccs_322_303", "[#15]@-[#9]");		
		queries.put("maccs_322_304", "[#15]@-[#35]");
		queries.put("maccs_322_305", "[#15]@-[#14]");
		queries.put("maccs_322_306", "[#15]@-[#53]");		
		queries.put("maccs_322_307", "[#15]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_308", "[#9]@-[#9]");		
		queries.put("maccs_322_309", "[#9]@-[#35]");
		queries.put("maccs_322_310", "[#9]@-[#14]");
		queries.put("maccs_322_311", "[#9]@-[#53]");		
		queries.put("maccs_322_312", "[#9]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_313", "[#35]@-[#35]");		
		queries.put("maccs_322_314", "[#35]@-[#14]");
		queries.put("maccs_322_315", "[#35]@-[#9]");
		queries.put("maccs_322_316", "[#35]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_317", "[#14]@-[#14]");		
		queries.put("maccs_322_318", "[#14]@-[#53]");
		queries.put("maccs_322_319", "[#14]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_320", "[#53]@-[#53]");
		queries.put("maccs_322_321", "[#53]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_322", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		
		return queries;
	}


	public static LinkedHashMap<String, String> getRINAtomFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic", "[*;a]");
		queries.put("sulfuric acid",
				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide", "[F,Cl,Br,I]");
		queries.put("acyl chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put(
				"organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[c;R1]1[c;R1][c;R1][o;R1][c;R1]1");
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put(
				"hydroxylamine", 
				"[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
//				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		// queries.put("hydroxylamine","[NX3H2,$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][OX2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
		queries.put("phenol",
				"[#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("primary carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("tertiary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl halide", "[F,Cl,Br,I][c]");
		queries.put("C ONS bond", "[#6]~[#7,#8,#16]");
		queries.put(
				"13-dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary aliph ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary arom ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]-1=[O;X1]");
		queries.put("p-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1](=O)-[#6;R1]=[#6;R1]-1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																											// this
																											// class
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																						// this
																						// class
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put(
				"thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put(
				"amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put(
				"n-hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7;X3](-[#6,#1])-[#7;X3](-[#6,#1])-[#6,#1]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1-(acyloxy)ethyl ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put(
				"(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine_monocyclic", "[#6;R1]1=,:[#6;R1][#7;R1]=,:[#6;R1][#7;R1]=,:[#6;R1]1");
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[$([#6]-[#7;X2]=[C;X2]=[S;X1]),$([SX2][CX2]#[NX1])]");
		queries.put("alpha_beta-unsaturated_system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
		queries.put(
				"CH-acidic",
				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put(
				"CH-acidic strong",
				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N-aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary aliph amine",
				"[NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary aliph amine",
				"[NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary aliph amine",
				"[NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("primary arom amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary arom amine",
				"[NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary arom amine",
				"[NX3H0+0,NX4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("peroxo group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put(
				"thioamide S-oxide derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put(
				"thioamide SS-dioxide derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n-nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c-nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic acid ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put(
				"sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put(
				"steroid",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene_monocyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#16;R1][#6;R1]=,:1");
		queries.put("thiophene", "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1");		
		queries.put("thiophene s-oxide", "[#8;X1-][S+]-,:1-,:[#6]=,:[#6][#6]=,:[#6]-,:1");
		queries.put("1,3-thiazole", "[#6]1=,:[#6][#16][#6]=,:[#7]1"); // SMARTCyp - a 2D-method for Prediction of Cytochrome P450 Mediated Drug Metabolism_Supplement
		queries.put("1,2-thiazole", "[#6]=,:1[#6]=,:[#7][#16][#6]=,:1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "[#6]1=,:[#6][#6]=,:[#7][#7]=,:[#6]1");
		queries.put("pyridazine", "c1ccnnc1");
		queries.put("13-dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put(
				"biphenyl",
				"[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put(
				"12-aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur compound", "[#6]~[#16]");
		queries.put("secondary aliphatic/aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put(
				"propargyl-type 13-dipolar organic compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put(
				"diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put(
				"phenol ether",
				"[$([#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");
		queries.put("heteroaromatic", "[a;!c]"); // 107 bits
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put(
				"hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put(
				"hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put(
				"steroid",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("phosphate monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3-acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n-acyl-aromatic alpha-amino-acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put(
				"n-acyl-aliphatic alpha-amino-acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[$([#6]-[#7;X3](-[#1,#6])-[#7;X3+](-[#8;X1-])=[O;X1]),$([#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3])]");
		queries.put("12-oxazole", "[#6]=,:1[#6]=,:[#7][#8][#6]=,:1");
		queries.put("13-oxazole", "[#6]1=,:[#6][#8][#6]=,:[#7]1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put(
				"coenzyme a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put(
				"2-N-linked ribose deriv.",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic aocohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("glycerol-3-phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");

//		Newly added ones
		queries.put("biguanide", 
				"[$([#7]!@-[#6](=[#7])!@-[#7]!@-[#6](!@-[#7])=[#7]),$([#7]!@-[#6](!@-[#7])=[#7]!@-[#6](!@-[#7])=[#7])]");
		queries.put("glycerolipid", 
				"[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]");
		queries.put("glycerophospholipid", 
				"[#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("glycerophosphonolipid", 
				"[#6;A]P([#8;A;X2H1,X1-])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]");
		queries.put("sphingolipid", 
				"[$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8]),$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X4]=[#6;A;H1X4]-[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8])]");
		queries.put("glycero-3-dithiophosphocholine", 
				"[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#8;X2]!@-P(!@-[#16])(!@=[S;X1])!@-[#8]!@-[#6;A;H2X4]!@-[#6;A;H1X4](!@-[#8])!@-[#6;A;H2X4]!@-[#8]");
		queries.put("glycero-3-thiophosphocholine", 
				"[#6;A;H3X4][N;R0+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#8]P([#8;A;X2H1,X1-])(=[S;X1])[#8;R0][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("bisphosphatidic acid", 
				"[#6]!@-[#6]!@-[#6]!@-[#6](!@=O)!@-[#8]!@-[#6]!@-[#6](!@-[#6]!@-[#8]!@-P(!@-[#8])(!@=O)!@-[#8]!@-[#6]!@-[#6](!@-[#6]!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6])!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6])!@-[#8]!@-[#6](!@=O)!@-[#6]!@-[#6]!@-[#6]");
		queries.put("lysobisphosphatidic acid", 
				"[$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0]),$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#8]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0]),$([#6;R0]-[#6;R0]-[#6;R0]-[#6;R0](=O)-[#8]-[#6;R0]-[#6;R0](-[$([OH]),$([O-])])-[#6;R0]-[#8]P([$([OH]),$([O-])])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[$([OH]),$([O-])])-[#8]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])]");
		queries.put("phosphatidylethanol", 
				"[#6;H3X4R0]-[#6;H2X4R0]-[#8;R0]P([!#1!#6;$([OX2H1]),$([O-])])(=O)[#8;R0]-[#6;H2X4R0]-[#6;H1X4R0](-[#6;H2X4R0]-[#8;R0]-[#6;R0](-[#6,#1])=[O;R0])-[#8;R0]-[#6;R0](-[#6,#1])=[O;R0]");		
//		queries.put("phosphoglycosphingolipid", 
//				"");
		queries.put("semilysobisphosphatidic acid", 
				"[H][#8;R0]-[#6;R0](-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])-[#6;R0]-[#8]P([#8])(=O)[#8]-[#6;R0]-[#6;R0](-[#6;R0]-[#8;R0]-[#6;R0](=O)-[#6;R0]-[#6;R0]-[#6;R0])-[#8]-[#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]");
		queries.put("3-benzazepine", 
				"N1C=CC2=CC=CC=C2C=C1");
		queries.put("dibenzo_p_dioxin", 
				"[#8]-1-c2ccccc2-[#8]-c2ccccc-12");		
		queries.put("13-oxazin-2-one",
				"[$(O=C1OC=CC=N1),$(O=C1NC=CCO1)]"); // dx.doi.org/10.1021/ml500297n | ACS Med. Chem. Lett. 2014, 5, 1156−1161
		queries.put("13-oxazine", 
				"[$(C1NC=C=CO1),$([#6]-1-[#8]-[#6]=[#7]-[#6]=[#6]-1),$([#6]-1-[#6]=[#6]-[#8]-[#6]=[#7]-1),$([#6]-1-[#8]-[#6]=[#6]-[#6]=[#7]-1),$([#6]-1-[#8]-[#6]-[#7]=C=[#6]-1)]"); // dx.doi.org/10.1021/ml500297n | ACS Med. Chem. Lett. 2014, 5, 1156−1161
		queries.put("benzo-14-diazepine", 
				"[$([#6]-,=1-,=[#6]-[#7]=[#6]-c2ccccc2-[#7]-,=1),$([#6]-,=1-,=[#6]-[#7]=[#6]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#7]-,=1)]");		
		queries.put("1-benzofuran", 
				"[$(c1cc2ccccc2o1),$(O1C=CC2=CC=CC=C12)]");
		queries.put("benzimidazole", 
				"[$(c1nc2ccccc2n1),$(N1C=NC2=CC=CC=C12)]");		
		queries.put("sulfonylurea", 
				"[$([#7;X3]!@-[#6](=O)!@-[#7]=S(=O)=O),$([#6]S(=O)(=O)[#7]-[#6](-[#7;X3])=O),$([#6]S(=O)(=O)-[#7]=[#6](/[#7;X3])-[#8]),$([#6]S(=O)(=O)[#7]-[#6](-[#8])=[#7;X2])]");	
		
		queries.put("phosphate diester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2]-[#6]");
		queries.put("alpha-amino acid or derivatives", 
				"[#7;A][#6;X4]-[#6;X3]([!#1!#6])=[O;X1]");  //modified from ClasyFire's alpha-amino-acid-erivative-1 by adding X4 on the C2 carbon.
		queries.put("alpha-amino acid", 
				"[#7;A;X3,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X4]-[#6;X3]([#8;A;X2H1,X1-])=[O;X1]"); //modified from ClasyFire's alpha-amino-acid-1 by adding X4 on the C2 carbon.
		queries.put("tryptamine", 
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6;R1]-1=[#6;R1]-[#7;R1]-[#6;R2]-2=[#6;R2]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("flavonol", 
				"[H][#8]-[c;R1]1[#6;R1](=O)c2-,:cc-,:cc-,:c2[#8;R1][c;R1]1-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1]([#8;A;H1X2])[c;R1]-,:1");
		queries.put("flavone", 
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:1)-[#6;R1]=,:1[#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1](=O)[#6;R1]=,:1-[!#8]"); // PMID: 15914008
		queries.put("coumarin", 
				"[$(O=[#6]-1-[#8]-c2ccccc2-[#6]=[#6]-1),$(O=[#6]-1-[#8]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]=[#6]-1)]");
		queries.put("alkoxy group", 
				"[#6;A;X4][#8;X2][#6;A;H2X4][#7;X3]");
		queries.put("n_isopropyl group", 
				"[#6;A;H3X4][#6;H2X4]([#6;A;H3X4])-[#7;X3]");
		queries.put("thiono group", 
				"[*]=[SX1]");
		queries.put("s_2_hydroxycarboxylate", 
				"[#6]-[#6@H](-[#8;H1X2])-[#6;R0](-[#8;X1-])=[O;X1]"); 
		queries.put("hexoside", 
				"[#6]-[#8;H0X2]-[#6;R1]-1-[#8]-[#6;R1](-[#6])-[#6;R1](-[#8])-[#6;R1](-[#8])-[#6;R1]-1-[#8]");
		queries.put("phenothiazine", 
				"[$([#7]-1-c2ccccc2-[#16]-c2ccccc-12),$(N1C2=C(SC3=C1C=CC=C3)C=CC=C2)]"); // doi:10.1016/j.pharep.2015.04.005
		queries.put("phenoxazine", 
				"[$([#7]-1-c2ccccc2-[#8]-c2ccccc-12),$(N1C2=C(OC3=C1C=CC=C3)C=CC=C2)]");			
		queries.put("pyrrole", 
				"[$(c1ccnc1),$(N1C=CC=C1)]");			
		queries.put("pyrroline", 
				"C1CC=NC1");			
		queries.put("thioxanthene", 
				"[$([#6]-1-c2ccccc2-[#16]-c2ccccc-12),$(C1C2=CC=CC=C2SC2=CC=CC=C12)]");		
		queries.put("phthalazine", 
				"[#6]-1[#7]-[#7][#6]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-1-2");			
		queries.put("isoflavan", 
				"[H][C;R1]-,:1([H])-,:[#6;R1]=,:[#6;R1]([#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]-,:1=,:2)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");	

		// Added April 26th, 2016
		queries.put("pyridine", 
				"[$([#6;R1]-1=[#6;R1]-[#6;R1]=[#7;R1]-[#6;R1]=[#6;R1]-1),$([c;R1]1[c;R1][c;R1][n;R1][c;R1][c;R1]1)]");		
		
		queries.put("xanthone", 
				"O=[#6]1[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");		
		
		queries.put("fluorene", 
				"[#6]1[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");			
		
		queries.put("naphthoquinone", 
				"O=[#6]1[#6]=,:[#6][#6](=O)[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]1=,:2");
		
		queries.put("urea", 
				"[#7;X3;!$([#7][!#6])][#6;X3]([#7;X3;!$([#7][!#6])])=[O;X1]");
		
		queries.put("thiourea", 
				"[#7;X3;!$([#7][!#6])][#6;X3]([#7;X3;!$([#7][!#6])])=[S;X1]");
		
		queries.put("guanidine", 
				"[#7;A;v3X3,v4X4+][#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		
		queries.put("nitrosothiol", 
				"[#6]-[#16;X2]-[#7;X2]=[O;X1]");
		
		queries.put("sugar_pattern_1", 
				"[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
		
		queries.put("sugar_pattern_2", 
				"[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");		
		
		queries.put("sugar_pattern_combi", 
				"[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");		
		
		queries.put("sugar_pattern_2_reducing", 
				"[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
		
		queries.put("sugar_pattern_2_alpha", 
				"[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");		
		
		queries.put("sugar_pattern_2_beta", 
				"[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
		
		queries.put("thiolactam", 
				"[#6;R][#6;X3R]([#7;X3;$([H1][#6;!$(C=[O,N,S])]),$([H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[S;X1]");
		
		queries.put("carbothiolactone", 
				"[#6][#6X3R](=[SX1])[#16X2][#6;!$(C=[O,N,S])]");
		
		queries.put("imidothiolactone", 
				"[#6R][#6X3R](=,:[#7X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])[SX2][#6;!$(C=[O,N,S])]");
		
		queries.put("barbiturate", 
				"[$([#8]-,=[#6]-1-,=[#7]-[#6](-,=[#8])-[#6]-[#6](-,=[#8])-[#7]-1),$([#8]-,=[#6]-1-[#6]-[#6](-,=[#8])-[#7]=[#6](-[#8])-[#7]-1),$([#8]-,=[#6]-1-[#6]-[#6](-,=[#8])-[#7]-[#6](=O)-[#7]-1),$([$([O-]C1=NC(=O)NC(=O)C1),$([O-]C1=NC([O-])=NC(=O)C1),$([O-]C1=NC(=O)N=C([O-])C1)])]");		

		queries.put("benzodioxole", "[$(C1Oc2ccccc2O1),$(C1OOc2ccccc12)]");
		
		queries.put("glucosinolate","[#6][#6](-[#16]-[#6;R1]-1-[#8]-[#6](-[#6]-[#8])-[#6;R1]([!#1!#6;O,$([O-])])-[#6;R1]([!#1!#6;O,$([O-])])-[#6;R1]-1[!#1!#6;O,$([O-])])=[#7]/[#8]S([!#1!#6;O,$([O-])])(=O)=O"); // Slightly modified from the ClassyFire SMARTS
		
		queries.put("psoralen", "O=C1Oc2cc3occc3cc2C=C1"); // Guo LQ et al / Acta Pharmacol Sin 2004 Feb; 25 (2): 129-136
	
		queries.put("carbonate_salt", "[$([+1,+2,+3,+4,+5,+6,+7].[#8;A;X2H1,X1-]!@-[#6](!@-[#8;X1-])=O),$([#8;A;X2H1,X1-]!@-[#6](=O)!@-[#8;X2]-[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86])]");
		
		queries.put("carbonic_acid_derivative","[O,X]-[#6](-[O,X])=O");		
		
		queries.put("dibenzocyclooctadiene","[#6;R1]-1-[#6;R1]-[#6;R1]-c2cccc[c;R2]2-[c;R2]2ccccc2-[#6;R1]-1"); // 10.1124/dmd.104.000646 DMD December 2004 vol. 32 no. 12 1351-1358
		
		queries.put("benzenesulfonamide","[#7]S(=O)(=O)!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1"); // 10.1124/dmd.104.000646 DMD December 2004 vol. 32 no. 12 1351-1358
		
		queries.put("3_phenylpyrazole","c1cc(nn1)-!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		
		queries.put("hetero_n_basic_h","[nX3H1+0]"); // doi: 10.3389/fphar.2015.00123 Modeling of interactions between xenobiotics and cytochrome P450 (CYP) enzymes
		
		queries.put("hetero_n_basic_no_h","[nX3H0+0]");	// doi: 10.3389/fphar.2015.00123 Modeling of interactions between xenobiotics and cytochrome P450 (CYP) enzymes

		queries.put("hetero_non_basic","[nX2,nX3+]");
		
		queries.put("piperidine_monocyclic","[#6;R1]-1-[#6;R1]-[#6;R1]-[#7;R1]-[#6;R1]-[#6;R1]-1");
		queries.put("piperidine","C1CCNCC1");
		queries.put("4-aminopiperidine_monocyclic","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6;R1]-1-[#6;R1]-[#6;R1]-[#7;R1]-[#6;R1]-[#6;R1]-1");
		queries.put("4-aminopiperidine","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6]-1-[#6]-[#6]-[#7]-[#6]-[#6]-1");
		queries.put("halobenzene","[F,Cl,Br,I]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("benzene","[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");

		// added on May 8, 2016	
		
		queries.put("14_dihydropyridine","C1C=CNC=C1");
		queries.put("dihydropyrimidinethione","[O;X1R0]=[#6;R1]-1-[#7]-[#6]-[#6]=[#6]-[#7]-1");
		queries.put("quinazoline","[#6]=,:1[#6]=,:[#6][#6]=,:2[#7][#6]=,:[#7][#6][#6]=,:2[#6]=,:1");
		queries.put("2_phenylpyrimidine_unfused_benzene","[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)!@-[#6]-1=[#7]-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("3_phenylpyrimidine_unfused_benzene","[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6]-1=[#7]-[#6]=[#7]-[#6]=[#6]-1");

		queries.put("13_thiazolidine","C1CSCN1");
		queries.put("carboxylic-acid-imide-N-substituted","[#6]-[#7;X3](-[#6;X3](-[#6,#1])=[O;X1])-[#6;X3](-[#6,#1])=[O;X1]");
		queries.put("carboxylic-acid-imide-N-unsubstituted","[H][#7;X3](-[#6;X3](-[#6,#1])=[O;X1])-[#6;X3](-[#6,#1])=[O;X1]");
		queries.put("vinylogous-halide","[#9,#17,#35,#53;A;X1]-,:[#6;X3]=,:[#6;X3][#6;X3]=[O;X1]"); // http://ccj.springeropen.com/articles/10.1186/1752-153X-7-49
		queries.put("vinylogous-amide","[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;X3]=,:[#6;X3][#6;X3]=[O;X1]");
		queries.put("tetralin","[#6]1=,:[#6][#6]=,:[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]2=,:[#6]1");
		queries.put("n_phenylurea","[#7]!@-[#6](=O)!@-[#7]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("piperazine","C1CNCCN1");
		queries.put("azacycle","[#6]@[#7]");
		queries.put("oxacycle","[#6]@[#8]");
		queries.put("organoheterocycle","[!C!c!R0]");
		queries.put("filter_bridged_rings_66_heterocyclic","*~1~*~*~*(~*~*~1)!@*~1~*~*~*~*~*~1");
		queries.put("filter_bridged_rings_65_heterocyclic","*~1~*~*~*(~*~1)!@*~1~*~*~*~*~*~1");
		queries.put("filter_bridged_rings_55_heterocyclic","*~1~*~*~*(~*~1)!@*~1~*~*~*~*~1");

		queries.put("filter_fused_rings_67_heterocyclic","*~1~*~*~*~2~*~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_66_heterocyclic","*~1~*~*~*~2~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_65_heterocyclic","*~1~*~*~2~*~*~*~*~*~2~*~1");
		queries.put("filter_fused_rings_55_heterocyclic","*~1~*~*~2~*~*~*~*~2~*~1");
		
		
		
		queries.put("aliphatic_chain_13","[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]");
		queries.put("aliphatic_chain_4","[*;A;R0]~[*;A;R0]~[*;A;R0]~[*;A;R0]");
		queries.put("ortho_substituted_ring","*-!:aa-!:*");
		queries.put("meta_substituted_ring","*-!:aaa-!:*");
		queries.put("para_substituted_ring","*-!:aaaa-!:*");
		queries.put("aryl_fluoride","F[$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("1-methyl-23-dihydro-1H-imidazole-2-thione_monocyclic","[#6;A;H3X4][#7;A;R1]1[#6;A;R1]=[#6;A;R1][#7;A;R1][#6;A;R1]1=[S;X1]");
		queries.put("23-dihydro-1H-imidazole-2-thione_monocyclic","[#6;A;R1]1[#7;A;R1][#6;A;R1]=[#6;A;R1][#7;A;R1]1");
		queries.put("terminal_alkyne","[H]C!@#C*");
		queries.put("w_1_acetylene","[#6;A;H3X4]C!@#C*");			
		queries.put("dithiocarboxylic_acid_amide","[#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])])][#16;X2]-[#6;X3]([#6,#1;A])=[S;v2]");
		queries.put("thiocarboxylic_acid_amide","[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put("terminal_cc_bond","[#6]-[#6]=[#6;A;H2X3]");	
		queries.put("thiophene_sulfoxide","O=S1[#6]=[#6]-[#6]=[#6]1");
		queries.put("trifluoroalkyl_group","[#6;A;X4][C;X4](F)(F)F");
		queries.put("12_dithiole_3_thione","S=C1SSC=C1");	
		queries.put("piperidine_monocyclic","[#6;A;R1]1[#6;A;R1][#6;A;R1][#7;A;R1][#6;A;R1][#6;A;R1]1");
		queries.put("piperidine","C1CCNCC1");
		queries.put("sulfenyl_group","[#6]-[#16;v2X2R0]-*");
		queries.put("sulfenyl_halide","[#6]-[#16;v2X2R0]-[F,Cl,Br,I]");	
		queries.put("thiophosphate_diester","[#6]-[#8;X2][P;X4]([#8;A;H1X2])(=[S;X1])[#8;X2]-[#6]");
		queries.put("45_dihydro_123_oxadiazol_5_one","O=C1CN=NO1");
		queries.put("alkyl_sulfanyl_carbothioyl_amine","[#7;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;X3R0](-[#16;X2])=[S;X1R0]");	
		queries.put("123_triazolidine","[#6]-1-[#6]-[#7]=[#7]-[#7]-1");
		queries.put("1_amino_123_triazolidine","NN1CCN=N1");
		queries.put("n_cyclobutylbenzylamine","[#6;R0](-[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("n-cyclopropylbenzylamine","[#6;R0](-[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("cyclopropylamine","[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-1");
		queries.put("cyclobutylamine","[#7;R0;$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])])]-[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-1");	
		queries.put("123_thiadiazole","c1csnn1");	
		queries.put("124_thiadiazole","c1ncsn1");	
		queries.put("125_thiadiazole","c1cnsn1");	
		queries.put("134_thiadiazole","c1nncs1");	
		queries.put("xanthate","[$([!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86;+].[#6;A;X4][#8;X2]-[#6;X3](!@-[#16-])=[SX1]),$([#6;A;X4][#8;X2]-[#6;X3](=[SX1])!@-[#16;X2]-[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86])]");	
		queries.put("pyrazine","c1cnccn1");	 // pp440,445 in Robert M. Rydzewski. Real World Drug Discovery: A Chemist's Guide to Biotech and Pharmaceutical Research
		queries.put("11_disubstituted_hydrazine","[#6]-[#7;X3](-[#6])!@-[#7;A;H2X3]");	// Inhibition of Cytochrome P450 Enzymes p.266
		queries.put("monosubstituted_hydrazine","[#6]-[#7;H1X3]!@-[#7;A;H2X3]");	// Inhibition of Cytochrome P450 Enzymes p.266
		queries.put("carboxylic_acid_hydrazide","[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)]-[#7](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#7](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#6](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])=O");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("12-dihydroquinoline","[#6;A;X4]1[#7]-[#6]-2=[#6](-[#6]=[#6]1)-[#6]=[#6]-[#6]=[#6]-2");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano 
		queries.put("alkylhydrazine","[#6;A;X4][#7;X3](-[#6,#1])-[#7;X3](-[#6,#1])-[#6,#1]");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("hydrazone","[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#7;X2]=[#6]");	// Inhibition of Cytochrome P450 Enzymes; Maria Almira Correia and Paul R. Ortiz de Monteflano
		queries.put("lactone","[#6;!$(C=[O,N,S])][#8;X2][#6;X3R]([#6])=[O;X1]");	
		queries.put("lactam","[#6;R][#6;X3R]([#7;X3;$([H1][#6;!$(C=[O,N,S])]),$([H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[O;X1]");	
		queries.put("linear_furanocoumarin_pattern_1","O=[#6]-1-[#8]-[#6]-2=[#6](-[#6]=[#6]-1)-[#6;R1]=[#6]-1-[#6]=[#6;R1]-[#8]-[#6]-1=[#6;R1]-2");	// J Pharm Pharm Sci. 2001 Sep-Dec;4(3):217-27.Inhibition of human CYP3A4 activity by grapefruit flavonoids, furanocoumarins and related compounds.
		queries.put("stilbene","[$(*=[#6]-1-[#6]=[#6]-[#6]=[#6]-[#6;R1]-1!@-[#6]!@=[#6]~[c;R1]1ccccc1),$([#6](!@=[#6]~/[c;R1]1ccccc1)!@-[#6;R1]-1=[#6]-[#6]=[#6]-[#6]=[#6]1)]");	// Mol Nutr Food Res. 2008 Jun;52 Suppl 1:S77-83. doi: 10.1002/mnfr.200700202.
		queries.put("anthocyanidin","[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#8;R1+]=,:1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("o_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6]-1-[#6;R0]=[#6;R0]-[#6;R0](-*)=O");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("m_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6;R1]-1=[#6]-[#6](-[#6;R0]=[#6;R0]-[#6;R0](-*)=O)=[#6;R1]-[#6;R1]=[#6;R1]-1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("p_hydroxycinnamic_acid_or_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6](-[#6;R0]=[#6;R0]-[#6;R0](-*)=O)-[#6]=[#6;R1]-1");	// Loai Basheer and Zohar Kerem; Interactions between CYP3A4 and Dietary Polyphenols; http://dx.doi.org/10.1155/2015/854015
		queries.put("o_hydroxybenzoic_acid_der","[#8;A;X2R0][#6]-1=[#6](-[#6]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("m_hydroxybenzoic_acid_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6]-[#6](=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("p_hydroxybenzoic_acid_der","[#8;A;H1X2R0][#6;R1]-1=[#6;R1]-[#6;R1]=[#6](-[#6]=[#6;R1]-1)-[#6;R0](!@-*)=[O;X2]");	
		queries.put("1_benzopyran","[#6]-1-,=[#6]-c2ccccc2-[#8]-[#6]-1");	
		queries.put("2_benzopyran","[#6]-1-,=[#6]-c2ccccc2-[#6]-[#8]-1");	
		queries.put("chalcone_or_der","[#8;R0]-,=[#6;R0](-,=[#6;R0]-,=[#6;R0]!@-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1)!@-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1");	
		queries.put("phenoxy_compound","*-[#8;X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("thienopyridine","[$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#6]=[#7]-[#6]=[#6]-2),$([#6]~1~[#7]~[#6]-[#6]-2=[#6](-[#6]~1)-[#16]-[#6]=[#6]-2),$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#6]=[#6]-[#7]=[#6]-2),$([#16]-1-[#6]=[#6]-[#6]-2=[#6]-1-[#7]=[#6]-[#6]=[#6]-2)]");  //http://www.sciencedirect.com/science/article/pii/S0968089614007147
		queries.put("1_benzothiophene","[#6]1=,:[#6][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#16]1"); //http://www.sciencedirect.com/science/article/pii/S0968089614007147
		queries.put("2_benzothiophene","[#6]=,:1[#16][#6]=,:[#6]2[#6]=,:[#6][#6]=,:[#6][#6]=,:12"); //http://www.sciencedirect.com/science/article/pii/S0968089614007147
//		queries.put("kavalactone","");	 // Robert S. FOTI and Jan L.WAHLSTROM; The role of dietary supplements in cytochrome P450-mediated drug interactions 
//		queries.put("ginkgolide_or_bilabolide","");		 // Robert S. FOTI and Jan L.WAHLSTROM; The role of dietary supplements in cytochrome P450-mediated drug interactions 
		queries.put("cyclic_carboxylic_ester","[#6;R;!$(C=[O,N,S])]-[#8;X2R][#6;A;X3R;$([R0][#6]),$([H1R0])]=[O;X1R0]");	// for macrolides
//		queries.put("macrocycle","[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11]");	// for macrolides	
		queries.put("macrocycle_r12_r40","[r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r40]");	// for macrolides
		queries.put("chloroethylene_deriv","[Cl]-[#6](-[#6,#1])=[#6](/Cl)[#17,#1;A]");
		queries.put("conjug_bond","[*,#1]-[#6]=[#6]-[#6]=[#6]-[*,#1]");
		
		queries.put("curcumin","[$([#8]~[#6;A](~[#6;A;R0]~[#6;A]~[#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A]~[#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1),$([#8]~[#6;A](~[#6;A]~[#6;A;R0][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)~[#6;A]~[#6;A](~[#8])~[#6;A;R0]~[#6;A][#6]-1=[#6]-[#6](-[#8]-[#6])=[#6](-[#8])-[#6]=[#6]-1)]");	
		queries.put("unfused_benzene","[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");	
		queries.put("s_in_aromatic_5_ring_with_lone_pair","[sX2r5]");	
		queries.put("o_in_aromatic_5_ring_with_lone_pair","[oX2r5]");			
		queries.put("n_in_aromatic_5_ring","[nX2r5]");	
		//	 1 However, no N-oxides of nitrogen atoms in five-membered aromatic rings have been found. Instead, such nitrogen atoms tend to bind directly to the heme iron atom, leading to CYP inhibition. http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=3412406&fileOId=3412409
	
		queries.put("s_in_aromatic_6_ring_with_lone_pair","[sX2r6]");	
		queries.put("o_in_aromatic_6_ring_with_lone_pair","[oX2r6]");			
		queries.put("n_in_aromatic_6_ring","[nX2r6]");
		queries.put("pattern_c1", 
				"[#6][#8,#16;A][#6;H2X4R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // Tioconazole Miconazole sulconazole Miconazole
		queries.put("pattern_c2","[#6;A]-,=1-,=[#6;A][#6;A;R2]2[#6;A]-,=[#6;A]-,=[#6;A]-,=[#6;A][#6;A;R2]2[#6;A]-,=1");
		queries.put("pattern_c3","*-[#6;X3](=O)!@-[#6;X3](-[F,Cl,Br,I])-[F,Cl,Br,I]");	
		queries.put("pattern_c4","[#6]-[#16;X2]!@-[#6;X3](-[#6])=O");	// from BioTransformer fingerprint-based clustering
		queries.put("pattern_c5","[#6]-[#7;R0]=[#6;R0]-[#7;R0]-[#8;H1X2R0]");	
		queries.put("cyclic","[*;R]");		
		queries.put("134_oxadiazole","c1nnco1");	
		queries.put("124_oxadiazole","c1ncno1");		
		queries.put("cc_double_bond","[#6]=[#6]");	
		queries.put("ring_of_size3","*~1~*~*~1");
		queries.put("ring_of_size5","*~1~*~*~*~*~1");	
		queries.put("ring_of_size6","*~1~*~*~*~*~*~1");
		queries.put("ring_of_size7","*~1~*~*~*~*~*~*~1");
		queries.put("heteroring_of_size5","*~1~*~*~[!#1!#6]~*~1");	
		queries.put("heteroring_of_size6","*~1~*~*~[!#1!#6]~*~*~1");	
		queries.put("heteroring_of_size7","*~1~*~*~*~[!#1!#6]~*~*~1");	
		queries.put("oxolane","C1CCOC1");
		queries.put("imidazoline","C1CN=CN1");
		queries.put("benzo_13_dioxolane","[#6]-1-[#8]-c2ccccc2-[#8]-1");

		
		// J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280
		queries.put("pattern_c6","*-*(-*)-[*;X2]-*(-*)-*");		//atom has two neighbors, each with three or more neighbors (including the central atom). (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("pattern_c7","[!#1!#6]*[!#1!#6]"); // atom with at least two heteroatom neighbors (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("pattern_c8","*!@-*!@-*"); //* atom with more than one chain bond (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("n_halide_bond","[#7]-[F,Cl,Br,I]"); // N-X bond (J. Chem. Inf. Comput. Sci. 2002, 42, 1273-1280)
		queries.put("nn_double_bond","[#7]=[#7]");
		queries.put("cn_double_bond","[#6]=[#7]");
		queries.put("cs_double_bond","[#6]=[#16]");
		queries.put("cn_single_bond","[#6]-[#7]");
		queries.put("co_single_bond","[#6]-[#8]"); // C-O single bond
		
		
		// Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition (PMID: 26108525). J Chem Inf Model. 2015 Jul 27;55(7):1426-34. doi: 10.1021/acs.jcim.5b00130. Epub 2015 Jul 2.
		queries.put("r5_1ht_sat","[#6]-1-[#6]-[#6][!#1!#6][#6]-1"); // saturated 5-membered ring with 1 heteroatom. Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition 
		queries.put("r5_2ht_sat","[$([#6]-1-[#6][!#1!#6][!#1!#6][#6]-1),$([#6]1-[#6][!#1!#6][#6][!#1!#6]1)]"); // saturated 5-membered ring with 2 heteroatoms
		queries.put("r5_3ht_sat","[$([#6]1-[#6][!#1!#6][!#1!#6][!#1!#6]1),$([#6]1[!#1!#6][#6][!#1!#6][!#1!#6]1),$([#6]1[!#1!#6][#6][!#1!#6][!#1!#6]1)]");  // saturated 5-membered ring with 3 heteroatoms
		queries.put("r5_4ht_sat","[#6]1[!#1!#6][!#1!#6][!#1!#6][!#1!#6]1");  // saturated 5-membered ring with 4 heteroatoms
		queries.put("r5_1ht_aro","c1cc[!#1!#6]c1"); // aromatic 5-membered ring with 1 heteroatom. Bo-Han Su et al. Rule-based Prediction Models of Cytochrome P450 Inhibition 
		queries.put("r5_2ht_aro","[$(c1c[!#1!#6][!#1!#6]c1),$(c1c[!#1!#6]c[!#1!#6]1)]"); // aromatic 5-membered ring with 2 heteroatoms
		queries.put("r5_3ht_aro","[$(c1c[!#1!#6][!#1!#6][!#1!#6]1),$(c1[!#1!#6]c[!#1!#6][!#1!#6]1),$(c1[!#1!#6]c[!#1!#6][!#1!#6]1)]");  // aromatic 5-membered ring with 3 heteroatoms
		queries.put("r5_4ht_aro","c1[!#1!#6][!#1!#6][!#1!#6][!#1!#6]1");  // aromatic 5-membered ring with 4 heteroatoms
		queries.put("carbon_atom","[#6]");
		queries.put("oxygen_atom","[#8]");
		queries.put("nitrogen_atom","[#7]");
		queries.put("sulfur_atom","[#16]");
		queries.put("phos_atom","[#15]");		

		queries.put("het_atom","[!#1!#6]");
		queries.put("organometallic","[#6]~[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86]");		

		queries.put("aliph_carboxylic","[#6;X4]-[#6]([#8;A;X2H1,X1-])=O");
		queries.put("sat_carb_r6","C1CCCCC1"); //  saturated carbon-only ring size 6

		queries.put("1_sulfenylpiperidine","O=[S][#7]-1-[#6]-[#6]-[#6]-[#6]-[#6]-1");
		queries.put("pattern_c10","CCC(C)N"); // *
		queries.put("n_o_single_bond","[#7]-[#8]");
		queries.put("ketone","[#6]-[#6;X3](-[#6])=[O;X1]");

		// Generate using SOCO Chepelev et al.
		
		queries.put("pattern_c9","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;X2;!$([#7][!#6])])[#8;X2;!$([#8][!#6]),OX1-]");
//		queries.put("pattern_c9","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;X2;!$([#7][!#6])])[#8;X2;!$([#8][!#6]),OX1-]");
		queries.put("pattern_c53","[#6;R][#6;R][#6;R][#6;R][#6;R]");
		queries.put("pattern_c11","[#6][#6][#6]([#6])[#6]");
		queries.put("pattern_c12","[#6]~[#6](~[#6])[#6]");
		queries.put("pattern_c13","[#6]-[#9,#17,#35,#53,#85]");
		queries.put("pattern_c14","[#6X3](=[OX1])[#6X3]=,:[#6X3][#7,#8,#16,F,Cl,Br,I]");
		queries.put("pattern_c15","[$([*@@](~*)(~*)(*)*),$([*@@H](*)(*)*),$([*@@](~*)(*)*),$([*@@H](~*)~*)]");
		queries.put("pattern_c16","C=CCCCC");
		queries.put("pattern_c17","C=CC=CCC");
		queries.put("pattern_c18","[#6X4H3]CN");
		queries.put("pattern_c19","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("pattern_c20","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("pattern_c21","[!#6][#6X3](=[!#6])[!#6]");
		queries.put("pattern_c22","[#7X3;!$([#7][!#6])][#6X3](=[OX1])[#7X3;!$([#7][!#6])]");
		queries.put("pattern_c23","[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put("pattern_c24","[R]~[*;!R]~[*;!R]~[*;!R]~[R]");
		queries.put("pattern_c25","[R]~[*;!R]~[*;!R]~[*;!R]~[*;!R]~[R]");
		queries.put("pattern_c26","[#6][#7][#6X4H3]");
		queries.put("pattern_c27","[#6X4H3][#6][#7]");
		queries.put("pattern_c28","[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("pattern_c29","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]");
		queries.put("pattern_c30","[#6X3](=[OX1])[#6X3]=,:[#6X3][#6;!$(C=[O,N,S])]");
		queries.put("pattern_c31","[$([#7X2,OX1,SX1]=,:**=,:*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=**=*),$([#7X3,OX2,SX2;!H0]*=,:**:n)]");
		queries.put("pattern_c32","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*)]");
		queries.put("pattern_c33","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*),$([*@@](~*)(~*)(*)*),$([*@@H](*)(*)*),$([*@@](~*)(*)*),$([*@@H](~*)~*)]");
		queries.put("pattern_c34","OCC");
		queries.put("pattern_c35","N(C)C");
		queries.put("pattern_c36","[*](=O)[!#6][*;!H]");
		queries.put("pattern_c37","[$([#6X3H0][#6]),$([#6X3H])](=[!#6])[!#6]");
		queries.put("pattern_c38","[$([#7X2,OX1,SX1]=*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=*),$([#7X3,OX2,SX2;!H0]*:n)]");
		queries.put("pattern_c39","[#6]~[#6]~[#7]");
		queries.put("pattern_c40","[#6]#[#7]");			
		queries.put("pattern_c41","[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("pattern_c42","[FX1][CX4;!$([H0][Cl,Br,I]);!$([F][C]([F])([F])[F])]([FX1])([FX1])");	
		queries.put("pattern_c42","[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("pattern_c43","[FX1][CX4;!$([H0][Cl,Br,I]);!$([F][C]([F])([F])[F])]([FX1])([FX1])");
		queries.put("pattern_c44","OCCCC");
		queries.put("pattern_c45","[#6][#6][#8][#6]");
		queries.put("pattern_c46","[#6]~[#8]~[#6]");
		queries.put("pattern_c47","[#6]~[#6]~[#6]~[#8]~[#6]");
		queries.put("pattern_c48","[#6][#6][#8][#6][#6]");
		queries.put("pattern_c49","[#6;R][#6;R][#6;R][#8;R]");
		queries.put("pattern_c50","[#6;R][#6;R][#6;R][#6;R][#8;R]");
		
// these should also be added to the fingerprint for SOM prediction
//		queries.put("acetanilide","[#6;X4]-[#6;X3](=[O;X1])-[#7;X3]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("aniline_der","[#7;X3;H1,H2]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("o_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1-[#8;H1X2]");
		queries.put("p_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1](-[#8;H1X2])[c;R1][c;R1]1");
		queries.put("aryl_urea","[#6;a]-[#7;X3]-[#6;X3](-[#7;X3;!$([#7][!#6])])=[O;X1]");
		queries.put("aryl_carbamate","[#6;a;R1]!@-[#7]!@-[#6](=O)-[#8]-[#6;!$(C=[O,N,S])]");
		queries.put("hydroquinone","[H][#8]!@-[c;R1]1[c;R1](-[*,#1;!$([OH])])[c;R1](-[*,#1;!$([OH])])[c;R1](!@-[#8][H])[c;R1](-[*,#1;!$([OH])])[c;R1]1-[*,#1;!$([OH])]");
		queries.put("n_nitrosourea","[#7]!@-[#6;X3](!@=[O;X1])!@-[#7;X3]!@-[#7;X2]!@=O");
		queries.put("n_mustard","[$(Cl[#6]-[#6]-[#7]-[#6]-[#6]Cl),$(F[#6]-[#6]-[#7]-[#6]-[#6]F),$(Br[#6]-[#6]-[#7]-[#6]-[#6]Br),$(I[#6]-[#6]-[#7]-[#6]-[#6]I)]");
		
// Contributions of Human Enzymes in Carcinogen Metabolism : Chem Res Toxicol. 2012 July 16; 25(7): 1316–1383. doi:10.1021/tx300132k.
		queries.put("aza_aromatic","[#7;a;R]");
		queries.put("oxazaphosphorin_2_amine_2_oxide","[#7][P;R1]1(=[O;X1])[#7;X3R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#8;R1]1");
		queries.put("nitrosamine","[#7;!$(N*=O)]-[#7;X2]=[O;X1]");
		queries.put("beta_lactam","O=[#6]-1-[#6]-[#6]-[#7]-1");
//		queries.put("michael_acceptor","");	 

//		Structural Alerts of Mutagens and Carcinogens: Current Computer-Aided Drug Design, 2006, Vol. 2, No. 2
//		
//		Alkyl phosphonate
//		Alkyl sulfonates
//		Alkyl hydrazines
//		Alkyl aldehydes
//		N-methylol derivatives
//		2-chloroethyl mustards
//		N-chloroamines
//		propiolactones
//		propiosultones
//		aromatic aziridinyl derivatives
//		aliphatic aziridinyl derivatives
//		alkyl N-nitrosoamines
//		cyclic N-nitrosamines
		
		queries.put("pattern_c51","[#6;H3X4]-[#6]-[#7;X3](-[#6;H3X4])-[#7;X2]=[O;X1]");
		queries.put("pattern_c52","[#6;A;H3X4][#6]-[#7;X3]([#6;A;H3X4])-[#6][#6;A;H3X4]");	
		

		// Electron donating group / Electron withdrawing group :  http://www.chem.ucalgary.ca/courses/350/Carey5th/Ch12/ch12-8b.html
		// Add those
		// Read this link: http://www.chem.ucalgary.ca/courses/350/Carey5th/Ch12/ch12-8d.html

//		Strongly activating EDG
		queries.put("phenoxide_ortho_subs","[#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*");
		queries.put("phenoxide_para_subs","[c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1");
		queries.put("prim_amin_ortho_subs","[#7;A;H2X3][c;R1]:[*;R1]!@-*");
		queries.put("prim_amin_para_subs","[#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*");
		queries.put("ter_amin_ortho_subs","*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]");
		queries.put("ter_amin_para_subs","*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]");	
		queries.put("phenol_ortho_subs","[#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*");
		queries.put("phenol_para_subs","[#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1");
		queries.put("phenolether_ortho_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenolether_para_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1");		
		queries.put("phenolether_ortho_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenolether_para_subs","[c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1");		

//		 Moderately activating EDG		
		queries.put("anilide_n_ortho_subs","[#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1!@-*"); 		// Kalgutkar AS. et al (PMID: 15975040; Current Drug Metabolism, 2005, 6, 161-225)
		queries.put("anilide_n_para_subs","[#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1][#6;R1]=,:1");	// Kalgutkar AS. et al (PMID: 15975040; Current Drug Metabolism, 2005, 6, 161-225)	
		queries.put("phenolester_ortho_subs","[#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1!@-*");
		queries.put("phenolester_para_subs","[#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1][#6;R1]=,:1");	

//		 Weakly activating EDG	
		queries.put("phenylalkyl_ortho_subs","[#7;X3]!@-[#6]([#6,#1;A])!@=O.[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*");
		queries.put("phenylalkyl_para_subs","[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1");
		queries.put("biphenyl_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("biphenyl_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");	
		queries.put("phenylvinyl_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]");
		queries.put("phenylvinyl_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1");
		
		
//		Weakly deactivating EWG
		queries.put("phenylhalide_ortho_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]");		
		queries.put("phenylhalide_para_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1");
		
		// Modelately deactivating EWG	
		queries.put("benzaldehyde_meta_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1");
		queries.put("phenylketone_meta_subs","[#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");		
		queries.put("benzoate_meta_subs","[#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzoic_acid_meta_subs","[#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzoyl_chloride_meta_subs","[Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");		

		// Strongly deactivating EWG		
		
		queries.put("benzotrifluoride_meta_subs","[F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");
		queries.put("benzotrichloride_meta_subs","[Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");			
		queries.put("benzonitrile_meta_subs","*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]");	
		queries.put("benzenesulfate_meta_subs","[#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("anilinium_meta_subs","[#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("phenyl_1qam_meta_subs","[#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	// benzene C1-substituted with a quaternary ammoniuam salt and C3- substituted any atom.
		queries.put("nitrobenzene_meta_subs","[#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	

		
		
		queries.put("edg","[$([#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1),$([#7;A;H2X3][c;R1]:[*;R1]!@-*),$([#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*),$(*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$(*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#7;X3]!@-[#6]([#6,#1;A])!@=O.[#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1)]");	
		queries.put("ewg","[$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1),$([#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]),$([#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1)]");	
		
		
		
		// CYP isoform specificity toward drug metabolism: analysis using common feature hypothesis
		// J Mol Model (2012) 18:709–720: DOI 10.1007/s00894-011-1105-5
		
		queries.put("isopropyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-*");
		queries.put("isobutyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-[#6;A;X4;H2,H3]");
//		queries.put("edg","");
//		queries.put("ewg","");		
		// http://www.ifm.liu.se/compchem/msi/doc/life/catalyst46/tutorials/11_excludeOrTool.doc.html
		// prim. amine OR sec. amine OR ter. amine OR amidine OR guadidino OR amidineH OR (positively charged center except dipole)
		queries.put("posionazable","[$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]),$([#7;AH1;v3X3,v4X4+][#6;X3]([#7;AH1;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]),$([#7;A;H1X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6]);!$([C][N])]=[#7;A;X2;!$(NC=[O,S])]);!$([+1,+2,+3,+4,+5,+6,+7]~[+0,+1,+2,+3,+4,+5,+6,+7])]");		
//		queries.put("","");
//		queries.put("","");		
		
		
		queries.put("acyl_carnitine","[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6](!@-[#6;A;H2X4]!@-[#6](!@-[#8;A;X1-,X2H1])!@=O)!@-[#8]!@-[#6](!@-[#6,#1;A])!@=O");
		queries.put("acyl_choline","[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#8]!@-[#6]([#6,#1;A])!@=O");
		queries.put("thiohemiacetal","[SX2]([#6;!$(C=[O,S,N])])[CX4;!$(C(S)(S)[!#6])][OX2H]");		
		
		queries.put("thioacetal","[#6]-[#16;X2][C;X4]([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#16;X2]-[#6]");
		queries.put("hemiaminal","[NX3v3;!$(NC=[#7,#8,#15,#16])]([#6])[CX4;!$(C(N)(N)[!#6])][OX2H]");
		queries.put("aminal","[*,#1;CX4,#1]-[#7](-[*,#1;CX4,#1])C([*,#1;CX4,#1])([*,#1;CX4,#1])[#7](-[*,#1;CX4,#1])-[*,#1;CX4,#1]");
		queries.put("hemiacetal","[H][#8]C([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[#6]");		
		queries.put("acetal","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]C([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("imidothioester","[#6;!$(C=[O,N,S])]-[#16;X2][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");	
		queries.put("imidothioacid","[$([SX2H]),$([SX1-])][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");
		queries.put("Imidothioacid_cyclic","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])[$([SX2H]),$([SX1-])]");		
		queries.put("imidothiolactone","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])-[#16;X2]-[#6;!$(C=[O,N,S])]");
		queries.put("imidolactam","[#6][#6;X3R;$([H0](=[NX2;!$(N(=[#6X3][#7X3])C=[O,S])])[#7X3;!$(N([#6X3]=[#7X2])C=[O,S])]),$([H0](-[NX3;!$(N([#6X3]=[#7X2])C=[O,S])])=,:[#7X2;!$(N(=[#6X3][#7X3])C=[O,S])])]");	
		queries.put("imidoyl_halide","[#9,#17,#35,#53;A;X1][#6;A;X3R0;$([H0][#6]),$([H1])]=[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])]");
		queries.put("imidoyl_halide_cyclic","[#6;R][#6;X3R](=,:[#7;A;X2;$([H1]),$([H0][#6;!$(C=[O,N,S])])])-,:[#9,#17,#35,#53;A;X1]");
		queries.put("carbxoxylic_acid_amidrazone","[$([$([#6X3][#6]),$([#6X3H])](=[#7X2v3])[#7X3v3][#7X3v3]),$([$([#6X3][#6]),$([#6X3H])]([#7X3v3])=[#7X2v3][#7X3v3])]");	
		queries.put("alpha_hydroxyacid","[H][#8;v2]-[#6](-[*,#1;#1,C,$([cX3](:*):*),$([cX2+](:*):*)])-[#6](=O)-[#8][H]");
		queries.put("alpha_hydroxyaldehyde","[#6]-,=[#6;R0](-[#8;R0][H])-[#6;R0]([H])=[O;R0]");
		queries.put("carboxylic_acid_orthoester","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]C([#1,C,$([cX3](:*):*),$([cX2+](:*):*)])([#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)])[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("carbonic_acid_ester_halide","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](=[O;X1])-[#8;X2][#9,#17,#35,#53;A;X1]");
		queries.put("carbonic_acid_monoester","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](-[$([OX2H]),$([OX1-])])=[O;X1]");	
		queries.put("carbonic_acid_diester","[#6;!$(C=[O,N,S])][#8;X2][#6;X3](=[O;X1])[#8;X2][#6;!$(C=[O,N,S])]");
		queries.put("thiocarbonic_acid_derivative","[$([*,#1;!$(C=[O,N,S])][#8;X2][#6;X3](=O)[#16;X2][*,#1;!$(C=[O,N,S])]),$([*,#1;!$(C=[O,N,S])][#8;X2][#6;X3](=[O;X1])[#16;X2][*,#1;!$(C=[O,N,S])])]");	
		queries.put("thiocarbonic_acid_ester_halide","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](=[S;X1])-[#8;X2][#9,#17,#35,#53;A;X1]");
		queries.put("thiocarbonic_acid_monoester","[#6;!$(C=[O,N,S])][#8;A;X2;!R][#6;X3](-[$([OX2H]),$([OX1-])])=[S;X1]");	
		queries.put("thiocarbonic_acid_diester","[#6;!$(C=[O,N,S])][#8;X2][#6;X3](=[S;X1])[#8;X2][#6;!$(C=[O,N,S])]");
		queries.put("isourea","[#7;X3;!$([#7][!#6])][#6;X3](=,:[#7;A;X2;!$([#7][!#6])])[#8X2!$([#8][!#6]),OX1-]");	
		queries.put("aziridinone","[O;X1]=[#6;X3]-1-[#6]-[#7]-1");
			
		queries.put("heteroaromatic_s","[sX2]");
		queries.put("heteroaromatic_o","[o]");	

		queries.put("nitrosamide","O=[#6]-[#7]-[#7;X2]=[O;X1]");	
		queries.put("sulfuric_acid_monoester","[#6]-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		queries.put("sulfuric_acid_diester","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#8]S(=O)(=O)[#8]-[C,$([cX3](:*):*),$([cX2+](:*):*)]");	
		queries.put("sulfuric_acid_diamide","[#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][S;X4]([#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])(=[O;X1])=[O;X1]");
		queries.put("sulfuric_acid_amide_ester","[#15;A;X4;$([H3]=[OX1]),$([H2](=[OX1])[#6]),$([H1](=[OX1])([#6])[#6]),$([H0](=[OX1])([#6])([#6])[#6])]");
		queries.put("sulfuric_acid_derivative","[#7,#8;A]S([#7,#8;A])(=O)=O");	
		queries.put("sulfonic_acid","[C,$([cX3](:*):*),$([cX2+](:*):*)][S;v6]([!#1!#6;!S])(=O)=O");
		queries.put("hydroxylamine_o_sulfonic_acid","[#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		queries.put("sulfonic_acid_ester","[#6;!$(C=[O,N,S])]-[#8][S;v6]([#6,#1;A])(=[O;X1])=[O;X1]");	
		queries.put("sulfonic_acid_derivative","[C,$([cX3](:*):*),$([cX2+](:*):*)][S;v6]([!#1!#6;!S])(=O)=O");
		queries.put("phosphine","[C,$([cX3](:*):*),$([cX2+](:*):*)]-[#15;v3X3](-[C,$([cX3](:*):*),$([cX2+](:*):*)])-[C,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("phosphine_imide","[$([#6]P([#6,#1;A])([#6,#1;A])=[#7;X2][#6,#1;A]),$([#6][P;X4+]([#6,#1;A])([#6,#1;A])[#7-][#6,#1;A])]");	
		queries.put("phosphine_sulfide","[$([#6]P([#6,#1;A])([#6,#1;A])=[S;X1]),$([#6][P;X4+]([#16;X1-])([#6,#1;A])[#6,#1;A])]");
		queries.put("phosphinoxide","[#15;A;X4;$([H3]=[OX1]),$([H2](=[OX1])[#6]),$([H1](=[OX1])([#6])[#6]),$([H0](=[OX1])([#6])([#6])[#6])]");	
		queries.put("phosphonium","[P+;!$([P]~[!#6]);!$([P]*~[#7,#8,#15,#16])]");
		queries.put("phosphonic_acid_derivative","[O,N,X]P([O,N,X])([C,$([cX3](:*):*),$([cX2+](:*):*)])=O");	
		queries.put("phosphonic_acid_ester","[#6;!$(C=[O,N,S])]-[#8;X2]P([#6;!$(C=[O,N,S])])([O,N,X])=[O;X1]");
		queries.put("phosphonic_acid_diester","[#6;!$(C=[O,N,S])]-[#8;X2][#15;A;X4;$([H1]),$([H0][#6])](=[O;X1])[#8;X2]-[#6;!$(C=[O,N,S])]");	
		queries.put("phosphonic_amide_ester","[#6;!$(C=[O,N,S])]-[#8;X2][#15;A;X4;$([H1]),$([H0][#6])]([#7;A;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])=[O;X1]");
		queries.put("phosphonic_dihydrazide","[#6]P(=O)([#7]([#6,#1;A])-[#7]([#6,#1;A])[#6,#1;A])[#7]([#6,#1;A])-[#7]([#6,#1;A])[#6,#1;A]");	
		queries.put("trialkylborane","[#6;X4]-[#5](-[#6;X4])-[#6;X4]");
//	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
////		queries.put("","");	
////		queries.put("","");
//
		//  From OpenBabel's FP4
		
		//hits chloromethylenethers and other reactive alkylating agents
		queries.put("Halogen_acetal_like","[NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1]");
		
		// # includes all of the above and other combinations (S-C-N, hydrates, ...), but still no aminomethylenesters and similar
		queries.put("acetal_like","[NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1,NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		// # also reactive alkylating agents. Acid does not have to be carboxylic acid, also S- and P-based acids allowed
		queries.put("halogenmethylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1]");
		
		// # Same as above, but N,O or S instead of halogen. Ester/amide allowed only on one side 

		queries.put("nos_methylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		// # Combination of the last two patterns
		queries.put("hetero_methylen_ester_and_similar","[NX3v3,SX2,OX2;$(**=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1,NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])]");
				
		queries.put("cyanhydrine","[NX1]#[CX2][CX4;$([CH2]),$([CH]([CX2])[#6]),$(C([CX2])([#6])[#6])][OX2H]");
		
		// # includes aminals, silylacetals, ketenesters, etc. C=C DB is not aromatic, everything else may be
		queries.put("Ketenacetal","[#7X2,#8X3,#16X2;$(*[#6,#14])][#6X3]([#7X2,#8X3,#16X2;$(*[#6,#14])])=[#6X3]");
				
		queries.put("spiro_compound","[D4R;$(*(@*)(@*)(@*)@*)]");
		
		// # two different rings sharing exactly two atoms
		queries.put("annelated_rings","[R;$(*(@*)(@*)@*);!$([R2;$(*(@*)(@*)(@*)@*)])]@[R;$(*(@*)(@*)@*);!$([R2;$(*(@*)(@*)(@*)@*)])]");
		
		
		//	# 5 or 6-membered ring containing one O and at least one (r5) or two (r6) oxygen-substituents.		
		queries.put("sugar_pattern_1","[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
 
		// # 5 or 6-membered ring containing one O and an acetal-like bond at postion 2.
		queries.put("sugar_pattern_2","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				 
		// # combination of the two above
		queries.put("sugar_pattern_combi","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_reducing","[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_alpha","[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				
		// # 5 or 6-membered cyclic hemi-acetal
		queries.put("sugar_pattern_2_beta","[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				
		// # pattern1 occours more than once (in same molecule, but moieties don't have to be adjacent!)
		queries.put("poly_sugar_1","([OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)].[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)])");
				
		// # pattern2 occours more than once (in same molecule, but moieties don't have to be adjacent!)
		queries.put("poly_sugar_2","([OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)].[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)])");
				
		queries.put("conjugated_double_bond","*=*[*]=,#,:[*]");
		queries.put("conjugated_tripple_bond","*#*[*]=,#,:[*]");
		
		//		# only one single-bonded substituent on each DB-atom. no aromats. 
		//		# only found when character of DB is explicitely stated.
//		queries.put("cis_double_bond","*/[D2]=[D2]\*");
		//		# only one single-bonded substituent on each DB-atom. no aromats. 
		//		# only found when character of DB is explicitely stated.		
		queries.put("trans_double_bond:","*/[D2]=[D2]/*");
		
		// # should hits all combinations of two acids
		queries.put("mixed_anhydride","[$(*=O),$([#16,#14,#5]),$([#7]([#6]=[OX1]))][#8X2][$(*=O),$([#16,#14,#5]),$([#7]([#6]=[OX1]))]");
		queries.put("halogen_on_hetero","[FX1,ClX1,BrX1,IX1][!#6]");
		queries.put("halogen_multi_subst","[F,Cl,Br,I;!$([X1]);!$([X0-])]");
		
		// # 1,3 migration of H allowed. Includes keto/enol and amide/enamide. 
		// # Aromatic rings must stay aromatic - no keto form of phenol 
		queries.put("_13-tautomerizable","[$([#7X2,OX1,SX1]=*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=*),$([#7X3,OX2,SX2;!H0]*:n)]");
		queries.put("_15-tautomerizable","[$([#7X2,OX1,SX1]=,:**=,:*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=**=*),$([#7X3,OX2,SX2;!H0]*=,:**:n)]");
		
		// # Hits atoms with tetrahedral chirality, if chiral center is specified in the SMILES string
		// # depictmach does not find oxonium, sulfonium, or sulfoxides!
		queries.put("chiral_center_specified","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*)]");
		
		// # Hits atoms with tetrahedral chirality, if chiral center is not specified in the SMILES string
		// # "@?" (unspecified chirality) is not yet supported in Open Babel Version 2.0 
		queries.put("chiral_center_unspecified","[$([*@?](~*)(~*)(*)*),$([*@?H](*)(*)*),$([*@?](~*)(*)*),$([*@?H](~*)~*)]");
//		queries.put("","");
//		queries.put("","");
//		queries.put("","");
//		queries.put("","");

		queries.put("pattern_c51", "ccc(c)O");
		queries.put("pattern_c52", "cccccO");	
		queries.put("pattern_c53", "cc(c)O");	
		queries.put("pattern_c54", "ccccO");	
		queries.put("pattern_c55", "ccccccO");	
		queries.put("pattern_c56", "cccc(O)cc");	
		queries.put("pattern_c57", "ccccc(c)O");	
		queries.put("pattern_c58", "Oc1ccccc1");	
		queries.put("pattern_c59", "cccc(c)O");	
		queries.put("pattern_c60", "ccnc");	
		queries.put("pattern_c61", "cnc");	
		queries.put("pattern_c62", "ccn");	
		queries.put("pattern_c63", "cccnc");	
		queries.put("pattern_c64", "cccn");	
		queries.put("pattern_c65", "ccN");	
		queries.put("pattern_c66", "cccN");	
		queries.put("pattern_c67", "ccC");	
		queries.put("pattern_c68", "Cc1ccccc1");	
		queries.put("pattern_c69", "cccC");
		queries.put("pattern_c70", "ccc(C)cc");	
		queries.put("pattern_c71", "cccc");	
		queries.put("pattern_c72", "ccccC");	
		queries.put("pattern_c73", "ccccc");	
		queries.put("pattern_c74", "cc(c)C");
		queries.put("pattern_c75", "ccc");	
		queries.put("pattern_c76", "c1ccccc1");	
		queries.put("pattern_c77", "cccccc");
		queries.put("pattern_c78", "cccc(C)cc");	
		queries.put("pattern_c79", "ccc(c)C");	
		queries.put("pattern_c80", "cccc(c)C");
		queries.put("pattern_c81", "cccccC");	
		queries.put("pattern_c82", "ccccc(c)C");	
		queries.put("pattern_c83", "ccccccC");
		queries.put("pattern_c84", "CCN");	
		queries.put("pattern_c85", "CNC");	
		queries.put("pattern_c86", "CCC");
		
		// CYP2A6
		queries.put("pattern_c87", "CCNC");	
		queries.put("pattern_c88", "CCCC");	
		queries.put("pattern_c89", "CCO");
		queries.put("pattern_c90", "cccccn");	
		queries.put("pattern_c91", "ccncc");	
		queries.put("pattern_c92", "ccccnc");
		queries.put("pattern_c93", "ccccn");	
		queries.put("pattern_c94", "c1ccncc1");	
		queries.put("pattern_c95", "cccncc");
		
		// CYP2B6
		queries.put("pattern_c96", "CCc1ccccc1");	
		queries.put("pattern_c97", "cccc(cc)CC");	
		queries.put("pattern_c98", "ccccccCC");
		queries.put("pattern_c99", "ccccc(c)CC");	
		queries.put("pattern_c100", "ccc(cc)CC");	
		queries.put("pattern_c101", "cccc(c)CC");
		queries.put("pattern_c102", "cccccCC");	
		queries.put("pattern_c103", "ccc(c)CC");	
		queries.put("pattern_c104", "ccccCC");
		queries.put("pattern_c105", "cc(c)CC");	
		queries.put("pattern_c106", "ccc(O)cc");	
		queries.put("pattern_c107", "NC=O");
		queries.put("pattern_c108", "cCN");	
		queries.put("pattern_c109", "ccCC");	
		queries.put("pattern_c110", "ccCN");		

		// MACCS 322, based on the definitions of the MDL keys, and what I found on the 
		// MayaChem tool website: http://www.mayachemtools.org/docs/modules/html/MACCSKeys.html
		
		queries.put("maccs_322_001", "*~*~*~*");
		queries.put("maccs_322_002", "[!#1!#6]");	
//		queries.put("maccs_322_003", "");	
		queries.put("maccs_322_004", "*~*~*~*~*");	
		queries.put("maccs_322_005", "[$(*~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6])]");
		queries.put("maccs_322_006", "[$(*~[!#1!#6]~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6]~[!#1!#6])]");	
		queries.put("maccs_322_007", "[!#1!#6][H]");	
		queries.put("maccs_322_008", "[$(*~*~[#6]([H])[H]),$(*~C(~*)([H])[H])]");		
		queries.put("maccs_322_009", "*~C([H])([H])[H]");
//		queries.put("maccs_322_010", "[F,Cl,Br,I]");	//halogen (already here)
		queries.put("maccs_322_011", "[$(*-*(-*)-*)]");	 // [$(),$(*-*(-*)-*)]
//		queries.put("maccs_322_012", "");	
		queries.put("maccs_322_013", "[$(*@-*(@-*)@-*)]"); // [$(),$(*@-*(@-*)@-*)]
		queries.put("maccs_322_014", "*@-*!@-*@-*");	
		queries.put("maccs_322_015", "*!:*:*!:*");	
		queries.put("maccs_322_016", "*!@-*!@-*");
		queries.put("maccs_322_017", "*!@-*@-*!@-*");
		queries.put("maccs_322_018", "*:*!:*:*");	
//		queries.put("maccs_322_019", "");	// Same as organoheterocycle
		queries.put("maccs_322_020", "[$(*-*(-*)(-*)(-*)-*),$(*@-*(@-*)(@-*)(@-*)@-*),$([!#1!#6!#7!#8!#16!#9!#17!#35!#53])]");	
		queries.put("maccs_322_021", "[!+0]");
//		queries.put("maccs_322_022", ""); // Same as nitrogen_atom	
//		queries.put("maccs_322_023", ""); // Same as sulfur_atom	
//		queries.put("maccs_322_024", ""); // Same as oxygen_atom
//		queries.put("maccs_322_025", "");
//		queries.put("maccs_322_026", "");	
		
		queries.put("maccs_322_027", "C(CC)");	
		queries.put("maccs_322_028", "C(CCC)");	
		queries.put("maccs_322_029", "C(CN)");
		queries.put("maccs_322_030", "C(CCN)");	
		queries.put("maccs_322_031", "C(NN)");	
		queries.put("maccs_322_032", "C(NNC)");
		queries.put("maccs_322_033", "C(NNN)");	
		queries.put("maccs_322_034", "C(CO)");	
		queries.put("maccs_322_035", "C(CCO)");
		queries.put("maccs_322_036", "C(NO)");
		queries.put("maccs_322_037", "C(NCO)");	
		queries.put("maccs_322_038", "C(NNO)");	
		queries.put("maccs_322_039", "C(OO)");	
		queries.put("maccs_322_040", "C(COO)");
		queries.put("maccs_322_041", "C(NOO)");	
		queries.put("maccs_322_042", "C(OOO)");	
		queries.put("maccs_322_043", "[!#1!#6](CC)");
		queries.put("maccs_322_044", "[!#1!#6](CCC)");	
		queries.put("maccs_322_045", "[!#1!#6](CN)");	
		queries.put("maccs_322_046", "[!#1!#6](CCN)");
		queries.put("maccs_322_047", "[!#1!#6](NN)");
		queries.put("maccs_322_048", "[!#1!#6](CNN)");	
		queries.put("maccs_322_049", "[!#1!#6](NNN)");	
		queries.put("maccs_322_050", "[!#1!#6](CO)");	
		queries.put("maccs_322_051", "[!#1!#6](CCO)");
		queries.put("maccs_322_052", "[!#1!#6](NO)");	
		queries.put("maccs_322_053", "[!#1!#6](CNO)");	
		queries.put("maccs_322_054", "[!#1!#6](NNO)");
		queries.put("maccs_322_055", "[!#1!#6](OO)");	
		queries.put("maccs_322_056", "[!#1!#6](COO)");	
		queries.put("maccs_322_057", "[!#1!#6](NOO)");
		queries.put("maccs_322_058", "[!#1!#6](OOO)");
		queries.put("maccs_322_059", "C-C");	
		queries.put("maccs_322_060", "C-N");	
		queries.put("maccs_322_061", "C-O");	
		queries.put("maccs_322_062", "C-S");
		queries.put("maccs_322_063", "C-[#17]");	
		queries.put("maccs_322_064", "C-P");	
		queries.put("maccs_322_065", "C-F");
		queries.put("maccs_322_066", "C-[Br]");
		queries.put("maccs_322_067", "C-[Si]");
		queries.put("maccs_322_068", "C-I");
		queries.put("maccs_322_069", "C-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_070", "N-N");
		queries.put("maccs_322_071", "N-O");
		queries.put("maccs_322_072", "N-S");
		queries.put("maccs_322_073", "N-[#17]");
		queries.put("maccs_322_074", "N-P");
		queries.put("maccs_322_075", "N-F");
		queries.put("maccs_322_076", "N-[Br]");
		queries.put("maccs_322_077", "N-[Si]");
		queries.put("maccs_322_078", "N-I");
		queries.put("maccs_322_079", "N-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_080", "O-O");
		queries.put("maccs_322_081", "O-S");
		queries.put("maccs_322_082", "O-[#17]");
		queries.put("maccs_322_083", "O-P");
		queries.put("maccs_322_084", "O-F");
		queries.put("maccs_322_085", "O-[Br]");
		queries.put("maccs_322_086", "O-[Si]");
		queries.put("maccs_322_087", "O-I");
		queries.put("maccs_322_088", "O-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_089", "S-S");
		queries.put("maccs_322_090", "S-[#17]");
		queries.put("maccs_322_091", "S-P");
		queries.put("maccs_322_092", "S-F");
		queries.put("maccs_322_093", "S-[Br]");
		queries.put("maccs_322_094", "S-[Si]");
		queries.put("maccs_322_095", "S-I");
		queries.put("maccs_322_096", "S-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_097", "[#17]-[#17]");
		queries.put("maccs_322_098", "[#17]-P");
		queries.put("maccs_322_099", "[#17]-F");
		queries.put("maccs_322_100", "[#17]-[Br]");
		queries.put("maccs_322_101", "[#17]-[Si]");
		queries.put("maccs_322_102", "[#17]-I");
		queries.put("maccs_322_103", "[#17]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_104", "P-P");
		queries.put("maccs_322_105", "P-F");
		queries.put("maccs_322_106", "P-[Br]");
		queries.put("maccs_322_107", "P-[Si]");
		queries.put("maccs_322_108", "P-I");
		queries.put("maccs_322_109", "P-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_110", "F-F");
		queries.put("maccs_322_111", "F-[Br]");
		queries.put("maccs_322_112", "F-[Si]");
		queries.put("maccs_322_113", "F-I");
		queries.put("maccs_322_114", "F-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_115", "[Br]-[Br]");
		queries.put("maccs_322_116", "[Br]-[Si]");
		queries.put("maccs_322_117", "[Br]-I");
		queries.put("maccs_322_118", "[Br]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_119", "[Si]-[Si]");
		queries.put("maccs_322_120", "[Si]-I");
		queries.put("maccs_322_121", "[Si]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_122", "I-I");
		queries.put("maccs_322_123", "I-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_124", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_125", "C=C");
		queries.put("maccs_322_126", "C=N");
		queries.put("maccs_322_127", "C=O");
		queries.put("maccs_322_128", "C=S");
		queries.put("maccs_322_129", "C=[#17]");
		queries.put("maccs_322_130", "C=P");
		queries.put("maccs_322_131", "C=F");
		queries.put("maccs_322_132", "C=[Br]");
		queries.put("maccs_322_133", "C=[Si]");
		queries.put("maccs_322_134", "C=I");
		queries.put("maccs_322_135", "C=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_136", "N=N");
		queries.put("maccs_322_137", "N=O");
		queries.put("maccs_322_138", "N=S");
		queries.put("maccs_322_139", "N=O");
		queries.put("maccs_322_140", "N=[#17]");
		queries.put("maccs_322_141", "N=P");
		queries.put("maccs_322_142", "N=F");
		queries.put("maccs_322_143", "N=[Br]");
		queries.put("maccs_322_144", "N=[Si]");
		queries.put("maccs_322_145", "N=I");
		queries.put("maccs_322_146", "N=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_147", "O=O");
		queries.put("maccs_322_148", "O=S");
		queries.put("maccs_322_149", "O=[#17]");
		queries.put("maccs_322_150", "O=P");
		queries.put("maccs_322_151", "O=[Br]");
		queries.put("maccs_322_152", "O=[Si]");
		queries.put("maccs_322_153", "O=I");
		queries.put("maccs_322_154", "O=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_155", "S=S");
		queries.put("maccs_322_156", "S=[#17]");
		queries.put("maccs_322_157", "S=P");
		queries.put("maccs_322_158", "S=F");
		queries.put("maccs_322_159", "S=[Br]");
		queries.put("maccs_322_160", "S=[Si]");
		queries.put("maccs_322_161", "S=I");
		queries.put("maccs_322_162", "S=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_163", "[#17]=[#17]");
		queries.put("maccs_322_164", "[#17]=P");
		queries.put("maccs_322_165", "[#17]=[#9]");
		queries.put("maccs_322_166", "[#17]=[Br]");
		queries.put("maccs_322_167", "[#17]=[Si]");
		queries.put("maccs_322_168", "[#17]=I");
		queries.put("maccs_322_169", "[#17]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_170", "P=P");
		queries.put("maccs_322_171", "P=[#9]");
		queries.put("maccs_322_172", "P=[Br]");
		queries.put("maccs_322_173", "P=[Si]");
		queries.put("maccs_322_174", "P=I");
		queries.put("maccs_322_175", "P=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_176", "[#9]=[#9]");
		queries.put("maccs_322_177", "[#9]=[Br]");
		queries.put("maccs_322_178", "[#9]=[Si]");
		queries.put("maccs_322_179", "[#9]=I");
		queries.put("maccs_322_180", "[#9]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_181", "[Br]=[Br]");
		queries.put("maccs_322_182", "[Br]=[Si]");
		queries.put("maccs_322_183", "[Br]=I");
		queries.put("maccs_322_184", "[Br]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_185", "[Si]=[Si]");
		queries.put("maccs_322_186", "[Si]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_187", "I-I");
		queries.put("maccs_322_188", "I=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_189", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]=[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_190", "C#C");
		queries.put("maccs_322_191", "C#N");
		queries.put("maccs_322_192", "C#O");
		queries.put("maccs_322_193", "C#S");
		queries.put("maccs_322_194", "C#S");
		queries.put("maccs_322_195", "C#[#17]");		
		queries.put("maccs_322_196", "C#P");
		queries.put("maccs_322_197", "C#F");
		queries.put("maccs_322_198", "C#[Br]");
		queries.put("maccs_322_199", "C#[Si]");			
		queries.put("maccs_322_200", "C#I");
		queries.put("maccs_322_201", "C#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_202", "N#N");
		queries.put("maccs_322_203", "N#O");	
		queries.put("maccs_322_204", "N#S");	
		queries.put("maccs_322_205", "N#[#17]");	
		queries.put("maccs_322_206", "N#P");	
		queries.put("maccs_322_207", "N#F");	
		queries.put("maccs_322_208", "N#[Br]");	
		queries.put("maccs_322_209", "N#[Si]");	
		queries.put("maccs_322_210", "N#I");
		queries.put("maccs_322_211", "N#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");	
		queries.put("maccs_322_212", "O#O");	
		queries.put("maccs_322_213", "O#S");	
		queries.put("maccs_322_214", "O#[#17]");	
		queries.put("maccs_322_215", "O#P");	
		queries.put("maccs_322_216", "O#F");	
		queries.put("maccs_322_217", "O#[Br]");	
		queries.put("maccs_322_218", "O#[Si]");	
		queries.put("maccs_322_219", "O#I");	
		queries.put("maccs_322_220", "O#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");	
		queries.put("maccs_322_221", "S#S");		
		queries.put("maccs_322_222", "S#[#17]");		
		queries.put("maccs_322_223", "S#P");
		queries.put("maccs_322_224", "S#[#9]");
		queries.put("maccs_322_225", "S#[Br]");
		queries.put("maccs_322_226", "S#[Si]");
		queries.put("maccs_322_227", "S#I");
		queries.put("maccs_322_228", "S#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_229", "[#17]#[#17]");
		queries.put("maccs_322_230", "[#17]#[#15]");
		queries.put("maccs_322_231", "[#17]#[#9]");
		queries.put("maccs_322_232", "[#17]#[Br]");
		queries.put("maccs_322_233", "[#17]#[Si]");
		queries.put("maccs_322_234", "[#17]#I");
		queries.put("maccs_322_235", "[#17]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_236", "[#15]#[#15]");
		queries.put("maccs_322_237", "[#15]#[#9]");
		queries.put("maccs_322_238", "[#15]#[Br]");
		queries.put("maccs_322_239", "[#15]#[Si]");
		queries.put("maccs_322_240", "[#15]#[#53]");
		queries.put("maccs_322_241", "[#15]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_242", "[#9]#[#9]");
		queries.put("maccs_322_243", "[#9]#[Br]");
		queries.put("maccs_322_244", "[#9]#[Si]");
		queries.put("maccs_322_245", "[#9]#[#53]");
		queries.put("maccs_322_246", "[#9]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_247", "[#35]#[Br]");
		queries.put("maccs_322_248", "[#35]#[Si]");
		queries.put("maccs_322_249", "[#35]#[#53]");
		queries.put("maccs_322_250", "[#35]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_251", "[Si]#[Si]");		
		queries.put("maccs_322_252", "[Si]#[#53]");
		queries.put("maccs_322_253", "[Si]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_254", "[#53]#[#53]");
		queries.put("maccs_322_255", "[#53]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_256", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]#[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_257", "[#6]@-[#6]");		
		queries.put("maccs_322_258", "[#6]@-[#7]");
		queries.put("maccs_322_259", "[#6]@-[#8]");		
		queries.put("maccs_322_260", "[#6]@-[#16]");
		queries.put("maccs_322_261", "[#6]@-[#17]");		
		queries.put("maccs_322_262", "[#6]@-[#15]");		
		queries.put("maccs_322_263", "[#6]@-[#9]");
		queries.put("maccs_322_264", "[#6]@-[#35]");		
		queries.put("maccs_322_265", "[#6]@-[#14]");		
		queries.put("maccs_322_266", "[#6]@-[#53]");
		queries.put("maccs_322_267", "[#6]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_268", "[#7]@-[#7]");		
		queries.put("maccs_322_269", "[#7]@-[#8]");
		queries.put("maccs_322_270", "[#7]@-[#16]");
		queries.put("maccs_322_271", "[#7]@-[#17]");
		queries.put("maccs_322_272", "[#7]@-[#15]");		
		queries.put("maccs_322_273", "[#7]@-[#9]");		
		queries.put("maccs_322_274", "[#7]@-[#35]");
		queries.put("maccs_322_275", "[#7]@-[#14]");		
		queries.put("maccs_322_276", "[#7]@-[#53]");
		queries.put("maccs_322_277", "[#7]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_278", "[#8]@-[#8]");		
		queries.put("maccs_322_279", "[#8]@-[#16]");
		queries.put("maccs_322_280", "[#8]@-[#17]");		
		queries.put("maccs_322_281", "[#8]@-[#15]");
		queries.put("maccs_322_282", "[#8]@-[#9]");		
		queries.put("maccs_322_283", "[#8]@-[#35]");		
		queries.put("maccs_322_284", "[#8]@-[#14]");
		queries.put("maccs_322_285", "[#8]@-[#53]");		
		queries.put("maccs_322_286", "[#8]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_287", "[#16]@-[#16]");
		queries.put("maccs_322_288", "[#16]@-[#17]");		
		queries.put("maccs_322_289", "[#16]@-[#15]");
		queries.put("maccs_322_290", "[#16]@-[#9]");
		queries.put("maccs_322_291", "[#16]@-[#35]");
		queries.put("maccs_322_292", "[#16]@-[#14]");		
		queries.put("maccs_322_293", "[#16]@-[#53]");		
		queries.put("maccs_322_294", "[#16]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_295", "[#17]@-[#17]");		
		queries.put("maccs_322_296", "[#17]@-[#15]");		
		queries.put("maccs_322_297", "[#17]@-[#9]");
		queries.put("maccs_322_298", "[#17]@-[#35]");		
		queries.put("maccs_322_299", "[#17]@-[#14]");
		queries.put("maccs_322_300", "[#17]@-[#53]");		
		queries.put("maccs_322_301", "[#17]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_302", "[#15]@-[#15]");
		queries.put("maccs_322_303", "[#15]@-[#9]");		
		queries.put("maccs_322_304", "[#15]@-[#35]");
		queries.put("maccs_322_305", "[#15]@-[#14]");
		queries.put("maccs_322_306", "[#15]@-[#53]");		
		queries.put("maccs_322_307", "[#15]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_308", "[#9]@-[#9]");		
		queries.put("maccs_322_309", "[#9]@-[#35]");
		queries.put("maccs_322_310", "[#9]@-[#14]");
		queries.put("maccs_322_311", "[#9]@-[#53]");		
		queries.put("maccs_322_312", "[#9]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_313", "[#35]@-[#35]");		
		queries.put("maccs_322_314", "[#35]@-[#14]");
		queries.put("maccs_322_315", "[#35]@-[#9]");
		queries.put("maccs_322_316", "[#35]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");
		queries.put("maccs_322_317", "[#14]@-[#14]");		
		queries.put("maccs_322_318", "[#14]@-[#53]");
		queries.put("maccs_322_319", "[#14]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_320", "[#53]@-[#53]");
		queries.put("maccs_322_321", "[#53]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");		
		queries.put("maccs_322_322", "[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]@-[!#1!#6!#7!#8!#14!#15!#16!#9!#17!#35!#53]");

		
		// PubChem keys 264 to 881 (The index xxx in the pubchem_xxx is based on a start from 0)
		// ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
		// The first 263 keys cannot be mapped to an atom (e.g. 261	>= 4 aromatic rings)
		
		queries.put("pubchem_263","[Li]-[H]");
		queries.put("pubchem_264","[Li]-[Li]");
		queries.put("pubchem_265","[Li]-B");
		queries.put("pubchem_266","[Li]-C");
		queries.put("pubchem_267","[Li]-O");
		queries.put("pubchem_268","[Li]-F");
		queries.put("pubchem_269","[Li]-P");
		queries.put("pubchem_270","[Li]-S");
		queries.put("pubchem_271","[Li]-[Cl]");
		queries.put("pubchem_272","B-[H]");
		queries.put("pubchem_273","B-B");
		queries.put("pubchem_274","B-C");
		queries.put("pubchem_275","B-N");
		queries.put("pubchem_276","B-O");
		queries.put("pubchem_277","B-F");
		queries.put("pubchem_278","B-[Si]");
		queries.put("pubchem_279","B-P");
		queries.put("pubchem_280","B-S");
		queries.put("pubchem_281","B-[Cl]");
		queries.put("pubchem_282","B-[Br]");
		queries.put("pubchem_283","C-[H]");
		queries.put("pubchem_284","C-C");
		queries.put("pubchem_285","C-N");
		queries.put("pubchem_286","C-O");
		queries.put("pubchem_287","C-F");
		queries.put("pubchem_288","C-[Na]");
		queries.put("pubchem_289","C-[Mg]");
		queries.put("pubchem_290","C-[Al]");
		queries.put("pubchem_291","C-[Si]");
		queries.put("pubchem_292","C-P");
		queries.put("pubchem_293","C-S");
		queries.put("pubchem_294","C-[Cl]");
		queries.put("pubchem_295","C-[As]");
		queries.put("pubchem_296","C-[Se]");
		queries.put("pubchem_297","C-[Br]");
		queries.put("pubchem_298","C-I");
		queries.put("pubchem_299","N-[H]");
		queries.put("pubchem_300","N-N");
		queries.put("pubchem_301","N-O");
		queries.put("pubchem_302","N-F");
		queries.put("pubchem_303","N-[Si]");
		queries.put("pubchem_304","N-P");
		queries.put("pubchem_305","N-S");
		queries.put("pubchem_306","N-[Cl]");
		queries.put("pubchem_307","N-[Br]");
		queries.put("pubchem_308","O-[H]");
		queries.put("pubchem_309","[O]-[O]");
		queries.put("pubchem_310","O-[Mg]");
		queries.put("pubchem_311","O-[Na]");
		queries.put("pubchem_312","O-[Al]");
		queries.put("pubchem_313","O-[Si]");
		queries.put("pubchem_314","O-P");
		queries.put("pubchem_315","O-[K]");
		queries.put("pubchem_316","F-P");
		queries.put("pubchem_317","F-S");
		queries.put("pubchem_318","[Al]-[H]");
		queries.put("pubchem_319","[Al]-[Cl]");
		queries.put("pubchem_320","[Si]-[H]");
		queries.put("pubchem_321","[Si]-[Si]");
		queries.put("pubchem_322","[Si]-[Cl]");
		queries.put("pubchem_323","P-[H]");
		queries.put("pubchem_324","P-P");
		queries.put("pubchem_325","[As]-[H]");
		queries.put("pubchem_326","[As]-[As]");
		queries.put("pubchem_327","C(~[Br])(~C)");
		queries.put("pubchem_328","C(~[Br])(~C)(~C)");
		queries.put("pubchem_329","C(~[Br])(~[H])");
		queries.put("pubchem_330","C(~[Br])(:C)");
		queries.put("pubchem_331","C(~[Br])(:N)");
		queries.put("pubchem_332","C(~C)(~C)");
		queries.put("pubchem_333","C(~C)(~C)(~C)");
		queries.put("pubchem_334","C(~C)(~C)(~C)(~C)");
		queries.put("pubchem_335","C(~C)(~C)(~C)(~[H])");
		queries.put("pubchem_336","C(~C)(~C)(~C)(~N)");
		queries.put("pubchem_337","C(~C)(~C)(~C)(~O)");
		queries.put("pubchem_338","C(~C)(~C)(~[H])(~N)");
		queries.put("pubchem_339","C(~C)(~C)(~[H])(~O)");
		queries.put("pubchem_340","C(~C)(~C)(~N)");
		queries.put("pubchem_341","C(~C)(~C)(~O)");
		queries.put("pubchem_342","C(~C)(~[Cl])");
		queries.put("pubchem_343","C(~C)(~[Cl])(~[H])");
		queries.put("pubchem_344","C(~C)(~[H])");
		queries.put("pubchem_345","C(~C)(~[H])(~N)");
		queries.put("pubchem_346","C(~C)(~[H])(~O)");
		queries.put("pubchem_347","C(~C)(~[H])(~O)(~O)");
		queries.put("pubchem_348","C(~C)(~[H])(~P)");
		queries.put("pubchem_349","C(~C)(~[H])(~S)");
		queries.put("pubchem_350","C(~C)(~I)");
		queries.put("pubchem_351","C(~C)(~N)");
		queries.put("pubchem_352","C(~C)(~O)");
		queries.put("pubchem_353","C(~C)(~S)");
		queries.put("pubchem_354","C(~C)(~[Si])");
		queries.put("pubchem_355","C(~C)(:C)");
		queries.put("pubchem_356","C(~C)(:C)(:C)");
		queries.put("pubchem_357","C(~C)(:C)(:N)");
		queries.put("pubchem_358","C(~C)(:N)");
		queries.put("pubchem_359","C(~C)(:N)(:N)");
		queries.put("pubchem_360","C(~[Cl])(~[Cl])");
		queries.put("pubchem_361","C(~[Cl])(~[H])");
		queries.put("pubchem_362","C(~[Cl])(:C)");
		queries.put("pubchem_363","C(~F)(~F)");
		queries.put("pubchem_364","C(~F)(:C)");
		queries.put("pubchem_365","C(~[H])(~N)");
		queries.put("pubchem_366","C(~[H])(~O)");
		queries.put("pubchem_367","C(~[H])(~O)(~O)");
		queries.put("pubchem_368","C(~[H])(~S)");
		queries.put("pubchem_369","C(~[H])(~[Si])");
		queries.put("pubchem_370","C(~[H])(:C)");
		queries.put("pubchem_371","C(~[H])(:C)(:C)");
		queries.put("pubchem_372","C(~[H])(:C)(:N)");
		queries.put("pubchem_373","C(~[H])(:N)");
		queries.put("pubchem_374","C(~[H])(~[H])(~[H])");
		queries.put("pubchem_375","C(~N)(~N)");
		queries.put("pubchem_376","C(~N)(:C)");
		queries.put("pubchem_377","C(~N)(:C)(:C)");
		queries.put("pubchem_378","C(~N)(:C)(:N)");
		queries.put("pubchem_379","C(~N)(:N)");
		queries.put("pubchem_380","C(~O)(~O)");
		queries.put("pubchem_381","C(~O)(:C)");
		queries.put("pubchem_382","C(~O)(:C)(:C)");
		queries.put("pubchem_383","C(~S)(:C)");
		queries.put("pubchem_384","C(:C)(:C)");
		queries.put("pubchem_385","C(:C)(:C)(:C)");
		queries.put("pubchem_386","C(:C)(:C)(:N)");
		queries.put("pubchem_387","C(:C)(:N)");
		queries.put("pubchem_388","C(:C)(:N)(:N)");
		queries.put("pubchem_389","C(:N)(:N)");
		queries.put("pubchem_390","N(~C)(~C)");
		queries.put("pubchem_391","N(~C)(~C)(~C)");
		queries.put("pubchem_392","N(~C)(~C)(~[H])");
		queries.put("pubchem_393","N(~C)(~[H])");
		queries.put("pubchem_394","N(~C)(~[H])(~N)");
		queries.put("pubchem_395","N(~C)(~O)");
		queries.put("pubchem_396","N(~C)(:C)");
		queries.put("pubchem_397","N(~C)(:C)(:C)");
		queries.put("pubchem_398","N(~[H])(~N)");
		queries.put("pubchem_399","N(~[H])(:C)");
		queries.put("pubchem_400","N(~[H])(:C)(:C)");
		queries.put("pubchem_401","N(~O)(~O)");
		queries.put("pubchem_402","N(~O)(:O)");
		queries.put("pubchem_403","N(:C)(:C)");
		queries.put("pubchem_404","N(:C)(:C)(:C)");
		queries.put("pubchem_405","O(~C)(~C)");
		queries.put("pubchem_406","O(~C)(~[H])");
		queries.put("pubchem_407","O(~C)(~P)");
		queries.put("pubchem_408","O(~[H])(~S)");
		queries.put("pubchem_409","O(:C)(:C)");
		queries.put("pubchem_410","P(~C)(~C)");
		queries.put("pubchem_411","P(~O)(~O)");
		queries.put("pubchem_412","S(~C)(~C)");
		queries.put("pubchem_413","S(~C)(~[H])");
		queries.put("pubchem_414","S(~C)(~O)");
		queries.put("pubchem_415","[Si](~C)(~C)");
		queries.put("pubchem_416","C=C");
		queries.put("pubchem_417","C#C");
		queries.put("pubchem_418","C=N");
		queries.put("pubchem_419","C#N");
		queries.put("pubchem_420","C=O");
		queries.put("pubchem_421","C=S");
		queries.put("pubchem_422","N=N");
		queries.put("pubchem_423","N=O");
		queries.put("pubchem_424","N=P");
		queries.put("pubchem_425","P=O");
		queries.put("pubchem_426","P=P");
		queries.put("pubchem_427","C(#C)(-C)");
		queries.put("pubchem_428","C(#C)(-[H])");
		queries.put("pubchem_429","C(#N)(-C)");
		queries.put("pubchem_430","C(-C)(-C)(=C)");
		queries.put("pubchem_431","C(-C)(-C)(=N)");
		queries.put("pubchem_432","C(-C)(-C)(=O)");
		queries.put("pubchem_433","C(-C)(-[Cl])(=O)");
		queries.put("pubchem_434","C(-C)(-[H])(=C)");
		queries.put("pubchem_435","C(-C)(-[H])(=N)");
		queries.put("pubchem_436","C(-C)(-[H])(=O)");
		queries.put("pubchem_437","C(-C)(-N)(=C)");
		queries.put("pubchem_438","C(-C)(-N)(=N)");
		queries.put("pubchem_439","C(-C)(-N)(=O)");
		queries.put("pubchem_440","C(-C)(-O)(=O)");
		queries.put("pubchem_441","C(-C)(=C)");
		queries.put("pubchem_442","C(-C)(=N)");
		queries.put("pubchem_443","C(-C)(=O)");
		queries.put("pubchem_444","C(-[Cl])(=O)");
		queries.put("pubchem_445","C(-[H])(-N)(=C)");
		queries.put("pubchem_446","C(-[H])(=C)");
		queries.put("pubchem_447","C(-[H])(=N)");
		queries.put("pubchem_448","C(-[H])(=O)");
		queries.put("pubchem_449","C(-N)(=C)");
		queries.put("pubchem_450","C(-N)(=N)");
		queries.put("pubchem_451","C(-N)(=O)");
		queries.put("pubchem_452","C(-O)(=O)");
		queries.put("pubchem_453","N(-C)(=C)");
		queries.put("pubchem_454","N(-C)(=O)");
		queries.put("pubchem_455","N(-O)(=O)");
		queries.put("pubchem_456","P(-O)(=O)");
		queries.put("pubchem_457","S(-C)(=O)");
		queries.put("pubchem_458","S(-O)(=O)");
		queries.put("pubchem_459","S(=O)(=O)");
		queries.put("pubchem_460","C-C-C#C");
		queries.put("pubchem_461","O-C-C=N");
		queries.put("pubchem_462","O-C-C=O");
		queries.put("pubchem_463","N:C-S-[#1]");
		queries.put("pubchem_464","N-C-C=C");
		queries.put("pubchem_465","O=S-C-C");
		queries.put("pubchem_466","N#C-C=C");
		queries.put("pubchem_467","C=N-N-C");
		queries.put("pubchem_468","O=S-C-N");
		queries.put("pubchem_469","S-S-C:C");
		queries.put("pubchem_470","C:C-C=C");
		queries.put("pubchem_471","S:C:C:C");
		queries.put("pubchem_472","C:N:C-C");
		queries.put("pubchem_473","S-C:N:C");
		queries.put("pubchem_474","S:C:C:N");
		queries.put("pubchem_475","S-C=N-C");
		queries.put("pubchem_476","C-O-C=C");
		queries.put("pubchem_477","N-N-C:C");
		queries.put("pubchem_478","S-C=N-[#1]");
		queries.put("pubchem_479","S-C-S-C");
		queries.put("pubchem_480","C:S:C-C");
		queries.put("pubchem_481","O-S-C:C");
		queries.put("pubchem_482","C:N-C:C");
		queries.put("pubchem_483","N-S-C:C");
		queries.put("pubchem_484","N-C:N:C");
		queries.put("pubchem_485","N:C:C:N");
		queries.put("pubchem_486","N-C:N:N");
		queries.put("pubchem_487","N-C=N-C");
		queries.put("pubchem_488","N-C=N-[#1]");
		queries.put("pubchem_489","N-C-S-C");
		queries.put("pubchem_490","C-C-C=C");
		queries.put("pubchem_491","C-N:C-[#1]");
		queries.put("pubchem_492","N-C:O:C");
		queries.put("pubchem_493","O=C-C:C");
		queries.put("pubchem_494","O=C-C:N");
		queries.put("pubchem_495","C-N-C:C");
		queries.put("pubchem_496","N:N-C-[#1]");
		queries.put("pubchem_497","O-C:C:N");
		queries.put("pubchem_498","O-C=C-C");
		queries.put("pubchem_499","N-C:C:N");
		queries.put("pubchem_500","C-S-C:C");
		queries.put("pubchem_501","[Cl]-C:C-C");
		queries.put("pubchem_502","N-C=C-[#1]");
		queries.put("pubchem_503","[Cl]-C:C-[#1]");
		queries.put("pubchem_504","N:C:N-C");
		queries.put("pubchem_505","[Cl]-C:C-O");
		queries.put("pubchem_506","C-C:N:C");
		queries.put("pubchem_507","C-C-S-C");
		queries.put("pubchem_508","S=C-N-C");
		queries.put("pubchem_509","[Br]-C:C-C");
		queries.put("pubchem_510","[#1]-N-N-[#1]");
		queries.put("pubchem_511","S=C-N-[#1]");
		queries.put("pubchem_512","C-[As]-O-[#1]");
		queries.put("pubchem_513","S:C:C-[#1]");
		queries.put("pubchem_514","O-N-C-C");
		queries.put("pubchem_515","N-N-C-C");
		queries.put("pubchem_516","[#1]-C=C-[#1]");
		queries.put("pubchem_517","N-N-C-N");
		queries.put("pubchem_518","O=C-N-N");
		queries.put("pubchem_519","N=C-N-C");
		queries.put("pubchem_520","C=C-C:C");
		queries.put("pubchem_521","C:N-C-[#1]");
		queries.put("pubchem_522","C-N-N-[#1]");
		queries.put("pubchem_523","N:C:C-C");
		queries.put("pubchem_524","C-C=C-C");
		queries.put("pubchem_525","[As]-C:C-[#1]");
		queries.put("pubchem_526","[Cl]-C:C-Cl");
		queries.put("pubchem_527","C:C:N-[#1]");
		queries.put("pubchem_528","[#1]-N-C-[#1]");
		queries.put("pubchem_529","[Cl]-C-C-Cl");
		queries.put("pubchem_530","N:C-C:C");
		queries.put("pubchem_531","S-C:C-C");
		queries.put("pubchem_532","S-C:C-[#1]");
		queries.put("pubchem_533","S-C:C-N");
		queries.put("pubchem_534","S-C:C-O");
		queries.put("pubchem_535","O=C-C-C");
		queries.put("pubchem_536","O=C-C-N");
		queries.put("pubchem_537","O=C-C-O");
		queries.put("pubchem_538","N=C-C-C");
		queries.put("pubchem_539","N=C-C-[#1]");
		queries.put("pubchem_540","C-N-C-[#1]");
		queries.put("pubchem_541","O-C:C-C");
		queries.put("pubchem_542","O-C:C-[#1]");
		queries.put("pubchem_543","O-C:C-N");
		queries.put("pubchem_544","O-C:C-O");
		queries.put("pubchem_545","N-C:C-C");
		queries.put("pubchem_546","N-C:C-[#1]");
		queries.put("pubchem_547","N-C:C-N");
		queries.put("pubchem_548","O-C-C:C");
		queries.put("pubchem_549","N-C-C:C");
		queries.put("pubchem_550","[Cl]-C-C-C");
		queries.put("pubchem_551","[Cl]-C-C-O");
		queries.put("pubchem_552","C:C-C:C");
		queries.put("pubchem_553","O=C-C=C");
		queries.put("pubchem_554","[Br]-C-C-C");
		queries.put("pubchem_555","N=C-C=C");
		queries.put("pubchem_556","C=C-C-C");
		queries.put("pubchem_557","N:C-O-[#1]");
		queries.put("pubchem_558","O=N-C:C");
		queries.put("pubchem_559","O-C-N-[#1]");
		queries.put("pubchem_560","N-C-N-C");
		queries.put("pubchem_561","[Cl]-C-C=O");
		queries.put("pubchem_562","[Br]-C-C=O");
		queries.put("pubchem_563","O-C-O-C");
		queries.put("pubchem_564","C=C-C=C");
		queries.put("pubchem_565","C:C-O-C");
		queries.put("pubchem_566","O-C-C-N");
		queries.put("pubchem_567","O-C-C-O");
		queries.put("pubchem_568","N#C-C-C");
		queries.put("pubchem_569","N-C-C-N");
		queries.put("pubchem_570","C:C-C-C");
		queries.put("pubchem_571","[#1]-C-O-[#1]");
		queries.put("pubchem_572","N:C:N:C");
		queries.put("pubchem_573","O-C-C=C");
		queries.put("pubchem_574","O-C-C:C-C");
		queries.put("pubchem_575","O-C-C:C-O");
		queries.put("pubchem_576","N=C-C:C-[#1]");
		queries.put("pubchem_577","C:C-N-C:C");
		queries.put("pubchem_578","C-C:C-C:C");
		queries.put("pubchem_579","O=C-C-C-C");
		queries.put("pubchem_580","O=C-C-C-N");
		queries.put("pubchem_581","O=C-C-C-O");
		queries.put("pubchem_582","C-C-C-C-C");
		queries.put("pubchem_583","[Cl]-C:C-O-C");
		queries.put("pubchem_584","C:C-C=C-C");
		queries.put("pubchem_585","C-C:C-N-C");
		queries.put("pubchem_586","C-S-C-C-C");
		queries.put("pubchem_587","N-C:C-O-[#1]");
		queries.put("pubchem_588","O=C-C-C=O");
		queries.put("pubchem_589","C-C:C-O-C");
		queries.put("pubchem_590","C-C:C-O-[#1]");
		queries.put("pubchem_591","[Cl]-C-C-C-C");
		queries.put("pubchem_592","N-C-C-C-C");
		queries.put("pubchem_593","N-C-C-C-N");
		queries.put("pubchem_594","C-O-C-C=C");
		queries.put("pubchem_595","C:C-C-C-C");
		queries.put("pubchem_596","N=C-N-C-C");
		queries.put("pubchem_597","O=C-C-C:C");
		queries.put("pubchem_598","[Cl]-C:C:C-C");
		queries.put("pubchem_599","[#1]-C-C=C-[#1]");
		queries.put("pubchem_600","N-C:C:C-C");
		queries.put("pubchem_601","N-C:C:C-N");
		queries.put("pubchem_602","O=C-C-N-C");
		queries.put("pubchem_603","C-C:C:C-C");
		queries.put("pubchem_604","C-O-C-C:C");
		queries.put("pubchem_605","O=C-C-O-C");
		queries.put("pubchem_606","O-C:C-C-C");
		queries.put("pubchem_607","N-C-C-C:C");
		queries.put("pubchem_608","C-C-C-C:C");
		queries.put("pubchem_609","[Cl]-C-C-N-C");
		queries.put("pubchem_610","C-O-C-O-C");
		queries.put("pubchem_611","N-C-C-N-C");
		queries.put("pubchem_612","N-C-O-C-C");
		queries.put("pubchem_613","C-N-C-C-C");
		queries.put("pubchem_614","C-C-O-C-C");
		queries.put("pubchem_615","N-C-C-O-C");
		queries.put("pubchem_616","C:C:N:N:C");
		queries.put("pubchem_617","C-C-C-O-[#1]");
		queries.put("pubchem_618","C:C-C-C:C");
		queries.put("pubchem_619","O-C-C=C-C");
		queries.put("pubchem_620","C:C-O-C-C");
		queries.put("pubchem_621","N-C:C:C:N");
		queries.put("pubchem_622","O=C-O-C:C");
		queries.put("pubchem_623","O=C-C:C-C");
		queries.put("pubchem_624","O=C-C:C-N");
		queries.put("pubchem_625","O=C-C:C-O");
		queries.put("pubchem_626","C-O-C:C-C");
		queries.put("pubchem_627","O=[As]-C:C:C");
		queries.put("pubchem_628","C-N-C-C:C");
		queries.put("pubchem_629","S-C:C:C-N");
		queries.put("pubchem_630","O-C:C-O-C");
		queries.put("pubchem_631","O-C:C-O-[#1]");
		queries.put("pubchem_632","C-C-O-C:C");
		queries.put("pubchem_633","N-C-C:C-C");
		queries.put("pubchem_634","C-C-C:C-C");
		queries.put("pubchem_635","N-N-C-N-[#1]");
		queries.put("pubchem_636","C-N-C-N-C");
		queries.put("pubchem_637","O-C-C-C-C");
		queries.put("pubchem_638","O-C-C-C-N");
		queries.put("pubchem_639","O-C-C-C-O");
		queries.put("pubchem_640","C=C-C-C-C");
		queries.put("pubchem_641","O-C-C-C=C");
		queries.put("pubchem_642","O-C-C-C=O");
		queries.put("pubchem_643","[#1]-C-C-N-[#1]");
		queries.put("pubchem_644","C-C=N-N-C");
		queries.put("pubchem_645","O=C-N-C-C");
		queries.put("pubchem_646","O=C-N-C-[#1]");
		queries.put("pubchem_647","O=C-N-C-N");
		queries.put("pubchem_648","O=N-C:C-N");
		queries.put("pubchem_649","O=N-C:C-O");
		queries.put("pubchem_650","O=C-N-C=O");
		queries.put("pubchem_651","O-C:C:C-C");
		queries.put("pubchem_652","O-C:C:C-N");
		queries.put("pubchem_653","O-C:C:C-O");
		queries.put("pubchem_654","N-C-N-C-C");
		queries.put("pubchem_655","O-C-C-C:C");
		queries.put("pubchem_656","C-C-N-C-C");
		queries.put("pubchem_657","C-N-C:C-C");
		queries.put("pubchem_658","C-C-S-C-C");
		queries.put("pubchem_659","O-C-C-N-C");
		queries.put("pubchem_660","C-C=C-C-C");
		queries.put("pubchem_661","O-C-O-C-C");
		queries.put("pubchem_662","O-C-C-O-C");
		queries.put("pubchem_663","O-C-C-O-[#1]");
		queries.put("pubchem_664","C-C=C-C=C");
		queries.put("pubchem_665","N-C:C-C-C");
		queries.put("pubchem_666","C=C-C-O-C");
		queries.put("pubchem_667","C=C-C-O-[#1]");
		queries.put("pubchem_668","C-C:C-C-C");
		queries.put("pubchem_669","[Cl]-C:C-C=O");
		queries.put("pubchem_670","[Br]-C:C:C-C");
		queries.put("pubchem_671","O=C-C=C-C");
		queries.put("pubchem_672","O=C-C=C-[#1]");
		queries.put("pubchem_673","O=C-C=C-N");
		queries.put("pubchem_674","N-C-N-C:C");
		queries.put("pubchem_675","[Br]-C-C-C:C");
		queries.put("pubchem_676","N#C-C-C-C");
		queries.put("pubchem_677","C-C=C-C:C");
		queries.put("pubchem_678","C-C-C=C-C");
		queries.put("pubchem_679","C-C-C-C-C-C");
		queries.put("pubchem_680","O-C-C-C-C-C");
		queries.put("pubchem_681","O-C-C-C-C-O");
		queries.put("pubchem_682","O-C-C-C-C-N");
		queries.put("pubchem_683","N-C-C-C-C-C");
		queries.put("pubchem_684","O=C-C-C-C-C");
		queries.put("pubchem_685","O=C-C-C-C-N");
		queries.put("pubchem_686","O=C-C-C-C-O");
		queries.put("pubchem_687","O=C-C-C-C=O");
		queries.put("pubchem_688","C-C-C-C-C-C-C");
		queries.put("pubchem_689","O-C-C-C-C-C-C");
		queries.put("pubchem_690","O-C-C-C-C-C-O");
		queries.put("pubchem_691","O-C-C-C-C-C-N");
		queries.put("pubchem_692","O=C-C-C-C-C-C");
		queries.put("pubchem_693","O=C-C-C-C-C-O");
		queries.put("pubchem_694","O=C-C-C-C-C=O");
		queries.put("pubchem_695","O=C-C-C-C-C-N");
		queries.put("pubchem_696","C-C-C-C-C-C-C-C");
		queries.put("pubchem_697","C-C-C-C-C-C(C)-C");
		queries.put("pubchem_698","O-C-C-C-C-C-C-C");
		queries.put("pubchem_699","O-C-C-C-C-C(C)-C");
		queries.put("pubchem_700","O-C-C-C-C-C-O-C");
		queries.put("pubchem_701","O-C-C-C-C-C(O)-C");
		queries.put("pubchem_702","O-C-C-C-C-C-N-C");
		queries.put("pubchem_703","O-C-C-C-C-C(N)-C");
		queries.put("pubchem_704","O=C-C-C-C-C-C-C");
		queries.put("pubchem_705","O=C-C-C-C-C(O)-C");
		queries.put("pubchem_706","O=C-C-C-C-C(=O)-C");
		queries.put("pubchem_707","O=C-C-C-C-C(N)-C");
		queries.put("pubchem_708","C-C(C)-C-C");
		queries.put("pubchem_709","C-C(C)-C-C-C");
		queries.put("pubchem_710","C-C-C(C)-C-C");
		queries.put("pubchem_711","C-C(C)(C)-C-C");
		queries.put("pubchem_712","C-C(C)-C(C)-C");
		queries.put("pubchem_713","Cc1ccc(C)cc1");
		queries.put("pubchem_714","Cc1ccc(O)cc1");
		queries.put("pubchem_715","Cc1ccc(S)cc1");
		queries.put("pubchem_716","Cc1ccc(N)cc1");
		queries.put("pubchem_717","Cc1ccc([Cl])cc1");
		queries.put("pubchem_718","Cc1ccc([Br])cc1");
		queries.put("pubchem_719","Oc1ccc(O)cc1");
		queries.put("pubchem_720","Oc1ccc(S)cc1");
		queries.put("pubchem_721","Oc1ccc(N)cc1");
		queries.put("pubchem_722","Oc1ccc([Cl])cc1");
		queries.put("pubchem_723","Oc1ccc([Br])cc1");
		queries.put("pubchem_724","Sc1ccc(S)cc1");
		queries.put("pubchem_725","Sc1ccc(N)cc1");
		queries.put("pubchem_726","Sc1ccc([Cl])cc1");
		queries.put("pubchem_727","Sc1ccc([Br])cc1");
		queries.put("pubchem_728","Nc1ccc(N)cc1");
		queries.put("pubchem_729","Nc1ccc([Cl])cc1");
		queries.put("pubchem_730","Nc1ccc([Br])cc1");
		queries.put("pubchem_731","Clc1ccc([Cl])cc1");
		queries.put("pubchem_732","Clc1ccc([Br])cc1");
		queries.put("pubchem_733","Brc1ccc([Br])cc1");
		queries.put("pubchem_734","Cc1cc(C)ccc1");
		queries.put("pubchem_735","Cc1cc(O)ccc1");
		queries.put("pubchem_736","Cc1cc(S)ccc1");
		queries.put("pubchem_737","Cc1cc(N)ccc1");
		queries.put("pubchem_738","Cc1cc([Cl])ccc1");
		queries.put("pubchem_739","Cc1cc([Br])ccc1");
		queries.put("pubchem_740","Oc1cc(O)ccc1");
		queries.put("pubchem_741","Oc1cc(S)ccc1");
		queries.put("pubchem_742","Oc1cc(N)ccc1");
		queries.put("pubchem_743","Oc1cc([Cl])ccc1");
		queries.put("pubchem_744","Oc1cc([Br])ccc1");
		queries.put("pubchem_745","Sc1cc(S)ccc1");
		queries.put("pubchem_746","Sc1cc(N)ccc1");
		queries.put("pubchem_747","Sc1cc([Cl])ccc1");
		queries.put("pubchem_748","Sc1cc([Br])ccc1");
		queries.put("pubchem_749","Nc1cc(N)ccc1");
		queries.put("pubchem_750","Nc1cc([Cl])ccc1");
		queries.put("pubchem_751","Nc1cc([Br])ccc1");
		queries.put("pubchem_752","Clc1cc([Cl])ccc1");
		queries.put("pubchem_753","Clc1cc([Br])ccc1");
		queries.put("pubchem_754","Brc1cc([Br])ccc1");
		queries.put("pubchem_755","Cc1c(C)cccc1");
		queries.put("pubchem_756","Cc1c(O)cccc1");
		queries.put("pubchem_757","Cc1c(S)cccc1");
		queries.put("pubchem_758","Cc1c(N)cccc1");
		queries.put("pubchem_759","Cc1c([Cl])cccc1");
		queries.put("pubchem_760","Cc1c([Br])cccc1");
		queries.put("pubchem_761","Oc1c(O)cccc1");
		queries.put("pubchem_762","Oc1c(S)cccc1");
		queries.put("pubchem_763","Oc1c(N)cccc1");
		queries.put("pubchem_764","Oc1c([Cl])cccc1");
		queries.put("pubchem_765","Oc1c([Br])cccc1");
		queries.put("pubchem_766","Sc1c(S)cccc1");
		queries.put("pubchem_767","Sc1c(N)cccc1");
		queries.put("pubchem_768","Sc1c([Cl])cccc1");
		queries.put("pubchem_769","Sc1c([Br])cccc1");
		queries.put("pubchem_770","Nc1c(N)cccc1");
		queries.put("pubchem_771","Nc1c([Cl])cccc1");
		queries.put("pubchem_772","Nc1c([Br])cccc1");
		queries.put("pubchem_773","Clc1c([Cl])cccc1");
		queries.put("pubchem_774","Clc1c([Br])cccc1");
		queries.put("pubchem_775","Brc1c([Br])cccc1");
		queries.put("pubchem_776","CC1CCC(C)CC1");
		queries.put("pubchem_777","CC1CCC(O)CC1");
		queries.put("pubchem_778","CC1CCC(S)CC1");
		queries.put("pubchem_779","CC1CCC(N)CC1");
		queries.put("pubchem_780","CC1CCC([Cl])CC1");
		queries.put("pubchem_781","CC1CCC([Br])CC1");
		queries.put("pubchem_782","OC1CCC(O)CC1");
		queries.put("pubchem_783","OC1CCC(S)CC1");
		queries.put("pubchem_784","OC1CCC(N)CC1");
		queries.put("pubchem_785","OC1CCC([Cl])CC1");
		queries.put("pubchem_786","OC1CCC([Br])CC1");
		queries.put("pubchem_787","SC1CCC(S)CC1");
		queries.put("pubchem_788","SC1CCC(N)CC1");
		queries.put("pubchem_789","SC1CCC([Cl])CC1");
		queries.put("pubchem_790","SC1CCC([Br])CC1");
		queries.put("pubchem_791","NC1CCC(N)CC1");
		queries.put("pubchem_792","NC1CCC([Cl])CC1");
		queries.put("pubchem_793","NC1CCC([Br])CC1");
		queries.put("pubchem_794","ClC1CCC([Cl])CC1");
		queries.put("pubchem_795","ClC1CCC([Br])CC1");
		queries.put("pubchem_796","BrC1CCC([Br])CC1");
		queries.put("pubchem_797","CC1CC(C)CCC1");
		queries.put("pubchem_798","CC1CC(O)CCC1");
		queries.put("pubchem_799","CC1CC(S)CCC1");
		queries.put("pubchem_800","CC1CC(N)CCC1");
		queries.put("pubchem_801","CC1CC([Cl])CCC1");
		queries.put("pubchem_802","CC1CC([Br])CCC1");
		queries.put("pubchem_803","OC1CC(O)CCC1");
		queries.put("pubchem_804","OC1CC(S)CCC1");
		queries.put("pubchem_805","OC1CC(N)CCC1");
		queries.put("pubchem_806","OC1CC([Cl])CCC1");
		queries.put("pubchem_807","OC1CC([Br])CCC1");
		queries.put("pubchem_808","SC1CC(S)CCC1");
		queries.put("pubchem_809","SC1CC(N)CCC1");
		queries.put("pubchem_810","SC1CC([Cl])CCC1");
		queries.put("pubchem_811","SC1CC([Br])CCC1");
		queries.put("pubchem_812","NC1CC(N)CCC1");
		queries.put("pubchem_813","NC1CC([Cl])CCC1");
		queries.put("pubchem_814","NC1CC([Br])CCC1");
		queries.put("pubchem_815","ClC1CC([Cl])CCC1");
		queries.put("pubchem_816","ClC1CC([Br])CCC1");
		queries.put("pubchem_817","BrC1CC([Br])CCC1");
		queries.put("pubchem_818","CC1C(C)CCCC1");
		queries.put("pubchem_819","CC1C(O)CCCC1");
		queries.put("pubchem_820","CC1C(S)CCCC1");
		queries.put("pubchem_821","CC1C(N)CCCC1");
		queries.put("pubchem_822","CC1C([Cl])CCCC1");
		queries.put("pubchem_823","CC1C([Br])CCCC1");
		queries.put("pubchem_824","OC1C(O)CCCC1");
		queries.put("pubchem_825","OC1C(S)CCCC1");
		queries.put("pubchem_826","OC1C(N)CCCC1");
		queries.put("pubchem_827","OC1C([Cl])CCCC1");
		queries.put("pubchem_828","OC1C([Br])CCCC1");
		queries.put("pubchem_829","SC1C(S)CCCC1");
		queries.put("pubchem_830","SC1C(N)CCCC1");
		queries.put("pubchem_831","SC1C([Cl])CCCC1");
		queries.put("pubchem_832","SC1C([Br])CCCC1");
		queries.put("pubchem_833","NC1C(N)CCCC1");
		queries.put("pubchem_834","NC1C([Cl])CCCC1");
		queries.put("pubchem_835","NC1C([Br])CCCC1");
		queries.put("pubchem_836","ClC1C([Cl])CCCC1");
		queries.put("pubchem_837","ClC1C([Br])CCCC1");
		queries.put("pubchem_838","BrC1C([Br])CCCC1");
		queries.put("pubchem_839","CC1CC(C)CC1");
		queries.put("pubchem_840","CC1CC(O)CC1");
		queries.put("pubchem_841","CC1CC(S)CC1");
		queries.put("pubchem_842","CC1CC(N)CC1");
		queries.put("pubchem_843","CC1CC([Cl])CC1");
		queries.put("pubchem_844","CC1CC([Br])CC1");
		queries.put("pubchem_845","OC1CC(O)CC1");
		queries.put("pubchem_846","OC1CC(S)CC1");
		queries.put("pubchem_847","OC1CC(N)CC1");
		queries.put("pubchem_848","OC1CC([Cl])CC1");
		queries.put("pubchem_849","OC1CC([Br])CC1");
		queries.put("pubchem_850","SC1CC(S)CC1");
		queries.put("pubchem_851","SC1CC(N)CC1");
		queries.put("pubchem_852","SC1CC([Cl])CC1");
		queries.put("pubchem_853","SC1CC([Br])CC1");
		queries.put("pubchem_854","NC1CC(N)CC1");
		queries.put("pubchem_855","NC1CC([Cl])CC1");
		queries.put("pubchem_856","NC1CC([Br])CC1");
		queries.put("pubchem_857","ClC1CC([Cl])CC1");
		queries.put("pubchem_858","ClC1CC([Br])CC1");
		queries.put("pubchem_859","BrC1CC([Br])CC1");
		queries.put("pubchem_860","CC1C(C)CCC1");
		queries.put("pubchem_861","CC1C(O)CCC1");
		queries.put("pubchem_862","CC1C(S)CCC1");
		queries.put("pubchem_863","CC1C(N)CCC1");
		queries.put("pubchem_864","CC1C([Cl])CCC1");
		queries.put("pubchem_865","CC1C([Br])CCC1");
		queries.put("pubchem_866","OC1C(O)CCC1");
		queries.put("pubchem_867","OC1C(S)CCC1");
		queries.put("pubchem_868","OC1C(N)CCC1");
		queries.put("pubchem_869","OC1C([Cl])CCC1");
		queries.put("pubchem_870","OC1C([Br])CCC1");
		queries.put("pubchem_871","SC1C(S)CCC1");
		queries.put("pubchem_872","SC1C(N)CCC1");
		queries.put("pubchem_873","SC1C([Cl])CCC1");
		queries.put("pubchem_874","SC1C([Br])CCC1");
		queries.put("pubchem_875","NC1C(N)CCC1");
		queries.put("pubchem_876","NC1C([Cl])CC1");
		queries.put("pubchem_877","NC1C([Br])CCC1");
		queries.put("pubchem_878","ClC1C([Cl])CCC1");
		queries.put("pubchem_879","ClC1C([Br])CCC1");
		queries.put("pubchem_880","BrC1C([Br])CCC1");		
		
		
		// MACCS 166 Keys from CDK 1.5.13 (maccs.txt) / RDKit (MACCSkeys.py)
		
//		queries.put("maccs_1","?");
		queries.put("maccs_2","[#104]");
		queries.put("maccs_3","[#32,#33,#34,#50,#51,#52,#82,#83,#84]");
		queries.put("maccs_4","[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]");
		queries.put("maccs_5","[Sc,Ti,Y,Zr,Hf]");
		queries.put("maccs_6","[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]");
		queries.put("maccs_7","[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]");
		queries.put("maccs_8","[!#6;!#1]1~*~*~*~1");
		queries.put("maccs_9","[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]");
		queries.put("maccs_10","[Be,Mg,Ca,Sr,Ba,Ra]");
		queries.put("maccs_11","*1~*~*~*~1");
		queries.put("maccs_12","[Cu,Zn,Ag,Cd,Au,Hg]");
		queries.put("maccs_13","[#8]~[#7](~[#6])~[#6]");
		queries.put("maccs_14","[#16]-[#16]");
		queries.put("maccs_15","[#8]~[#6](~[#8])~[#8]");
		queries.put("maccs_16","[!#6;!#1]1~*~*~1");
		queries.put("maccs_17","[#6]#[#6]");
		queries.put("maccs_18","[#5,#13,#31,#49,#81]");
		queries.put("maccs_19","*1~*~*~*~*~*~*~1");
		queries.put("maccs_20","[#14]");
		queries.put("maccs_21","[#6]=[#6](~[!#6;!#1])~[!#6;!#1]");
		queries.put("maccs_22","*1~*~*~1");
		queries.put("maccs_23","[#7]~[#6](~[#8])~[#8]");
		queries.put("maccs_24","[#7]-[#8]");
		queries.put("maccs_25","[#7]~[#6](~[#7])~[#7]");
		queries.put("maccs_26","[#6]=;@[#6](@*)@*");
		queries.put("maccs_27","[I]");
		queries.put("maccs_28","[!#6;!#1]~[CH2]~[!#6;!#1]");
		queries.put("maccs_29","[#15]");
		queries.put("maccs_30","[#6]~[!#6;!#1](~[#6])(~[#6])~*");
		queries.put("maccs_31","[!#6;!#1]~[F,Cl,Br,I]");
		queries.put("maccs_32","[#6]~[#16]~[#7]");
		queries.put("maccs_33","[#7]~[#16]");
		queries.put("maccs_34","[CH2]=*");
		queries.put("maccs_35","[Li,Na,K,Rb,Cs,Fr]");
		queries.put("maccs_36","[#16R]");
		queries.put("maccs_37","[#7]~[#6](~[#8])~[#7]");
		queries.put("maccs_38","[#7]~[#6](~[#6])~[#7]");
		queries.put("maccs_39","[#8]~[#16](~[#8])~[#8]");
		queries.put("maccs_40","[#16]-[#8]");
		queries.put("maccs_41","[#6]#[#7]");
		queries.put("maccs_42","F");
		queries.put("maccs_43","[!#6;!#1;!H0]~*~[!#6;!#1;!H0]");
		queries.put("maccs_44","[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]");
		queries.put("maccs_45","[#6]=[#6]~[#7]");
		queries.put("maccs_46","Br");
		queries.put("maccs_47","[#16]~*~[#7]");
		queries.put("maccs_48","[#8]~[!#6;!#1](~[#8])(~[#8])");
		queries.put("maccs_49","[!+0]");
		queries.put("maccs_50","[#6]=[#6](~[#6])~[#6]");
		queries.put("maccs_51","[#6]~[#16]~[#8]");
		queries.put("maccs_52","[#7]~[#7]");
		queries.put("maccs_53","[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]");
		queries.put("maccs_54","[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]");
		queries.put("maccs_55","[#8]~[#16]~[#8]");
		queries.put("maccs_56","[#8]~[#7](~[#8])~[#6]");
		queries.put("maccs_57","[#8R]");
		queries.put("maccs_58","[!#6;!#1]~[#16]~[!#6;!#1]");
		queries.put("maccs_59","[#16]!:*:*");
		queries.put("maccs_60","[#16]=[#8]");
		queries.put("maccs_61","*~[#16](~*)~*");
		queries.put("maccs_62","*@*!@*@*");
		queries.put("maccs_63","[#7]=[#8]");
		queries.put("maccs_64","*@*!@[#16]");
		queries.put("maccs_65","c:n");
		queries.put("maccs_66","[#6]~[#6](~[#6])(~[#6])~*");
		queries.put("maccs_67","[!#6;!#1]~[#16]");
		queries.put("maccs_68","[!#6;!#1;!H0]~[!#6;!#1;!H0]");
		queries.put("maccs_69","[!#6;!#1]~[!#6;!#1;!H0]");
		queries.put("maccs_70","[!#6;!#1]~[#7]~[!#6;!#1]");
		queries.put("maccs_71","[#7]~[#8]");
		queries.put("maccs_72","[#8]~*~*~[#8]");
		queries.put("maccs_73","[#16]=*");
		queries.put("maccs_74","[CH3]~*~[CH3]");
		queries.put("maccs_75","*!@[#7]@*");
		queries.put("maccs_76","[#6]=[#6](~*)~*");
		queries.put("maccs_77","[#7]~*~[#7]");
		queries.put("maccs_78","[#6]=[#7]");
		queries.put("maccs_79","[#7]~*~*~[#7]");
		queries.put("maccs_80","[#7]~*~*~*~[#7]");
		queries.put("maccs_81","[#16]~*(~*)~*");
		queries.put("maccs_82","*~[CH2]~[!#6;!#1;!H0]");
		queries.put("maccs_83","[!#6;!#1]1~*~*~*~*~1");
		queries.put("maccs_84","[NH2]");
		queries.put("maccs_85","[#6]~[#7](~[#6])~[#6]");
		queries.put("maccs_86","[C;H2,H3][!#6;!#1][C;H2,H3]");
		queries.put("maccs_87","[F,Cl,Br,I]!@*@*");
		queries.put("maccs_88","[#16]");
		queries.put("maccs_89","[#8]~*~*~*~[#8]");
		queries.put("maccs_90","[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_91","[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_92","[#8]~[#6](~[#7])~[#6]");
		queries.put("maccs_93","[!#6;!#1]~[CH3]");
		queries.put("maccs_94","[!#6;!#1]~[#7]");
		queries.put("maccs_95","[#7]~*~*~[#8]");
		queries.put("maccs_96","*1~*~*~*~*~1");
		queries.put("maccs_97","[#7]~*~*~*~[#8]");
		queries.put("maccs_98","[!#6;!#1]1~*~*~*~*~*~1");
		queries.put("maccs_99","[#6]=[#6]");
		queries.put("maccs_100","*~[CH2]~[#7]");
		queries.put("maccs_101","[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]");
		queries.put("maccs_102","[!#6;!#1]~[#8]");
		queries.put("maccs_103","Cl");
		queries.put("maccs_104","[!#6;!#1;!H0]~*~[CH2]~*");
		queries.put("maccs_105","*@*(@*)@*");
		queries.put("maccs_106","[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]");
		queries.put("maccs_107","[F,Cl,Br,I]~*(~*)~*");
		queries.put("maccs_108","[CH3]~*~*~*~[CH2]~*");
		queries.put("maccs_109","*~[CH2]~[#8]");
		queries.put("maccs_110","[#7]~[#6]~[#8]");
		queries.put("maccs_111","[#7]~*~[CH2]~*");
		queries.put("maccs_112","*~*(~*)(~*)~*");
//		queries.put("maccs_125","?");
		queries.put("maccs_126","*!@[#8]!@*");
		queries.put("maccs_127","*@*!@[#8]");
		queries.put("maccs_128","[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_129","[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_130","[!#6;!#1]~[!#6;!#1]");
		queries.put("maccs_131","[!#6;!#1;!H0]");
		queries.put("maccs_132","[#8]~*~[CH2]~*");
		queries.put("maccs_133","*@*!@[#7]");
		queries.put("maccs_134","[F,Cl,Br,I]");
		queries.put("maccs_135","[#7]!:*:*");
		queries.put("maccs_136","[#8]=*");
		queries.put("maccs_137","[!C;!c;R]");
		queries.put("maccs_138","[!#6;!#1]~[CH2]~*");
		queries.put("maccs_139","[O;!H0]");
		queries.put("maccs_140","[#8]");
		queries.put("maccs_141","[CH3]");
		queries.put("maccs_142","[#7]");
		queries.put("maccs_143","*@*!@[#8]");
		queries.put("maccs_144","*!:*:*!:*");
		queries.put("maccs_145","*1~*~*~*~*~*~1");
		queries.put("maccs_146","[#8]");
		queries.put("maccs_147","[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]");
		queries.put("maccs_148","*~[!#6;!#1](~*)~*");
		queries.put("maccs_149","[C;H3,H4]");
		queries.put("maccs_150","*!@*@*!@*");
		queries.put("maccs_151","[#7;!H0]");
		queries.put("maccs_152","[#8]~[#6](~[#6])~[#6]");
		queries.put("maccs_153","[!#6;!#1]~[CH2]~*");
		queries.put("maccs_154","[#6]=[#8]");
		queries.put("maccs_155","*!@[CH2]!@*");
		queries.put("maccs_156","[#7]~*(~*)~*");
		queries.put("maccs_157","[#6]-[#8]");
		queries.put("maccs_158","[#6]-[#7]");
		queries.put("maccs_159","[#8]");
		queries.put("maccs_160","[C;H3,H4]");
		queries.put("maccs_161","[#7]");
		queries.put("maccs_162","a");
		queries.put("maccs_163","*1~*~*~*~*~*~1");
		queries.put("maccs_164","[#8]");
		queries.put("maccs_165","[R]");
//		queries.put("maccs_166","?");	
		
		
		
		
		
		
		
		
		
		
		
		return queries;
	}

	
	public static LinkedHashMap<String, String> getTestFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> fprint = new LinkedHashMap<String, String>();
		
		// CYP1A2
		fprint.put("ccc(c)O", "ccc(c)O");
		fprint.put("cccccO", "cccccO");	
		fprint.put("cc(c)O", "cc(c)O");	
		fprint.put("ccccO", "ccccO");	
		fprint.put("ccccccO", "ccccccO");	
		fprint.put("cccc(O)cc", "cccc(O)cc");	
		fprint.put("ccccc(c)O", "ccccc(c)O");	
		fprint.put("Oc1ccccc1", "Oc1ccccc1");	
		fprint.put("cccc(c)O", "cccc(c)O");	
		fprint.put("ccnc", "ccnc");	
		fprint.put("cnc", "cnc");	
		fprint.put("ccn", "ccn");	
		fprint.put("cccnc", "cccnc");	
		fprint.put("cccn", "cccn");	
		fprint.put("ccN", "ccN");	
		fprint.put("cccN", "cccN");	
		fprint.put("ccC", "ccC");	
		fprint.put("Cc1ccccc1", "Cc1ccccc1");	
		fprint.put("cccC", "cccC");
		fprint.put("ccc(C)cc", "ccc(C)cc");	
		fprint.put("cccc", "cccc");	
		fprint.put("ccccC", "ccccC");	
		fprint.put("ccccc", "ccccc");	
		fprint.put("cc(c)C", "cc(c)C");
		fprint.put("ccc", "ccc");	
		fprint.put("c1ccccc1", "c1ccccc1");	
		fprint.put("cccccc", "cccccc");
		fprint.put("cccc(C)cc", "cccc(C)cc");	
		fprint.put("ccc(c)C", "ccc(c)C");	
		fprint.put("cccc(c)C", "cccc(c)C");
		fprint.put("cccccC", "cccccC");	
		fprint.put("ccccc(c)C", "ccccc(c)C");	
		fprint.put("ccccccC", "ccccccC");
		fprint.put("CCN", "CCN");	
		fprint.put("CNC", "CNC");	
		fprint.put("CCC", "CCC");
		
		// CYP2A6
		fprint.put("CCNC", "CCNC");	
		fprint.put("CCCC", "CCCC");	
		fprint.put("CCO", "CCO");
		fprint.put("cccccn", "cccccn");	
		fprint.put("ccncc", "ccncc");	
		fprint.put("ccccnc", "ccccnc");
		fprint.put("ccccn", "ccccn");	
		fprint.put("c1ccncc1", "c1ccncc1");	
		fprint.put("cccncc", "cccncc");
		
		// CYP2B6
		fprint.put("CCc1ccccc1", "CCc1ccccc1");	
		fprint.put("cccc(cc)CC", "cccc(cc)CC");	
		fprint.put("ccccccCC", "ccccccCC");
		fprint.put("ccccc(c)CC", "ccccc(c)CC");	
		fprint.put("ccc(cc)CC", "ccc(cc)CC");	
		fprint.put("cccc(c)CC", "cccc(c)CC");
		fprint.put("cccccCC", "cccccCC");	
		fprint.put("ccc(c)CC", "ccc(c)CC");	
		fprint.put("ccccCC", "ccccCC");
		fprint.put("cc(c)CC", "cc(c)CC");	
		fprint.put("ccc(O)cc", "ccc(O)cc");	
		fprint.put("NC=O", "NC=O");
		fprint.put("cCN", "cCN");	
		fprint.put("ccCC", "ccCC");	
		fprint.put("ccCN", "ccCN");
		
		
		// MACCS 322
		
		
		
		
//		fprint.put("", "");	
//		fprint.put("", "");
//		fprint.put("", "");	
//		fprint.put("", "");	
//		fprint.put("", "");
//		fprint.put("", "");	
//		fprint.put("", "");
//		fprint.put("", "");	
//		fprint.put("", "");	
//		fprint.put("", "");
		
		return fprint;
	}
	
	
	public static LinkedHashMap<String, String> getPhaseIIfingerpint(){
		LinkedHashMap<String, String> phaseIIPatterns = new LinkedHashMap<String, String>();
		
		phaseIIPatterns.put("CATECHOL_O_METHYLATION_PATTERN1","[H][#8;A;X2][#6;R1]1=,:[#6;R1](-[*,#1;!$([OX2H1])])[#6;R1](-[*,#1])=,:[#6;R1](-[*,#1])[#6;R1](-[*,#1;!$([OX2H1])])=,:[#6;R1]1[#8;A;X2]C([H])([H])[H]");
		phaseIIPatterns.put("CATECHOL_O_METHYLATION_PATTERN2","[H][#8;A;X2][#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1](-[!#8])[#6;R1]([H])=,:[#6;R1]1[#8;A;X2]C([H])([H])[H]");
		phaseIIPatterns.put("CATECHOL_O_METHYLATION_PATTERN3","[H][#8]-[#6]1=,:[#6]([H])[#6]([H])=,:[#6](-[!#8])[#6]([H])=,:[#6]1-[#8]C([H])([H])[H]");
		phaseIIPatterns.put("GLUCURONIDATION","[#8]-[#6]-1-[#8]-[#6;R1](-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#6;R0](-[#8])=[O;R0]");
		phaseIIPatterns.put("GLUTATHIONE","NC(CCC(=O)NC(CS)C(=O)NCC(O)=O)C(O)=O");
		phaseIIPatterns.put("ARYLAMINE_N_ACETYLATION_PATTERN1","[H]C([H])([H])[#6](=O)-[#8][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;a]");
		phaseIIPatterns.put("ARYLAMINE_N_ACETYLATION_PATTERN2","[H][#7;A;X3]([#6;a])[#7;A;X3+0,X4+;!$([N]~[!#6])][#8]-[#6](=O)C([H])([H])[H]");
		phaseIIPatterns.put("ARALKYLAMINE_N_ACETYLATION","[#6;a][#6;A;H2X4][#6;A;H2X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])][#6](-[#6])=O");
		phaseIIPatterns.put("N_HYDROXYARYLAMINE_O_ACETYLATION","[#6]-[#6](=O)-[#8][#7;A;H1X3]!@-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
//		phaseIIPatterns.put("N_HYDROXYL_N_ARYLACETAMIDE_O_ACETYLATION","");
		phaseIIPatterns.put("SULFATE","[#6]-[#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
//		phaseIIPatterns.put("SULFATION_OF_SECONDARY_ALCOHOL","");
//		phaseIIPatterns.put("HYDRAZINE","");
		phaseIIPatterns.put("GLYCINATION_OF_ALIPHATIC_ACID","[H][#8]-[#6](=O)-[#6]-[#7]([H])-[#6;X3](-[#6;$([#7;A;R0]-,=[#6;A;R0]-,=[#6;A;R0]-,=[#6;A;R0]),$([#7]-[#6;A;X3]=[#6;A]-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)])=[O;X1]");
		phaseIIPatterns.put("GLYCINATION_OF_ARYL_ACID","[H][#8;A][#6;X3](=O)[#6;A;X4]([H])([H])[#7;X3]([H])-[#6;X3](-[#6;a;R1])=O");
//		phaseIIPatterns.put("","");
//		phaseIIPatterns.put("","");
//		phaseIIPatterns.put("","");
//		phaseIIPatterns.put("","");
//		phaseIIPatterns.put("","");
//		phaseIIPatterns.put("","");

		
		return phaseIIPatterns;
		
	}
	
	
	/**
	 * 
	 * @param mole
	 *            : A molecule of interest
	 * @param queries
	 *            : A HashMap with the functional groups and patterns to detect,
	 *            with their SMARTS patterns
	 * @return : A list containing atom-based list of the bits (1 or 0) that
	 *         describe for every functional group/pattern, whether the atom of
	 *         interest is part of a match.
	 * @throws Exception
	 */

	public ArrayList<ArrayList<Integer>> generateClassyfireAtomFingeprint(
			IAtomContainer mole, LinkedHashMap<String, String> queries) throws Exception {
		ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> allBits = new ArrayList<ArrayList<Integer>>();
		LinkedHashMap<Integer, int[]> atomFingerprints = new LinkedHashMap<Integer, int[]>();

		for (Map.Entry<String, String> item : queries.entrySet()) {
			// System.out.println(item.getKey());
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);
			List<List<Integer>> matchingAtoms = getMatchingAtoms(mole, smartsPattern);
			ArrayList<Integer> bits = new ArrayList<Integer>();

			if (matchingAtoms.size() > 0) {
				System.out.println(item.getKey());
				System.out.println("matching_atoms: " + matchingAtoms);
				for (int i = 0; i < matchingAtoms.size(); i++) {
					// System.out.println(matching_atoms.get(i));
					List<Integer> indexes = matchingAtoms.get(i);
					for (int j = 0; j < indexes.size(); j++) {
						bits.add(indexes.get(j));
						// results.add(bits);
					}
				}
			}
			if (bits.size() > 0) {
				System.out.println("Bits: " + bits);
			}

			allBits.add(new ArrayList<Integer>(new HashSet<Integer>((bits)))); // to
																				// remove
																				// duplicates
		}

		System.out.println(allBits);
		System.out.println(allBits.size());

		for (int k = 0; k < mole.getAtomCount(); k++) {
			ArrayList<Integer> nilArraylist = new ArrayList<Integer>();
			for (int a = 0; a < queries.size(); a++) {
				nilArraylist.add(0);
			}
			results.add(nilArraylist);
		}

		for (int l = 0; l < allBits.size(); l++) {
			if (allBits.get(l).size() > 0) {
				for (int m = 0; m < allBits.get(l).size(); m++) {
					int index = allBits.get(l).get(m);
					results.get(index).set(l, 1);
				}
			}
		}

		return results;
	}

	/**
	 * This function generates a unique fingerprint for a molecule.
	 * 
	 * @param molecule
	 *            : The molecule of interest
	 * @return : An array list with the fingerprint of the molecule
	 * @throws Exception
	 */

	public ArrayList[] generateClassyfireFingerprint(IAtomContainer molecule)
			throws Exception {

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
//		LinkedHashMap<String, String> queries = cs.getFingerprintPatterns();
		LinkedHashMap<String, String> queries = cs.getRINFingerprintPatterns();
		
		LinkedHashMap<String, Integer> bits = new LinkedHashMap<String, Integer>();
		int nr_of_present_groups = 0;

		ArrayList<String> fingerprint_bits = new ArrayList<String>();
		ArrayList<String> features = new ArrayList<String>();
		ArrayList[] results = new ArrayList[2];

		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		for (Map.Entry<String, String> item : queries.entrySet()) {
			// iterate over the items of a LinkedHashMap
			// System.out.println(item.getKey());
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);

			if (getMatchingAtoms(molecule, smartsPattern).size() > 0) {
				bits.put(item.getKey(), 1);
				nr_of_present_groups = nr_of_present_groups + 1;
				fingerprint_bits.add("1");
				features.add(item.getKey());
			}
			else {
				fingerprint_bits.add("0");
			}
		}
		// for (Map.Entry<String,Integer> bit : bits.entrySet()) {
		// fingerprint = fingerprint + bit.getValue().toString();
		// }
		// fingerprint = bits.get("carboxyl").toString() +
		// bits.get("hydroxyl").toString() + bits.get("aromatic").toString() +
		// bits.get("sulfuric acid").toString() +
		// bits.get("carbonyl").toString() + bits.get("aldehyde").toString() +
		// bits.get("indole").toString();
		// System.out.println("Fingerprint: " + fingerprint);
		// System.out.println("Fingerprint size: " + fingerprint.length());
		// System.out.println("Number of on bits: " + nr_of_present_groups);

		results[0] = features;
		results[1] = fingerprint_bits;
		return results;
	}


	
	
	/**
	 * This function generates a unique fingerprint for a molecule.
	 * 
	 * @param molecule
	 *            : The molecule of interest
	 * @param queries
	 *            : The dictionary of structural patterns forming the fingerprint
	 * @return : An SFingerprint representation of the molecule
	 * @throws Exception
	 */

	public SFingerprint generateClassyfireFingerprintAsDouble(IAtomContainer molecule, LinkedHashMap<String, String> queries)
			throws Exception {

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
//		LinkedHashMap<String, String> queries = cs.getFingerprintPatterns();
		int nr_of_present_groups = 0;

		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		ArrayList<String> features = new ArrayList<String>();

		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		int i = 0;
		for (Map.Entry<String, String> item : queries.entrySet()) {

			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);

			if (getMatchingAtoms(molecule, smartsPattern).size() > 0) {
				nr_of_present_groups = nr_of_present_groups + 1;
				fingerprint_bits.add(1.0);
				features.add(item.getKey());
				
			}
			else {
				fingerprint_bits.add(0.0);
			}
			i = i+1;
		}

//		System.out.println("Length -----> " + fingerprint_bits.size());
		return 	new SFingerprint(fingerprint_bits, queries);
	}


	public SFingerprint generateClassyfireRinFingerprintAsDouble(IAtomContainer molecule)
			throws Exception {

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
		LinkedHashMap<String, String> queries = cs.getRINFingerprintPatterns();
		int nr_of_present_groups = 0;

		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		ArrayList<String> features = new ArrayList<String>();

		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		Aromaticity arom = new Aromaticity(ElectronDonation.piBonds(),
                Cycles.all());
		int i = 0;
		for (Map.Entry<String, String> item : queries.entrySet()) {

//			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);

			Pattern ptrn = VentoFoggia.findSubstructure(SMARTSParser.parse(item.getValue(), builder));
			arom.apply(molecule);
//			fingerprint_bits.add(Double.valueOf(ptrn.matchAll(molecule).uniqueAtoms().count()));
			if(ptrn.matches(molecule)){
				fingerprint_bits.add(1.0);
				features.add(item.getKey());
				nr_of_present_groups = nr_of_present_groups + 1;
			} else {
				fingerprint_bits.add(0.0);
			}

			i = i+1;
		}

//		System.out.println("Length -----> " + fingerprint_bits.size());
		return 	new SFingerprint(fingerprint_bits, queries);
	}

	
	public SFingerprint generateClassyfireCountFingerprintAsDouble(IAtomContainer molecule, LinkedHashMap<String, String> queries)
			throws Exception {

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
		int nr_of_present_groups = 0;

		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		ArrayList<String> features = new ArrayList<String>();

		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		Aromaticity arom = new Aromaticity(ElectronDonation.piBonds(),
                Cycles.all());
		int i = 0;
		for (Map.Entry<String, String> item : queries.entrySet()) {

//			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);

			Pattern ptrn = VentoFoggia.findSubstructure(SMARTSParser.parse(item.getValue(), null));
			arom.apply(molecule);
			SmartsMatchers.prepare(molecule, true);
			
			int nm = ptrn.matchAll(molecule).uniqueAtoms().count();
			if(nm>0){
				fingerprint_bits.add(Double.valueOf(nm));
				features.add(item.getKey());
				nr_of_present_groups = nr_of_present_groups + 1;
			} else {
				fingerprint_bits.add(0.0);
			}

			i = i+1;
		}

//		System.out.println("Length -----> " + fingerprint_bits.size());
		return 	new SFingerprint(fingerprint_bits, queries);
	}

	
	
	public static LinkedHashMap<String,Integer> generateBTRailsCountFingerprint(IAtomContainer molecule)
			throws Exception {

		ChemStructureFingerprinter cs = new ChemStructureFingerprinter();
		int nr_of_present_groups = 0;
		LinkedHashMap<String,Integer> counts = new LinkedHashMap<String,Integer>();
		LinkedHashMap<String, String> queries = getMiniFingerprintPatterns();
		
		
		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		ArrayList<String> features = new ArrayList<String>();

		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		Aromaticity arom = new Aromaticity(ElectronDonation.piBonds(),
                Cycles.all());
		int i = 0;
		for (Map.Entry<String, String> item : queries.entrySet()) {

			Pattern ptrn = VentoFoggia.findSubstructure(SMARTSParser.parse(item.getValue(), null));
			arom.apply(molecule);
			SmartsMatchers.prepare(molecule, true);
			
			int nm = ptrn.matchAll(molecule).uniqueAtoms().count();
			counts.put(item.getKey(), nm);
			 if(nm > 0){
				 i = i+1;
			 }
		}

		System.out.println("Number of on bits" + i);
		return counts;
	}
	
	
	
	
	
	/**
	 * 
	 * @param sdfInput
	 *            : A SDF file
	 * @throws Exception
	 * 			: Throws an Exception
	 */

	public void serial_fingerprinter_sdf(File sdfInput) throws Exception {
		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		LinkedHashMap<String, Integer> bit_stats = new LinkedHashMap<String, Integer>();

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			ArrayList<String>[] results = generateClassyfireFingerprint(molecule);
			StringBuilder fingerprint_bits = new StringBuilder();
			for (String bit : results[1]) {
				fingerprint_bits.append(bit).append("\t"); // separating
															// contents using
															// commas
			}

			if (fingerprint_bits.toString().length() > 0) {
				System.out.println(i
						+ "\t"
						+ fingerprint_bits.toString().substring(0,
								(fingerprint_bits.toString().length() - 1)));
			}
			else {
				System.out.println("nil");
			}
		}
		System.out.println("\nNr of compounds: " + i);
		reader.close();
	}

	/**
	 * This function creates a tsv file with the atom fingerprints
	 * 
	 * @param sdfInput
	 *            : a SDF file
	 * @param queries
	 *            : a HashMap with fingerprint patterns and their SMARTS
	 *            expressions
	 * @throws Exception
	 * 			: Throws an Exception
	 */
	public void serialAtomfingerprinterToSdf(File sdfInput,
			LinkedHashMap<String, String> queries) throws Exception {
		String absolutePath = sdfInput.getAbsolutePath();
		String destinationDirectory = absolutePath.substring(0,
				absolutePath.lastIndexOf(File.separator));

		String[] substrings = sdfInput.getName().split("/");
		String name = substrings[substrings.length - 1].split(".sdf")[0];
		File sdfOutput = new File(destinationDirectory + "/atom_fingerprint_" + name
				+ ".tsv");
		FileOutputStream sdfOutputOs = new FileOutputStream(sdfOutput);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(sdfOutputOs));

		StringBuilder patternNames = new StringBuilder();
		for (Map.Entry<String, String> item : queries.entrySet()) {
			patternNames.append(item.getKey()).append("\t");
		}
		System.out.println(patternNames);
		bw.write("molecule ID\tatom ID\t" + patternNames + "\n"); // add header

		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		ArrayList<ArrayList<ArrayList<Integer>>> totalBits = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ElectronDonation model = ElectronDonation.cdk();
		CycleFinder cycles = Cycles.cdkAromaticSet();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule); // 1
			aromaticity.apply(molecule); // 1
			// Kekulization.kekulize(molecule); // 2

			// CDKHydrogenAdder adder =
			// CDKHydrogenAdder.getInstance(molecule.getBuilder()); // 3
			// adder.addImplicitHydrogens(molecule); // 3
			// CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
			// // 4 -
			// //http://efficientbits.blogspot.ca/2013/12/new-smiles-behaviour-parsing-cdk-154.html

			ArrayList<ArrayList<Integer>> atom_fingerprint = generateClassyfireAtomFingeprint(
					molecule, queries);
			totalBits.add(atom_fingerprint);
		}

		for (int j = 0; j < totalBits.size(); j++) {
			System.out.println("Number of atoms of molecule " + (j + 1) + ": "
					+ totalBits.get(j).size());
			for (int k = 0; k < totalBits.get(j).size(); k++) {
				ArrayList<Integer> bits = totalBits.get(j).get(k);
				// System.out.println("Number of patterns" + "\t" +
				// bits.size());
				StringBuilder fingerprintBits = new StringBuilder();

				for (int bit : bits) {
					fingerprintBits.append(bit).append("\t"); // separating
																// contents
																// using commas
				}
				bw.write((j + 1)
						+ "\t"
						+ (k + 1)
						+ "\t"
						+ fingerprintBits.deleteCharAt((fingerprintBits.length() - 1))
								.toString() + "\n");
				System.out.println((j + 1) + "\t" + (k + 1) + "\t" + fingerprintBits);
				// sdf_output.

			}
		}
		// return total_bits
		bw.close();
	}


	/**
	 * 
	 * @param sdfInput
	 *            : a SDF file
	 * @param queries
	 *            : a HashMap with fingerprint patterns and their SMARTS
	 *            expressions
	 * @return An array list of fingerprints for every atom of each molecule in
	 *         the SDF file
	 * @throws Exception
	 * 			  : Throws an Exception
	 */
	public ArrayList<ArrayList<ArrayList<Integer>>> geterateSerialAtomfingerprintToArraylist(
			File sdfInput, LinkedHashMap<String, String> queries) throws Exception {
		String absolutePath = sdfInput.getAbsolutePath();
		String destinationDirectory = absolutePath.substring(0,
				absolutePath.lastIndexOf(File.separator));

		String[] substrings = sdfInput.getName().split("/");
		String name = substrings[substrings.length - 1].split(".sdf")[0];
		File sdf_output = new File(destinationDirectory + "/atom_fingerprint_" + name
				+ ".tsv");
		FileOutputStream sdfOutputOs = new FileOutputStream(sdf_output);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(sdfOutputOs));

		StringBuilder patternNames = new StringBuilder();
		for (Map.Entry<String, String> item : queries.entrySet()) {
			patternNames.append(item.getKey()).append("\t");
		}
		System.out.println(patternNames);
		bw.write("molecule ID\tatom ID\t" + patternNames + "\n"); // add header

		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		ArrayList<ArrayList<ArrayList<Integer>>> totalBits = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ElectronDonation model = ElectronDonation.cdk();
		CycleFinder cycles = Cycles.cdkAromaticSet();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule); // 1
			aromaticity.apply(molecule); // 1
			// Kekulization.kekulize(molecule); // 2

			 CDKHydrogenAdder adder =
			 CDKHydrogenAdder.getInstance(molecule.getBuilder()); // 3
			 adder.addImplicitHydrogens(molecule); // 3
			// CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
			// // 4 -
			// //http://efficientbits.blogspot.ca/2013/12/new-smiles-behaviour-parsing-cdk-154.html
			
			System.out.println("\n\n\n\n----------------------> " + molecule.getProperty(CDKConstants.TITLE));
			ArrayList<ArrayList<Integer>> atom_fingerprint = generateClassyfireAtomFingeprint(
					molecule, queries);
			totalBits.add(atom_fingerprint);
		}

		for (int j = 0; j < totalBits.size(); j++) {
			System.out.println("Number of atoms of molecule " + (j + 1) + ": "
					+ totalBits.get(j).size());
			for (int k = 0; k < totalBits.get(j).size(); k++) {
				ArrayList<Integer> bits = totalBits.get(j).get(k);
				// System.out.println("Number of patterns" + "\t" +
				// bits.size());
				StringBuilder fingerprintBits = new StringBuilder();

				for (int bit : bits) {
					fingerprintBits.append(bit).append("\t"); // separating
																// contents
																// using commas
				}
				bw.write((j + 1)
						+ "\t"
						+ (k + 1)
						+ "\t"
						+ fingerprintBits.deleteCharAt((fingerprintBits.length() - 1))
								.toString() + "\n");
				System.out.println((j + 1) + "\t" + (k + 1) + "\t" + fingerprintBits);
				// sdf_output.

			}
		}
		bw.close();
		reader.close();
		return totalBits;

	}



	/**
	 * 
	 * @param sdfInput
	 *            : a SDF file
	 * @param queries
	 *            : a HashMap with fingerprint patterns and their SMARTS
	 *            expressions
	 * @return : A string with the atom-based fingerprints of every molecule. A
	 *         .csv file is also saved.
	 * @throws Exception
	 *  		  : Throws an Exception
	 */

	public String saveSerialAtomFingerprinterToCSV(File sdfInput,
			LinkedHashMap<String, String> queries) throws Exception {
		String absolutePath = sdfInput.getAbsolutePath();
		String destinationDirectory = absolutePath.substring(0,
				absolutePath.lastIndexOf(File.separator));
		String[] substrings = sdfInput.getName().split("/");
		String name = substrings[substrings.length - 1].split(".sdf")[0];
		String output = destinationDirectory + "/" + name + "_atom_fingerprint" + ".csv";
		File sdfOutput = new File(output);
		FileOutputStream sdfOutputOs = new FileOutputStream(sdfOutput);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(sdfOutputOs));

		StringBuilder patternNames = new StringBuilder();
		for (Map.Entry<String, String> item : queries.entrySet()) {
			patternNames.append(item.getKey()).append(",");
		}
		// System.out.println(pattern_names);
		// bw.write("molecule ID\tatom ID\t" + pattern_names + "\n"); //add
		// header

		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		ArrayList<ArrayList<ArrayList<Integer>>> totalBits = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ElectronDonation model = ElectronDonation.cdk();
		CycleFinder cycles = Cycles.cdkAromaticSet();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule); // 1
			aromaticity.apply(molecule); // 1
			// Kekulization.kekulize(molecule); // 2

			// CDKHydrogenAdder adder =
			// CDKHydrogenAdder.getInstance(molecule.getBuilder()); // 3
			// adder.addImplicitHydrogens(molecule); // 3
			// CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
			// // 4 -
			// //http://efficientbits.blogspot.ca/2013/12/new-smiles-behaviour-parsing-cdk-154.html

			ArrayList<ArrayList<Integer>> atomFingerprint = generateClassyfireAtomFingeprint(
					molecule, queries);
			totalBits.add(atomFingerprint);
		}

		for (int j = 0; j < totalBits.size(); j++) {
			// System.out.println("Number of atoms of molecule " + (j+1) + ": "
			// + total_bits.get(j).size());
			for (int k = 0; k < totalBits.get(j).size(); k++) {
				ArrayList<Integer> bits = totalBits.get(j).get(k);
				// System.out.println("Number of patterns" + "\t" +
				// bits.size());
				StringBuilder fingerprintBits = new StringBuilder();

				for (int bit : bits) {
					fingerprintBits.append(bit).append(","); // separating
																// contents
																// using commas
				}
				bw.write((j + 1)
						+ ","
						+ (k + 1)
						+ ","
						+ fingerprintBits.deleteCharAt((fingerprintBits.length() - 1))
								.toString() + "\n");
				// System.out.println( (j+1) + "," +(k+1) + "," +
				// fingerprint_bits);
				// sdf_output.

			}
		}
		bw.close();
		reader.close();
		return output;
		// return total_bits;

	}

	
	public static LinkedHashMap<String, String> getDEREP_NPFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		
		queries.put("CH3 (All)" , "[CX4H3]"); // [CH3]
		queries.put("CH3 singlet (All)" , "[CH3][*H0]");
		queries.put("CH3-Cq (singlet)" , "[CH3][CX4H0]");
		queries.put("CH3-CH (doublet)" , "[CH3][CX4H]");
		queries.put("CH3-CH2 (triplet)" , "[CH3][CX4H2]");
		queries.put("CH3 (Isopropyl)" , "[CHX4]([CH3X4])[CH3X4]");
		queries.put("CH3-Ar" , "[CH3]c");
		queries.put("CH3-Csp2" , "[CH3][CX3]");
		queries.put("CH3-C=C (Vinyl)" , "[CH3][CX3]=[CX3]");
		queries.put("CH3-CO-O (Acetyl)" , "[CX4H3][CX3](=[OX1])[OX2]");
		queries.put("CH3-CO-N (acetamide)" , "[CH3]C(=O)[#7]");
		queries.put("CH3-C-CO-" , "[CH3][#6][CX3](=[OX1])[*]");
		queries.put("CH3-O (All)" , "[CH3][OX2]");
		queries.put("CH3-O-Ar (Methoxy)" , "[CH3][O]c");
		queries.put("CH3-N (All)" , "[CH3][#7]");
		queries.put("CH3-N-CH2-" , "[CH3][$([N][CH2])][!$([CH2])]");
		queries.put("CH3-N-CH-" , "[CH3][$([N][CH])][!$([CH])]");
		queries.put("CH3-N-CO-" , "[CH3][N]([CX3]=[OX1])[!$([CX3]=[OX1])]");
		queries.put("CH2 (All)" , "[#6;H2]");
		queries.put("CH2 (sp3)" , "[CX4H2]");
		queries.put("CH2 (sp2)" , "[CX3H2]");
		queries.put("CH2-O (All)" , "[CX4H2]([!#8])[OX2]");
		queries.put("CH2-(O)O" , "[OX2][CX4H2][OX2]"); // originally CH2 (dioxy) in the paper
		queries.put("CH2-N (All)" , "[CH2][#7]");
		queries.put("CH2-N-CH2" , "[CH2][#7][CH2]");
		queries.put("CH2-N-CH" , "[CH2][#7][CH]");
		queries.put("CH2=Cq di-subst" , "[CX3H2]=[CX3H0]");  // CH2=Cq 1,1-di-subst alk
		queries.put("CH2=CH- (vinyl)" , "[CX3H2]=[CX3H1]");
		queries.put("CH (All)" , "[#6;H1]");
		queries.put("CH sp3" , "[CX4H]");
		queries.put("CH sp2 (non arom)" , "[CX3H]");
		queries.put("CH sp2 (arom)" , "[H]c");
//		queries.put("CH (arom)" , "");
//		queries.put("CH arom (singlet)" , "");
		queries.put("CH sp" , "[CX2H]");
		queries.put("CH-X (All)" , "[#6H]([#6])([#6])[#7,#8,#16]");
//		queries.put("CH-x sp3" , "");
		queries.put("CH-x (arom)" , "[H][c]([n,o,s])[c]");		
		queries.put("CH-N (All)" , "[CH]([#7])([#6])[#6]");
		queries.put("CH-O (All)" , "[#6H]([#6])([#6])[#8]");
		queries.put("CH-O (sp3)" , "[CX4H]([OX2])([#6])[#6]");
		queries.put("CH-(X)X (sp3)" , "[#6H]([#6])([#7,#8,#16])[#7,#8,#16]");
		queries.put("CH-(x)x (arom)" , "[H][c]([n,o,s])[n,o,s]");
		queries.put("CH-(O)O (All)" , "[#6H]([#8])([#8])[#6]");
		queries.put("CH-(O)N" , "[CH]([#8])([#7])[#6]");
		queries.put("CH-(O)(O)O" , "[#6H]([#8])([#8])[#8]");
		queries.put("CH (Peptide)" , "[CH]([CX3]=O)[NX3]([CX3]=O)");
		queries.put("CH (Aldehyde)" , "[CX3H](=[OX1,NX2])[#6]");
		queries.put("CH=CH cis" , "*/[CH]=[CH]\\*");
		queries.put("CH=CH trans" , "*/[CH]=[CH]/*");
//		queries.put("CH=Cq 1,1,2-tri-subst" , "");
		queries.put("Cq (All)" , "[#6;H0]");
		queries.put("Cq sp3 (All)" , "[CX4;H0]");
		queries.put("Cq sp3 (spiro)" , "[X4;R2;r4,r5,r6](@[r4,r5,r6])(@[r4,r5,r6])(@[r4,r5,r6])@[r4,r5,r6]");
		queries.put("Cq sp2 (All)" , "[CX3,c;H0]");
		queries.put("Cq sp2 (non arom)" , "[CX3H0]");
		queries.put("Cq sp2 (arom)" , "[cH0]");
		queries.put("Cq sp (alkyne)" , "[CX2;H0]");
		queries.put("Cq sp (nitrile)" , "[CX2]#[NX1]");
		queries.put("Cq sp (isonitrile)" , "[CX1-]#[NX2+]");
		queries.put("Cq sp (allene)" , "[$([CX2](=C)=C)]");
		queries.put("CO carbonyl (All)" , "[CX3]=[O]");
		queries.put("CO carbonyl (Ester/Lactone)" , "[#6][CX3](=O)[OX2H0][#6]");
		queries.put("COOH" , "[OX2H][CX3]=[OX1] ");
		queries.put("OH (All)" , "[OX2H]");
		queries.put("OH (Alcohol)" , "[CX4][OH]");
		queries.put("OH (Phenol)" , "[c]:[cX3]([OX2H]):[c]");
		queries.put("OH (Acidic)" , "[$([OH]-*=[!#6])]");
		queries.put("NH (All)" , "[#7H]");
		queries.put("NH (Arom)" , "[nH]");
		queries.put("amides" , "[#6][CX3](=O)[N;H0,H1,H2]");
		queries.put("NH2 (All)" , "[#7H2]");
		queries.put("benz 1-monosubst" , "[cX3!H]1[cX3H][cX3H][cX3H][cX3H][cX3H]1");
		queries.put("benz 1,2-disubst" , "[cX3!H]1[cX3!H][cX3H][cX3H][cX3H][cX3H]1");
		queries.put("benz 1,3-disubst" , "[cX3!H]1[cX3H][cX3!H][cX3H][cX3H][cX3H]1");
		queries.put("benz 1,4-disubst" , "[cX3!H]1[cX3H][cX3H][cX3!H][cX3H][cX3H]1");
		queries.put("benz 1,2,3-trisubst" , "[cX3!H]1[cX3!H][cX3!H][cX3H][cX3H][cX3H]1");
		queries.put("benz 1,3,5-trisubst" , "[cX3!H]1[cX3H][cX3!H][cX3H][cX3!H][cX3H]1");
		queries.put("benz 1,2,4-trisubst" , "[cX3!H]1[cX3H][cX3!H][cX3!H][cX3H][cX3H]1");
		queries.put("benz 1,2,3,4-tetrasubst" , "[cX3!H]1[cX3!H][cX3!H][cX3!H][cX3H][cX3H]1");
		queries.put("benz 1,2,3,5-tetrasubst" , "[cX3!H]1[cX3H][cX3!H][cX3!H][cX3!H][cX3H]1");
		queries.put("benz 1,2,4,5-tetrasubst" , "[cX3!H]1[cX3!H][cX3H][cX3!H][cX3!H][cX3H]1");
		queries.put("benz pentasubst" , "[cX3!H]1[cX3!H][cX3!H][cX3!H][cX3!H][cX3H]1");

		
////		queries.put("CH-X sp2" , "[CX3H][OX2,N,S] ");
//		queries.put("CH=Cq tri-subst" , "[CX3H]=[CX3!H]");
		
//		queries.put("benz hexasubst" , "[cX3!H]1[cX3!H][cX3!H][cX3!H][cX3!H][cX3!H]1");
//		queries.put("vinyl Me doublet" , "[CH3][CX3H1]=[CX3]");
//		queries.put("1,2 disub alkene All" , "[CH]=[CH]");
//		queries.put("CH3-O-CO-" , "[CH3][O][CX3](=[OX1])[#6]");
//		queries.put("imides" , "[NH]([CX3]=[O])[CX3]=[O]");
//		queries.put("NMe imide" , "[CH3][N]([CX3]=[O])[CX3]=[O]");
//		queries.put("CH3-N(CH2)CH2" , "[CH3][N]([CH2])[CH2]");
//		queries.put("CH3-N-(CH)CH" , "[CH3][N]([CH])[CH]");
//		queries.put("CH3-N-(Cq)Cq" , "[CH3][N]([CH0])[CH0]");
		
		return queries;
	}

	
	
	public LinkedHashMap<String,Integer> generateDEREP_NPFingerprint(IAtomContainer atc, boolean preprocess) throws Exception{
		LinkedHashMap<String,Integer> fingerprint = new LinkedHashMap<String,Integer>();
		LinkedHashMap<String, String> derepnp = ChemStructureFingerprinter.getDEREP_NPFingerprintPatterns();
		
		for(Map.Entry<String, String> fp : derepnp.entrySet()){
//			System.out.println(fp.getKey());
			SmartsPatternCDK scdk = new SmartsPatternCDK(fp.getValue());
			IAtomContainer atcp = atc;
			if(preprocess){
				atcp = ChemStructureManipulator.preprocessContainer(atc);			
			}
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atcp);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(atc);
//			List<List<Integer>> l = scdk.getUniqueMatchingAtoms(atcp);
			if(scdk.hasSMARTSPattern(atcp) <= 0){
				fingerprint.put(fp.getKey(),0);
			} else{
				fingerprint.put(fp.getKey(),scdk.getUniqueMatchingAtoms(atcp).size());
			}
				
			
			
		}
				
		return fingerprint;
	}
	

	public void generateDEREP_NPFingerprint(String tsvStructureFileName) throws Exception{
		BufferedReader bReader = new BufferedReader(new FileReader(tsvStructureFileName));
		System.out.println(FileUtils.getDirname(tsvStructureFileName) + "/" + FilenameUtils.getBaseName(tsvStructureFileName).toString());
		BufferedWriter bWriter = new BufferedWriter(new FileWriter( FileUtils.getDirname(tsvStructureFileName) + "/" + FilenameUtils.getBaseName(tsvStructureFileName).toString() + "_DEREP_NP.tsv"));
		
		String line = null;
		
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		String header = StringUtils.join(ChemStructureFingerprinter.getDEREP_NPFingerprintPatterns().keySet().toArray(), "\t"); 
		
//		System.out.println("Header\n" + "ID\tSMILES\t"+ header);
		bWriter.write("ID\tSMILES\t"+ header);
		bWriter.newLine();
		int errors = 0;
		
		while((line = bReader.readLine()) != null){
		
			try {
//				System.out.println(line);
				String[] sline = line.split("\t");
				String id = sline[0];
				String smiles = sline[1];
				
				IAtomContainer atc = sp.parseSmiles(smiles);
				IAtomContainer atcp = ChemStructureManipulator.preprocessContainer(atc);
				
				
				if (ChemStructureExplorer.checkConnectivity(atcp).getAtomContainerCount()==1){
					LinkedHashMap<String,Integer> c = generateDEREP_NPFingerprint(atcp, false);
					bWriter.write(id + "\t" + smiles + "\t" + StringUtils.join(c.values().toArray(),"\t"));
					bWriter.newLine();
				}else{
					bWriter.write(id + "\t" + smiles);
					bWriter.newLine();
					System.err.println("Error for: " + line);
					System.err.println("\tThe compound is a mixture");
					errors++;
					}
			} 
			catch(Exception e){
				System.err.println("Error for: " + line);
				System.err.println("\t" + e.getLocalizedMessage());
				errors++;
			}
		}
		
		System.out.println("No. of errors: " + errors);
		bWriter.close();
		bReader.close();
	
	}
	
}
