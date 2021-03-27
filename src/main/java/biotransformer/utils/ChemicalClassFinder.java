/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;

public class ChemicalClassFinder {
	
	public ChemicalClassFinder() {
		// TODO Auto-generated constructor stub
	}

	public enum ChemicalClassName {
		ALPHA_HYDROXYLATED_FATTY_ACID, BETA_HYDROXYLATED_FATTY_ACID,		
		BILE_ACID, ETHER_LIPID, 
		FATTY_ACID,
		GLYCEROLIPID, GLYCEROPHOSPHOLIPID, 
		OMEGA_HYDROXYLATED_FATTY_ACID,
		SPHINGOLIPID, 
		UNSATURATED_FATTY_ACID,
		UNFUNCTIONALIZED_UNSATURATED_FATTY_ACID,
		SATURATED_FATTY_ACID,
		UNFUNCTIONALIZED_SATURATED_FATTY_ACID,
		GLYCEROL_3_PHOSPHATE_INOSITOL,
		C23_BILE_ACID, C24_BILE_ACID,
		SULFATE_ESTER, STILBENOID, CURCUMINOID,
		TETRAPYRROLE, COENZYME_A_DERIVATIVE,SACCHARIDE,
		GLYCOSYLATED_COMPOUND, SULFATED_COMPOUND, GLUTATHIONE_CONJUGATE, 
		ACYL_CoA_CONJUGATE, GLYCINATED_COMPOUND, OLIGO_OR_POLYSACCHARIDE,
		TAURINATED_COMPOUND, GLUTAMATE_CONJUGATE, CYSTEINYLGLYCINE_S_CONJUGATE,
		ACYLCARNITINE_CONJUGATE
		
	}
	
	public static LinkedHashMap<ChemicalClassName, LinkedHashMap<String, String[]>>			chemicalClassDefinitions;
	static {
		chemicalClassDefinitions = new LinkedHashMap<ChemicalClassName, LinkedHashMap<String, String[]>>();
		
		chemicalClassDefinitions.put(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.BILE_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.BILE_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.ETHER_LIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.FATTY_ACID).put("smarts", new String[]{
				"[#6]!@-,!@=[#6]!@-,!@=[#6]-[#6;X3]([#8;A;X1-,X2H1])=[O;X1]",
				"C=C"

		});
		chemicalClassDefinitions.get(ChemicalClassName.FATTY_ACID).put("negativeSmarts", new String[]{

		});
		
		chemicalClassDefinitions.put(ChemicalClassName.SATURATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).put("smarts", new String[]{
				"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=[O;X1]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).put("negativeSmarts", new String[]{
				"C=C",
				"C#C"
		});

		chemicalClassDefinitions.put(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID).put("smarts", new String[]{
				"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=[O;X1]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID).put("negativeSmarts", new String[]{
				"C=C",
				"C#C",
				"[#6;X4][#8;A;H1X2]"
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).put("smarts", new String[]{
				"[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]"	
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).put("negativeSmarts", new String[]{
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROPHOSPHOLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).put("smarts", new String[]{
				"["
				+ "$([#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]),"
				+ "$([#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;X3R0](=O)[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).put("negativeSmarts", new String[]{

		});
		
		chemicalClassDefinitions.put(ChemicalClassName.SPHINGOLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).put("smarts", new String[]{
				"["
//				+ "$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
//				+ "$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"

				+ "$([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
//				+ "$([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"
				+ "$([H][#6]([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4])-,=[#6]([H])[#6;A;H2X4]C([H])([#1,OX2H1])[#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
				// 	5-hydroxy,3E-sphingosine (LMSP01080004)
				+ "$([H][#8]C([H])([#6;A;H2X4][#6]([H])-,=[#6]([H])[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4])[#6;H1X3]=[#6;H1X3][#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).put("negativeSmarts", new String[]{
		});		
		
		chemicalClassDefinitions.put(ChemicalClassName.ETHER_LIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("smarts", new String[]{
				"[$([#8;X2][#6;A;H2X4]!@-[#6;A;X4](!@-[!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8]!@-[#6;A;X4](!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8;X2][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("negativeSmarts", new String[]{
		});
		
			
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).put("smarts", new String[]{
				"[#8][#6;A;H1X4R1]1[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8]P([#8;X2H1,X1-])(=O)[#8]-[#6;H2X4]-[#6;H1X4](-[#6;H2X4]-[#8]-[#6]([#6,#1;A])=O)-[#8]-[#6]([#6,#1;A])=O)[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]1[#8]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).put("negativeSmarts", new String[]{
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.C24_BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).put("smarts", new String[]{
				"["
				+ "$([#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4][#6]([!#1!#6;OX1,OX2H1,NX3H1])=O)[#6]1-[#6]-,=[#6]-[#6]2-[#6]-,=3-[#6;CX4H2,$([CX4H1]-[OX2H1])]-[#6]-,=[#6]4-,=[#6]-,=[#6]([#8;A;X2H1,X1-])-[#6]-[#6]C4([#6;A;H3X4])[#6]-,=3-[#6]-[#6;CX4H2,$([CX4H1]-[OX2H1])]C12[#6;A;H3X4])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).put("negativeSmarts", new String[]{
		});		
		
		chemicalClassDefinitions.put(ChemicalClassName.C23_BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).put("smarts", new String[]{
				"["
				+ "$([#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6]([!#1!#6;OX1,OX2H1,NX3H1])=O)[#6]1-[#6]-,=[#6]-[#6]2-[#6]-,=3-[#6;CX4H2,$([CX4H1]-[OX2H1])]-[#6]-,=[#6]4-,=[#6]-,=[#6]([#8;A;X2H1,X1-])-[#6]-[#6]C4([#6;A;H3X4])[#6]-,=3-[#6]-[#6;CX4H2,$([CX4H1]-[OX2H1])]C12[#6;A;H3X4])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).put("negativeSmarts", new String[]{
		});		

		chemicalClassDefinitions.put(ChemicalClassName.SULFATE_ESTER, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).put("smarts", new String[]{
				"["
				+ "$([#6]-[#8;X2][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).put("negativeSmarts", new String[]{
		});
		chemicalClassDefinitions.put(ChemicalClassName.STILBENOID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.STILBENOID).put("smarts", new String[]{
				"["
				+ "$([#6;R0](=[#6;R0]/[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)"
				+ "]"
		});		
		chemicalClassDefinitions.get(ChemicalClassName.STILBENOID).put("negativeSmarts", new String[]{
		});		
		
		chemicalClassDefinitions.put(ChemicalClassName.CURCUMINOID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.CURCUMINOID).put("smarts", new String[]{
				"["
				+ "$([#8]~[#6](~[#6;R0]~[#6]~[#6;R1]-1=[#6;R1]-[#6](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)~[#6]~[#6](~[#8])~[#6;R0]~[#6]~[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)"
				+ ",$([#8]~[#6](~[#6;R0]~[#6]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)~[#6]~[#6](~[#8])~[#6;R0]~[#6]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)"
				+ ",$([#8]~[#6](~[#6;R0]~[#6]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)~[#6]~[#6](~[#8])~[#6;R0]~[#6]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)"
				+ ",$([#8]~[#6](~[#6]~[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)~[#6]~[#6](~[#8])~[#6;R0]~[#6]-[#6;R1]-1=[#6;R1]-[#6;R1](-[#8]-[#6])=[#6;R1](-[#8])-[#6;R1]=[#6;R1]-1)"
				+ "]"
		});		
		chemicalClassDefinitions.get(ChemicalClassName.CURCUMINOID).put("negativeSmarts", new String[]{
		});		
	
		chemicalClassDefinitions.put(ChemicalClassName.TETRAPYRROLE, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.TETRAPYRROLE).put("smarts", new String[]{
				"["
				+ "$([#6]-1-[#6]-[#7]=[#6](-[#6]-1)\\[#6]=[#6]-1\\[#6]-[#6]-[#6](=[#7]-1)-[#6]-1-[#6]-[#6]-[#6](\\[#6]=[#6]-2\\[#6]-[#6]-[#6]=[#7]-2)=[#7]-1)"
				+ ",$([#6](~c1cccn1)~c1ccc(~[#6]~c2ccc(~[#6]~c3cccn3)n2)n1)"
				+ ",$([#6](~[#6]~1~[#6]~[#6]~[#6]~[#7]~1)~[#6]~1~[#6]~[#6]~[#6](~[#6]~[#6]~2~[#6]~[#6]~[#6](~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#7]~3)~[#7]~2)~[#7]~1)"
				+ ",$([#6]-1-[#6]-[#7+]=[#6](-[#6]-1)-[#6]-1=[#6]-2-[#7]\\[#6](=[#6]/[#6]-3=[#7+]/[#6](/[#6]-[#6]3)=[#6]\\[#6]-3=[#6]-[#6]=[#6]-[#7]3)-[#6]=[#6]-2-[#6]-[#6]-1)"
				+ "]"
		});		
		chemicalClassDefinitions.get(ChemicalClassName.TETRAPYRROLE).put("negativeSmarts", new String[]{
		});
	

		chemicalClassDefinitions.put(ChemicalClassName.COENZYME_A_DERIVATIVE, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.COENZYME_A_DERIVATIVE).put("smarts", new String[]{
				"["
				+ "$([#8]-,=[#6](-[#6]-,=[#6]-[#7]-,=[#6](-,=[#8])-[#6](-[#8])-,=C([#6])([#6])-,=[#6]-,=[#8]P([#8;A;X2H1,X1-])(=O)[#8]P([#8;A;X2H1,X1-])(=O)[#8]-[#6]-[#6]-1-[#8]-[#6](-[#7]~2~[#6]~[#7]~[#6]~3~[#6]~2~[#7]~[#6]~[#7]~[#6]~3-[#7])-[#6]([#8;A;X2H1,X1-])-[#6]-1-[#8]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=O)-,=[#7]-[#6]-,=[#6]-[#16])"
				+ "]"
		});		
		chemicalClassDefinitions.get(ChemicalClassName.COENZYME_A_DERIVATIVE).put("negativeSmarts", new String[]{
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.SACCHARIDE, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SACCHARIDE).put("smarts", new String[]{
				"["
				+ "$([#6]-[#8;R0][#6;A;H1X4R1]1[#8][#6;A;H1X4R1]([#6;R0][#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]1[#8,#7,#16;A])"
				+ ",$([#6]-[#8;R0][#6;A;H1X4R1]1[#8][#6;A;H1X4R1]([#6;R0][#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]1[#8,#7,#16;A])"
//				+ "$()"
//				+ ",$()"
				+ "]"
		});		
		chemicalClassDefinitions.get(ChemicalClassName.SACCHARIDE).put("negativeSmarts", new String[]{
		});		
		
//		chemicalClassDefinitions.put(ChemicalClassName., new LinkedHashMap<String, String[]>());		
//		chemicalClassDefinitions.get(ChemicalClassName.).put("smarts", new String[]{
//				"["
//				+ "$()"
//				+ "]"
//		});		
//		chemicalClassDefinitions.get(ChemicalClassName.).put("negativeSmarts", new String[]{
//		});			

		
//		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROPHOSPHOLIPID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.OMEGA_HYDROXYLATED_FATTY_ACID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.SPHINGOLIPID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.UNSATURATED_FATTY_ACID, "");
	
	}

	public static boolean compoundMatchesClassSmartsDefinitions(IAtomContainer molecule, ChemicalClassName className) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(className).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(className).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}	
		return b;
	}
	
	public static ChemicalClassName findChemicalClass(IAtomContainer molecule) throws SMARTSException, CloneNotSupportedException, CDKException{
		ChemicalClassName chemClass = null;
//		IAtomContainer molecule = mol.clone();
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		if(isAlphaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID;
		} else if(isBetaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID;
		} else if(isOmegaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.OMEGA_HYDROXYLATED_FATTY_ACID;
		} else if(isUnsaturatedFattyAcid(molecule)){
			chemClass = ChemicalClassName.UNSATURATED_FATTY_ACID;
		} else if(isFattyAcid(molecule)){
			chemClass = ChemicalClassName.FATTY_ACID;
		}  		
		else {
			for(Map.Entry<ChemicalClassName, LinkedHashMap<String, String[]>> cc : chemicalClassDefinitions.entrySet()){
				if( 
						cc.getKey() == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL || 
						cc.getKey() == ChemicalClassName.GLYCEROLIPID || 
						cc.getKey() == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
						cc.getKey() == ChemicalClassName.ETHER_LIPID || 
						cc.getKey() == ChemicalClassName.SPHINGOLIPID ||
						cc.getKey() == ChemicalClassName.C24_BILE_ACID ||
						cc.getKey() == ChemicalClassName.C23_BILE_ACID ||
						cc.getKey() == ChemicalClassName.CURCUMINOID ||
						cc.getKey() == ChemicalClassName.STILBENOID ||
						cc.getKey() == ChemicalClassName.SACCHARIDE){
					
//					boolean b = true;
					boolean b = compoundMatchesClassSmartsDefinitions(molecule, cc.getKey());
					
//					for(String smart : chemicalClassDefinitions.get(cc.getKey()).get("smarts")){
//						SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//						if(!(pattern.hasSMARTSPattern(molecule)>0)){
//							b = false;
//							break;
//						}
//					}					
//					for(String negativeSmart : chemicalClassDefinitions.get(cc.getKey()).get("negativeSmarts")){
//						SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//						if(npattern.hasSMARTSPattern(molecule)>0){
//							b = false;
//							break;
//						}
//					}
					
					if(b == true){
						chemClass = cc.getKey();
					}					
				}
			}
		}
		

		
		return chemClass;
	}
	

	public static ArrayList<ChemicalClassName> AssignChemicalClasses(IAtomContainer molecule) throws SMARTSException, CloneNotSupportedException, CDKException{
		ArrayList<ChemicalClassName> chemClasses = new ArrayList<ChemicalClassName>();
//		IAtomContainer molecule = mol.clone();
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		if(isAlphaHydroxyFattyAcid(molecule)){
			chemClasses.add(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID) ;
		}
		if(isBetaHydroxyFattyAcid(molecule)){
			chemClasses.add(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID);
		} 
		if(isOmegaHydroxyFattyAcid(molecule)){
			chemClasses.add(ChemicalClassName.OMEGA_HYDROXYLATED_FATTY_ACID);
		}
		if(isUnsaturatedFattyAcid(molecule)){
			chemClasses.add(ChemicalClassName.UNSATURATED_FATTY_ACID);
		}
		if(isFattyAcid(molecule)){
			chemClasses.add(ChemicalClassName.FATTY_ACID);
		}			
		if(isGlycosylatedCompound(molecule)){
			chemClasses.add(ChemicalClassName.GLYCOSYLATED_COMPOUND);
		} 
		if(isSulfatedCompound(molecule)){
			chemClasses.add(ChemicalClassName.SULFATED_COMPOUND);
		} 
		if(isGlutathioneConjugate(molecule)){
			chemClasses.add(ChemicalClassName.GLUTATHIONE_CONJUGATE);
		} 
		if(isAcylCoAConjugate(molecule)){
			chemClasses.add(ChemicalClassName.ACYL_CoA_CONJUGATE);
		} 
		if(isGlycinatedCompound(molecule)){
			chemClasses.add(ChemicalClassName.GLYCINATED_COMPOUND);
		} 
		if(isTaurinatedCompound(molecule)){
			chemClasses.add(ChemicalClassName.TAURINATED_COMPOUND);
		} 		
		if(isGlutamateConjugate(molecule)){
			chemClasses.add(ChemicalClassName.GLUTAMATE_CONJUGATE);
		} 		
		if(isCysteinylGlycineConjugate(molecule)){
			chemClasses.add(ChemicalClassName.CYSTEINYLGLYCINE_S_CONJUGATE);
		} 			
		if(isAcylcarnitineConjugate(molecule)){
			chemClasses.add(ChemicalClassName.ACYLCARNITINE_CONJUGATE);
		} 			
		
		

		for(Map.Entry<ChemicalClassName, LinkedHashMap<String, String[]>> cc : chemicalClassDefinitions.entrySet()){
			if( 
					cc.getKey() == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL || 
					cc.getKey() == ChemicalClassName.GLYCEROLIPID || 
					cc.getKey() == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
					cc.getKey() == ChemicalClassName.ETHER_LIPID || 
					cc.getKey() == ChemicalClassName.SPHINGOLIPID ||
					cc.getKey() == ChemicalClassName.C24_BILE_ACID ||
					cc.getKey() == ChemicalClassName.C23_BILE_ACID ||
					cc.getKey() == ChemicalClassName.CURCUMINOID ||
					cc.getKey() == ChemicalClassName.STILBENOID ||
					cc.getKey() == ChemicalClassName.SACCHARIDE){

				if(compoundMatchesClassSmartsDefinitions(molecule, cc.getKey())) {
					chemClasses.add(cc.getKey());
				}
					
			}
		}
			
		return chemClasses;
	}
	

	
	public static boolean isUnsubstitutedSatudatedFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	
	public static boolean isFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	
	public static boolean isUnsaturatedFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	public static boolean isSaturatedFattyAcid(IAtomContainer molecule) throws SMARTSException{	
		boolean sfa = true;

//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				sfa = false;
//				break;
//			}
//		}	
//		
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				sfa = false;
//				break;
//			}
//		}
//		
//		// As per observation, 
		
		
		
		return sfa;
	}
	
	public static boolean isAlphaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isBetaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isOmegaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isEtherLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.ETHER_LIPID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}	
		return b;
	}	

	public static boolean isGlyceroLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.GLYCEROLIPID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}
	
	public static boolean isGlycerophosphoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.GLYCEROPHOSPHOLIPID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}	
	
	public static boolean isSphingoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.SPHINGOLIPID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}	

	public static boolean isGlycerol_3_PhosphateInositol(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}	

	public static boolean isC23BileAcid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.C23_BILE_ACID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}
	
	public static boolean isC24BileAcid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.C24_BILE_ACID);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}
		
		return b;

	}	
	
	
	public static boolean isSulfateEster(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.SULFATE_ESTER);
		
//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				b = false;
//				break;
//			}
//		}					
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				b = false;
//				break;
//			}
//		}	
		return b;
	}

	public static boolean isTetrapyrrole(IAtomContainer molecule) throws SMARTSException{	
		
		return compoundMatchesClassSmartsDefinitions(molecule, ChemicalClassName.TETRAPYRROLE);
	}		
	
	public static boolean isOligoOrPolysaccharide(IAtomContainer molecule) throws SMARTSException, CDKException{	
		boolean is_oligo_or_poly_saccahride = false;
		
//		if( ChemStructureExplorer.checkConnectivity(molecule).getAtomContainerCount()==1 ) {
			
			SmartsPatternCDK patterns = new SmartsPatternCDK(chemicalClassDefinitions.get(ChemicalClassName.SACCHARIDE).get("smarts")[0]);
			SmartsPatternCDK negativePatterns = null;
			if(chemicalClassDefinitions.get(ChemicalClassName.SACCHARIDE).get("negativeSmarts").length>0) {
				negativePatterns = new SmartsPatternCDK(chemicalClassDefinitions.get(ChemicalClassName.SACCHARIDE).get("negativeSmarts")[0]);
				negativePatterns.match(molecule);
			}
					
			patterns.match(molecule);
			if(patterns.getUniqueMatchingAtoms().size() > 2) {
				if(negativePatterns == null) {
					is_oligo_or_poly_saccahride = true;
				}
				else if(negativePatterns.hasSMARTSPattern(molecule)<=0) {
					is_oligo_or_poly_saccahride = true;
				}
			}
//			if(patterns.hasSMARTSPattern(molecule)>0) {
//				is_oligo_or_poly_saccahride = true;
//			}
//		}
		
		return is_oligo_or_poly_saccahride;
		
	}		

	
	public static boolean isGlycosylatedCompound(IAtomContainer molecule) throws SMARTSException {		
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;A;X2H1,X1-][#6](=O)-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6,#7,#8,#16])-[#6;R1](-[#1,#7,#8,#16;R0])-[#6;R1](-[#1,#7,#8,#16;R0])-[#6;R1]-1-[#1,#7,#8,#16;R0]");
		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[*,#1;OX2H1,$(NC(=O)C)])-[#6](!@-[#8])-[#6]-1!@-[#8]),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[OX2H1,$(NC(=O)C)])-[#6]-[#6]-1!@-[#8])"
				+ "]"
				);
		
		glucuronidePattern.match(molecule);
		glycosylMoietyPattern.match(molecule);		
		boolean isGlycosylated = glucuronidePattern.hasSMARTSPattern(molecule)>0 || glycosylMoietyPattern.hasSMARTSPattern(molecule)>0 ;
			
		return isGlycosylated;
	}
	
	public static boolean isSulfatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");	
		return sulfatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
		
	public static boolean isGlutathioneConjugate(IAtomContainer molecule) throws SMARTSException{
		SmartsPatternCDK glutathioneConjugatePattern = new SmartsPatternCDK("[H][#7]([#6;A;H2X4][#6](-[#8])=O)-[#6](=O)[#6;A;H1X4]([#6;A;H2X4][#16][#6,#8,#16;A])[#7]([H])-[#6](=O)[#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#7])[#6](-[#8])=O");	
		return glutathioneConjugatePattern.hasSMARTSPattern(molecule)>0;	
	}
	
	public static boolean isAcylCoAConjugate(IAtomContainer molecule) throws SMARTSException{
		SmartsPatternCDK glutathioneConjugatePattern = new SmartsPatternCDK("[#6;R0]-[#6;R0](=O)-[#16]-[#6]-[#6]-[#7]-[#6](=O)-[#6]-[#6]-[#7]-[#6](=O)-[#6](-[#8])C([#6])([#6])[#6]-[#8]P([#8])(=O)[#8]P([#8])(=O)[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([#8])([#8])=O)-[#7]-1-[#6]=[#7]-[#6]-2=[#6]-1-[#7]=[#6]-[#7]=[#6]-2-[#7]");	
		return glutathioneConjugatePattern.hasSMARTSPattern(molecule)>0;	
	}
	
	public static boolean isGlycinatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK glycinatedRadicalPattern = new SmartsPatternCDK("[#6][#7;A;H1X3]!@-[#6;A;H2X4]!@-[#6;X3](!@-[#8;A;X2H1,X1-])=[O;X1]");	
		return glycinatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	
	
	public static boolean isTaurinatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK taurinatedRadicalPattern = new SmartsPatternCDK("[#7;A;H1X3]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-S(!@-[#8;A;X2H1,X1-])(!@=O)!@=O");	
		return taurinatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	
	
	public static boolean isGlutamateConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK glutamatedRadicalPattern = new SmartsPatternCDK("[#7;A;H1X3][#6;A;H1X4](!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6]([#8;A;X2H1,X1-])=O)!@-[#6]([#8;A;X2H1,X1-])=O");	
		return glutamatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	
	
	public static boolean isCysteinylGlycineConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK cysteinylglycineRadicalPattern = new SmartsPatternCDK("[#6]-[#16;X2]-[#6]-[#6;X4](-[#7;X3])-[#6](=O)[#7;A;H1X3][#6;X4]-[#6]([#8;A;X2H1,X1-])=O");	
		return cysteinylglycineRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	
	public static boolean isAcylcarnitineConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK acylcarnitineRadicalPattern = new SmartsPatternCDK("[#6]-[#6](=O)!@-[#8]!@-[#6;A;H1X4]([#6;A;H2X4][#6]([#8;A;X2H1,X1-])=O)[#6;A;H2X4][N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4]");	
		return acylcarnitineRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	
	public static void main(String[] args) throws SMARTSException, CloneNotSupportedException, CDKException{
		ChemicalClassFinder ccf = new ChemicalClassFinder();
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesParser	smiParser		= new SmilesParser(builder);
		IAtomContainer molecule = smiParser.parseSmiles("CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)OC(=O)CCCCCCCCCCCCCCC");
		System.out.println(ccf.findChemicalClass(molecule));
	}
	
}
