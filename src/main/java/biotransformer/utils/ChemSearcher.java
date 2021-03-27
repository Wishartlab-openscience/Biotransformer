/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

public class ChemSearcher {
	
//	public ChemSearcher() {
//		// TODO Auto-generated constructor stub
//	}

	
	public static LinkedHashMap<String, String> getCustomAtomFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
	
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic", "[*;a]");
		queries.put("sulfuric_acid", "[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl_aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
						"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide", "[F,Cl,Br,I]");
		queries.put("acyl_chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put("organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#8;R1][#6;R1]=,:1"); // [c;R1]1[c;R1][c;R1][o;R1][c;R1]1
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic_alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put("hydroxylamine",
				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		queries.put("phenol",
				"[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");
		queries.put("primary_carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary_carboxamide",
				"[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])]-,:[#7;H1X3;$([H1][#6;!$(C=[O,N,S])])][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("tertiary_carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl_halide", "[F,Cl,Br,I][c]");
		queries.put("cons_bond", "[#6]~[#7,#8,#16]");
		queries.put("_13_dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary_aliph_ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary_arom_ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o_quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]1=[O;X1]");
		queries.put("p_quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1](=[O;X1])[#6;R1]=,:[#6;R1]1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put("thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put("amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put("n_hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7H1]-[#7H1]-[#6;a]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1_acyloxy_ethyl_ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl_carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine_monocyclic", "[#6;R1]1=,:[#6;R1][#7;R1]=,:[#6;R1][#7;R1]=,:[#6;R1]1");
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[SX2][CX2]#[NX1]");
		queries.put("alpha_beta_unsaturated system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic_alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
		queries.put("CH_acidic",
				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put("CH_acidic strong",
				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N_aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary_aliph_amine",
				"[#6;A][#7;X3H2+0,X4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary_aliph_amine",
				"[#6;A][#7;A;X3H1+0,X4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])][#6;A]");
		queries.put("tertiary_aliph_amine",
				"[#6;A][NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]([#6;A])[#6;A]");
		queries.put("primary_arom_amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary_arom_amine",
				"[#6;a][#7;A;X3H1+0,X4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])][#6;a]");
		queries.put("tertiary_arom_amine",
				"[#6;a][NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]([#6;a])[#6;a]");
		queries.put("peroxo_group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put("thioamide_s_oxide_derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put("thioamide_ss_dioxide_derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary_alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary_alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary_alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n_nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c_nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic_acid_ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic_acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put("steroid_1",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("steroid_2",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("steroid_3", 
				"[H]C1([H])[#6]-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]-3~[#6,#8,#7,#16;A;R2](~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]-[#6,#7]=,:4[#6,#7]=,:[#6,#7][#6,#7]=,:[#6,#7][#6,#7]-3=,:4)-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6]1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene_monocyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#16;R1][#6;R1]=,:1");
		queries.put("thiophene", "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1");		
		queries.put("thiophene_s_oxide", "[#8;X1-][S+]-,:1-,:[#6]=,:[#6][#6]=,:[#6]-,:1");
		queries.put("13_thiazole", "[#6]1=,:[#6][#16][#6]=,:[#7]1"); // SMARTCyp - a 2D-method for Prediction of Cytochrome P450 Mediated Drug Metabolism_Supplement
		queries.put("12_thiazole", "[#6]=,:1[#6]=,:[#7][#16][#6]=,:1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "[#6]1=,:[#6][#6]=,:[#7][#7]=,:[#6]1");
		queries.put("13_dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put("biphenyl",
				"[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("_12_aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur_compound", "[#6]~[#16]");
		queries.put("secondary_aliphatic_or_aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("propargyl_type_13_dipolar_organic_compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put("diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put("phenol_ether",
				"[#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("heteroaromatic", "[a;!c]");
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put("hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put("hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic_ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put("phosphate_monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3_acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n_acyl_aromatic_alpha_amino_acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("n_acyl_aliphatic_alpha_amino_acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl_alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[$([#6]-[#7;X3](-[#1,#6])-[#7;X3+](-[#8;X1-])=[O;X1]),$([#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3])]");
		queries.put("12_oxazole", "[#8]-1-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("13_oxazole", "[#8]-1-[#6]=[#6]-[#7]=[#6]-1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("glycerol_3_phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");
		queries.put("coenzyme_a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty_acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put("2N_linked_ribose_deriv",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic_alcohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("alpha_amino_acid_or_derivative", 
				"[#7;A][#6;X4]-[#6;X3]([!#1!#6])=[O;X1]");  //modified from ClasyFire's alpha-amino-acid-erivative-1 by adding X4 on the C2 carbon.
		queries.put("alpha_amino_acid", 
				"[#7;A;X3,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X4]-[#6;X3]([#8;A;X2H1,X1-])=[O;X1]"); //modified from ClasyFire's alpha-amino-acid-1 by adding X4 on the C2 carbon.
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put("tryptamine", 
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6;R1]1=,:[#6;R1][#7;R1][#6;R2]=,:2-[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]1=,:2");	
		queries.put("flavonol", 
				"[H][#8]-[c;R1]1[#6;R1](=O)c2-,:cc-,:cc-,:c2[#8;R1][c;R1]1-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1]([#8;A;H1X2])[c;R1]-,:1");
		queries.put("flavone", 
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:1)-[#6;R1]=,:1[#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1](=O)[#6;R1]=,:1-[!#8]"); // PMID: 15914008
		queries.put("coumarin", 
				"[$(O=[#6]-1-[#8]-c2ccccc2-[#6]=[#6]-1),$(O=[#6]-1-[#8]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]=[#6]-1)]");

		//	these should also be added to the fingerprint for SOM prediction
		queries.put("acetanilide","[#6;X4]-[#6;X3](=[O;X1])-[#7;X3]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("aniline_der","[#7;X3;H1,H2]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("o_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1-[#8;H1X2]");
		queries.put("p_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1](-[#8;H1X2])[c;R1][c;R1]1");
		queries.put("aryl_urea","[#6;a]-[#7;X3]-[#6;X3](-[#7;X3;!$([#7][!#6])])=[O;X1]");
		queries.put("aryl_carbamate","[#6;a;R1]!@-[#7]!@-[#6](=O)-[#8]-[#6;!$(C=[O,N,S])]");
		queries.put("hydroquinone","[H][#8]!@-[c;R1]1[c;R1](-[*,#1;!$([OH])])[c;R1](-[*,#1;!$([OH])])[c;R1](!@-[#8][H])[c;R1](-[*,#1;!$([OH])])[c;R1]1-[*,#1;!$([OH])]");
		queries.put("n_nitrosourea","[#7]!@-[#6;X3](!@=[O;X1])!@-[#7;X3]!@-[#7;X2]!@=O");
		queries.put("n_mustard","[$(Cl[#6]-[#6]-[#7]-[#6]-[#6]Cl),$(F[#6]-[#6]-[#7]-[#6]-[#6]F),$(Br[#6]-[#6]-[#7]-[#6]-[#6]Br),$(I[#6]-[#6]-[#7]-[#6]-[#6]I)]");

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
		queries.put("benzotrifluoride_meta_subs","[F;X1][C;X4]([F;X1])([F;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");
		queries.put("benzotrichloride_meta_subs","[Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](!@-*)[#6;R1]=,:1");			
		queries.put("benzonitrile_meta_subs","*!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1]1)!@-[C;X2]#[N;X1]");	
		queries.put("benzenesulfate_meta_subs","[#8;A;X1-,X2H1][S;X4](=[O;X1])(=[O;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	
		queries.put("anilinium_meta_subs","[#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("phenyl_1qam_meta_subs","[#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	// benzene C1-substituted with a quaternary ammoniuam salt and C3- substituted any atom.
		queries.put("nitrobenzene_meta_subs","[#8;X1-]-[#7;X3+](=[O;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	

		queries.put("edg","[$([#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1),$([#7;A;H2X3][c;R1]:[*;R1]!@-*),$([#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*),$(*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$(*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#7;X3]!@-[#6]([#6,#1;A])!@=O),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1)]");
		queries.put("ewg","[$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1),$([#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]),$([#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1)]");	
		
		
		
		// CYP isoform specificity toward drug metabolism: analysis using common feature hypothesis
		// J Mol Model (2012) 18:709–720: DOI 10.1007/s00894-011-1105-5		
		queries.put("isopropyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-*");
		queries.put("isobutyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-[#6;A;X4;H2,H3]");

		// http://www.ifm.liu.se/compchem/msi/doc/life/catalyst46/tutorials/11_excludeOrTool.doc.html
		// prim. amine OR sec. amine OR ter. amine OR amidine OR guadidino OR amidineH OR (positively charged center except dipole)
		queries.put("pos_ionizable","[$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]),$([#7;AH1;v3X3,v4X4+][#6;X3]([#7;AH1;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]),$([#7;A;H1X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6]);!$([C][N])]=[#7;A;X2;!$(NC=[O,S])]);!$([+1,+2,+3,+4,+5,+6,+7]~[+0,+1,+2,+3,+4,+5,+6,+7])]");		
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

		//  ********** From OpenBabel's FP4		
		// hits chloromethylenethers and other reactive alkylating agents
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
				
		// # 5 or 6-membered ring containing one O and at least one (r5) or two (r6) oxygen-substituents.		
		queries.put("sugar_pattern_1","[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
 
		// # 5 or 6-membered ring containing one O and an acetal-like bond at postion 2.
		queries.put("sugar_pattern_2","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
				 
		// # combination of the two above
		queries.put("sugar_pattern_combi","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");
				
		// # _5_or_6_membered_cyclic_hemi-acetal
		queries.put("sugar_pattern_2_reducing","[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
				
		// # _5_or_6_membered_cyclic_hemi_acetal
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

		queries.put("alkene_epoxidation_pattern2", "[#6]@-[#6;X3R1](-[#1,#6])=[#6;X3R1](@/[#6])-[#1,#6]");
		queries.put("aryldimethylurea", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("_3_hydroxylation_of_coumarins", "[H][#6;R1]1=,:[#6;R1][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]1=O");
		queries.put("organophosphorodithioate", "[#6][#16;A;X2][P;X4]([#8])([#8])=[S;v2X1]");
		queries.put("hydroxylation_of_methyl_carbon_adjacent_to_aliphatic_ring", "[H]C([H])([H])[#6;A;R;!$([CX4H3])]");
		queries.put("sp3_adjacent_to_n_in_nitrosamines", "[H][#6;X3]-[#7;X3;!$(N=O)]-[#7;X2]=[O;X1]");
		queries.put("aliphatic_tertiary_amine_pattern1", "[#6;A][#7;A;H0X3;+0!$([N][#6]=[!#6;!#1])!$([N]#[!#6,!#1])]([#6;A])[#6;A]");
		queries.put("hydroxylation_of_benzene_para_to_strongly_edg", "[H][#6;R1]=,:1[#6;R1](-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])])=,:[#6;R1][#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])])])=,:[#6;R1][#6;R1]=,:1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("hydroxylation_of_benzene_meta_to_ewg", "[H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[$(C([Br])([Br])[Br]),$(C([Cl])([Cl])[Cl]),$(C([F])([F])[F]),$(C([I])([I])[I]),$(C#[NX1]),$([NX1+]),$([#16;A](=O)(=O)[OX2H1,OX1-]),$([NX4H3+1]),$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])]),$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]),$([$([CX3H][#6]),$([CX3H2])]=[OX1]),$([#6][CX3](=[OX1])[#6]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])])=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1])!$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])!$([$([CX3H][#6]),$([CX3H2])]=[OX1])!$([#6][CX3](=[OX1])[#6])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])]");
		queries.put("halotrifluoroethane", "[H][C;X4]([F,Cl,Br,I])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("carbon_alpha_to_secondary_or_tertiary_alkyl_n", "[H][C;X4]([#6])([#6,#1;A])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("n_deisopropylation", "[#6;A;H3X4][#6;H1X4]([#6;A;H3X4])-[#7;X3]");
		queries.put("hydroxylation_of_aromatic_carbon_meta_to_halide_group", "[H]c:*:c-[F,Cl,Br,I]");
		queries.put("hydroxylation_of_alicyclic_secondary_carbon_pattern1", "[H][#6;A]1([#1,#6])[#6](-[*,#1])-,=[#6]1-[*,#1]");
		queries.put("formation_of_pyridinium_from_4_substituted_piperidine", "[H][#6;R1]-1-[#6;R1]([H])-[#7;X3R1](-[#6;R0])-[#6;R1]([H])-[#6;R1]([H])-[#6;R1]-1-[OX2H1,F]");
		queries.put("_2_hydroxylation_of_14_disubstituted_benzenes_pattern3", "[H][#6]1=,:[#6]([H])[#6](-[$([O][C]([F])([F])[F]),$([O][C]([H])([H])[C]([H])([H])[H]),$([N;X3+](=O)([OX1-])),$(S[C]([H])([H])[H]),$([S](=O)(=O)[N]([H])[H]),$(C(=O)[C]([H])([H])[C]([H])([H])[H]),$(C([H])([H])[H]),$(C([F])([F])[F]),$(C(=O)[OX2H1]),$(C#[N]),$(C([H])([H])C(=O)[O][H]),F,Cl,Br,I])=,:[#6]([H])[#6]([H])=,:[#6]1-[$([*;D1]),$([O][H]),$(O[C]([H])([H])[H]),$(N([H])[H]),$(N([H])[C]([H])([H])[H]),$(N([C]([H])([H])[H])[C]([H])([H])[H]),$(N([H])[C][O][C]([H])([H])[H]),$(N([H])[S](=O)(=O)[C]([H])([H])[H])]");
		queries.put("aromatic_hydroxylation_of_fused_benzene_ring_pattern3", "[#6]@[#6;R2]1=,:[#6;R2](@[#8,#7,#16])[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1]1[H]");
		queries.put("alpha_hydroxylation_of_carbonyl_group", "[H][C;X4]([#1,#6])([#1,#6])[#6;X3]([#8,#7,#6,#16,#1,#17,#9,#35,#53;A])=[O;X1]");
		queries.put("aromatic_hydroxylation_of_fused_benzene_ring_pattern4", "[H][#6;R1]=,:1[#6;R2]=,:2[#6]=,:[#6][#6][#6]3=,:[#6][#6]=,:[#6][#6;R2]([#6]=,:23)=,:[#6;R1](!@-[*,#1])[#6;R1]=,:1[H]");
		queries.put("c_next_to_amide", "[H][#6;X4R0]-[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("acyclic_aliphatic_secondary_carbon", "[H][#6;A;X4R0]([H])([#6;X4](-[*,#1])-[*,#1])[#6;X4](-[*,#1])-[*,#1]");
		queries.put("arylhydroxylamine_pattern1", "[#8;H1X2]-[#7;H1X3R0]-[#6;R1]1=,:[#6][#6]=,:[#6][#8;R1]1");
		queries.put("sterol_7_hydroxylation_pattern1", "[#6;A;H3X4]C12[#6;R1;$([CX4]-[CX3]=O),$([CX3]=O),$([CX4H1]-C(-[CX4H3])CCC(-[#1,CX4H2,CX4H3])C(-[CX4H3])[CX4H3,$([CX4H2]-[OX2H1])])]-[#6;R1]-[#6][#6;A;H1X4]1[#6]-1[#6;A;H2X4][#6;A;X4H2,X3H1]-,=[#6]3[#6;A;H2X4][#6;$([CX4H1]-[OX2H1]),$([CX3]=O)][#6;A;H2X4][#6;A;H2X4]C3([#6;A;H3X4])[#6]-1[#6;A;H2X4][#6;A;H2X4]2");
		queries.put("hydroxylation_of_benzene_ortho_to_strongly_edg", "[H][#6;R1]1=,:[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])])])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("alkyl_methine", "[H][#6;A;R0]([#6;A;X4])([#6;A;X4])[#6;A;X4]");
		queries.put("arene_epoxidation_pattern1", "[H][#6;X3R1]=,:1[#6;X3]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;X3]=,:[#6;X3](-[*,#1;!$([OX2H1])])[#6;X3R1]=,:1[H]");
		queries.put("acyclic_tertiary_amine", "[H][#6;A;X4][#7;A;R0;NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$(N*=[#7,#8,#15,#16])]([#6])[#6]");
		queries.put("pyrrolidine", "[H][#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#7;R1]-1!@-[#6]");
		queries.put("fused_benzene_ring_pattern1", "[H][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1");
		queries.put("n_alkylnitrosamine", "[H][#6;A;X4][#7;!$(N*=O)]-[#7;X2]=[O;X1]");
		queries.put("arylhydroxylamine_pattern4", "[#8;H1X2]-[#7;H1X3]!@-[#6;R1]1=,:[#6][#6][#6][#6]=,:[#6;H1X3R1]1");
		queries.put("formation_of_imminium_ion_from_n_substituted_piperidine", "[H][#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#7;X3R1]-1!@-[#6;R0,R1]");
		queries.put("aliphatic_azaheterocycle_pattern1", "[H][#6;R1]1-[#6;R1]([H])[C;R1]([H])(!@-[#6])[#6;R1]([H])-[#6;R1]([H])[#7;A;X3R1;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]1!@-[#6]");
		queries.put("o_aryl_dealkylation_not_adjacent_to_substituted_carbon", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#8;X2R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1[H]");
		queries.put("alicyclic_tertiary_amine", "[#6][#7;A;X3R1,X3R2;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])](@-[#6;A])@-[#6;A]");
		queries.put("_4_aryl_substituted_14_dihydropyridines", "[H][#7;X3R1]1[#6;R1]=,:[#6;R1]-,:[C;R1]([H])(!@-[#6;a])-,:[#6;R1]=,:[#6;R1]1");
		queries.put("cycloguanide_formation", "[H][#7;X3R0]-[#6;X3R0](-[#7;X2R0]-[#6;X3R0](!@=[#7;X3]-[*,#1])-[#7;X3R0]([H])-[#6;X4R0])!@=[#7;X3]-[*,#1]");
		queries.put("alpha_beta_unsaturated_amide", "[#6]-[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;A;X3](=[O;X1])[#6]=[#6]");
		queries.put("alicyclic_secondary_carbon_pattern3", "[H][#6;A;X4]1([H])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("secondary_arylalkyl_amine", "[H][#6;A;X4][#7;A;H1X3;+0!$(N*=[#7,#8,#15,#16])][#6;a]");
		queries.put("secondary_arylamide", "[H][#7;R0](-[*;a;R1])-[#6;R0](-[#1,#6])=[O;R0]");
		queries.put("aliphatic_secondary_antepenultimate_carbon_pattern1", "[H][C;X4R0]([H])([H])[C;R0]([H])([H])[C;R0]([H])([H])[#6]");
		queries.put("halotrifluoroethane_pattern1", "[H][C;X4]([#17,#9])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("prim_sec_aliphatic_amines", "[H][#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("allyl_group", "[*,#1;#1,#6,Br,Cl,F,I][#6](-[*,#1;#1,#6,Br,Cl,F,I])=,:[#6;X3;R0,R1,R2]([*,#1;#1,#6,Br,Cl,F,I])[#6;A;X4;R0,R1]([H])([CX4,#1])[CX4,#1]");
		queries.put("hydrazine", "[H][#6;A;X4][#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])]");
		queries.put("_14_disubstituted_benzenes_pattern1", "[H][#6]1=,:[#6]([H])[#6](-[*;D1])=,:[#6]([H])[#6]([H])=,:[#6]1-[$([O][H]),$(O[C]([H])([H])[H]),$(N([H])[H]),$(N([H])[C]([H])([H])[H]),$(N(C([H])([H])[H])[C]([H])([H])[H]),$(N([H])[C][O][C]([H])([H])[H]),$(N([H])[S](=O)(=O)[C]([H])([H])[H])]");
		queries.put("alicyclic_secondary_carbon_pattern4", "[H][#6;A]1([H])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("mixed_tertiary_amine", "[#6;a][#7;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$(N*=[#7,#8,#15,#16])]([#6])-[#6;X4][H]");
		queries.put("phosphoramide_pattern2", "[#6;A;X4;H1,H2,H3][#7;A;X3][P;X4]([#7;A;X3])([#7;A;X3])=[X1;O,S]");
		queries.put("n_demethoxymethylation", "[H][C;X4]([H])([H])[#8;X2]C([H])([H])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])][#6;a]");
		queries.put("methyl_carbon_adjacent_to_alkyne_group", "[H][#6;A;X4]([H])([H])[#6;A]#[#6;A]");
		queries.put("alicyclic_tertiary_amines_pattern1", "[H][#6;A;X4]!@-[#7;A;X3;+0!$([N]~[!#6])!$(N*=[#7,#8,#15,#16])](@-[#6;A])@-[#6;A]");
		queries.put("aryldimethylformamidine", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])-[#7;X2][#6;A;H1X3]=[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("carbon_adjacent_to_halogen_group", "[H][#6;X4]-[F,Cl,Br,I]");
		queries.put("phosphoramide_pattern1", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]P([#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])(=[O;X1])[#8]-[#1,CX4]");
		queries.put("ns_cleavage", "*-[#7;X3]!@-[#16;X2]-*");
		queries.put("_24_thiazolidinedione", "[H][#7;X3R1]-1-[#6;R1](=O)-[#6;R1]-[#16;R1;SX2,$([SX3]=[OX1]),$([SX3+]-[OX1-])]-[#6;R1]-1=O");
		queries.put("aryldimethylformamidine", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])-[#7;X2][#6;A;H1X3]=[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("arylhydroxylamine_pattern3", "[#8;H1X2]!@-[#7;H1X3]!@-[#6;R1]1=,:[#6][#7;R1]=,:[#6;X3R1][#7;R1]1");
		queries.put("hydroxylation_of_benzene_para_to_edg", "[H][#6;R1]=,:1[#6;R1](-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])])=,:[#6;R1][#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*])])=,:[#6;R1][#6;R1]=,:1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("arylhydroxylamine_pattern2", "[#8;H1X2]-[#7;H1X3]-[#6;X3R1]=,:1[#7;R1][#6]=,:[#6][#7;R1]=,:1");
		queries.put("arene_epoxidation_pattern2", "[H][c;X3R1]1-,:[c;X3R1](-[*,#1;!$([OX2H1])])[c;X3R1](-[*,#1])-,:[c;X3R1](-[*,#1])[c;R1](-[*,#1;!$([OX2H1])])-,:[c;X3R1]1-[F,Cl,Br,I,#1]");
		queries.put("hydroxylation_of_aliphatic_tertiary_penultimate_carbon", "[H][C;X4]([H])([H])[C;X4]([H])([#6])[#6]");
		queries.put("dicarboximide", "[#6,#1;A]-,:[#6;X3](=[O;R0])[#7;X3](!@-[#6;X4][H])[#6;X3]=[O;R0]");
		queries.put("m_hydroxylation_of_monosubstituted_benzene", "[H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[$(C([Br])([Br])[Br]),$(C([Cl])([Cl])[Cl]),$(C([F])([F])[F]),$(C([I])([I])[I]),$(C#[NX1]),$([NX1+]),$([#16;A](=O)(=O)[OX2H1,OX1-]),$([NX4H3+1]),$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])]),$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]),$([$([CX3H][#6]),$([CX3H2])]=[OX1]),$([#6][CX3](=[OX1])[#6]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])])=,:[#6;R1]([H])[#6;R1]=,:1[H]");
		queries.put("organophosphosulfur_compounds_pattern2", "[#6;a]-[#8;X2R0][P;X4]([#8])([#8])=[!#1!#6;X1;R0]");
		queries.put("hydroxylation_of_alicyclic_secondary_carbon_pattern5", "[H][#6;A]1([#1,!#6])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("propargylated_sec_and_tert_amines", "[H]C#C[C;R0]([H])([H])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("aryldimethylamine_pattern2", "[#6;A;H3X4][#7;A;X3]([#6;A;H3X4])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#7;R1]=,:[#6;R1][#6;R1]=,:1");
		queries.put("aryldimethylamine_pattern1", "[#6;A;H3X4][#7;A;X3]([#6;A;H3X4])[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6,#7;R1]1");
		queries.put("aliphatic_non_terminal_carbon_adjacent_to_aromatic_ring", "[H][#6;A;!$([CX4H3])](!@-[#6;a])(!@-[#6,#1;A;R0])!@-*");
		queries.put("thiocarbamate", "[#6]-[#16;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#1,#6])-[#1,#6]");
		queries.put("reduction_of_ketone_to_alcohol", "[#6]-[#6;X3](-[#6])=[O;X1]");
		queries.put("sterol_17_alpha_hydroxylation_pattern2", "[H][#6;A;X4]1([#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]2-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1][#6,#8,#7,#16]-[#6]-4=[#6]-[#6](=O)-[#6,#8,#7,#16]-[#6,#8,#7,#16]-[#6,#8,#7,#16]-4([#1,#6;A])~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-[#6,#8,#7,#16]12[#1,#6;A])[#6;R0]([#6;A;H3X4])=O");
		queries.put("sterol_17_alpha_hydroxylation_pattern1", "[H][#8]-[#6]-1-[#6,#8,#7,#16]-[#6,#8,#7,#16]-[#6,#8,#7,#16]-2([#1,#6;A])~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-[#6,#8,#7,#16]4([#1,#6;A])[#6,#8,#7,#16;A;R2](-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16][#6;A;H1X4]4[#6;R0]([#6;A;H3X4])=O)-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1][#6]=[#6]-2-[#6,#8,#7,#16]-1");
		queries.put("sterol_21_hydroxylation", "[#6;A;H3X4R0][#6;R0][#6;A;H1X4]1[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~4~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-,=[#6,#8,#7,#16]1-,=2");
		queries.put("dithiocarbamate", "[#1,#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3](-[#1,#6])-[#1,#6]");
		queries.put("o_hydroxylation_of_monosubstituted_benzene", "[H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*]),F,Cl,Br,I])=,:[#6;R1]([H])[#6;R1]=,:1[H]");
		queries.put("aliphatic_hydroxylation_of_acyclic_12_disubstituted_alkene_pattern1", "[H][#6;A;X4]([#6;X3]=*)[#6;A;X3](-[H])=[#6;A;X3](-[H])[#6;X4]");
		queries.put("n_dealkylation_of_hydrazide", "[#1,#6][#7;X3](-[#6;X4][H])[#7;X3](-[#1,#6])[#6;X3]([#1,#6])=[O;X1]");
		queries.put("_4p_hydroxylation_of_oxazaphosphorine", "[#7;X3R0;NX3H1+0,NX3H2+0,$([NX3+0;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])])][P;X4]1(=[O;X1])[#7;NX3H1,$([NX3]-[#6])]-[#6;R1;CX4H2,$([CX4H1]-[#6])]-[#6;R1]-[#6;R1]-[#8]1");
		queries.put("benzene_ortho_to_edg", "[H][#6;R1]1=,:[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*])])[#6;R1]=,:[#6;R1][#6;R1](-[*,#1;!$([OX2H1])!$([OX2]-[CX4H3])])=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("o_alkylated_group", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#8;X2R0]-[#6;$([CX4;!$(C([OX2])[O,S,#7,#15])]),c]");
		queries.put("methyl_carbon_adjacent_to_aromatic_ring", "[H][#6;A;X4]([H])([H])c:*");
		queries.put("s_adjacent_to_aliphatic_nitrogens", "[*,#1][#7;A]([*,#1])[#16;X2][#7;A]([*,#1])[*,#1]");
		queries.put("alkene_epoxidation_pattern1", "[#6][#6;A;X3R0]([#1,#6])=[#6;A;X3R0](/[#1,#6])[#1,#6]");
		queries.put("n_alkoxy_n_aryl_n_chloroacetamides", "[H][#6;X4]-[#8;X2][C;R0]([H])([H])[#7;X3R0](-[#6;R1]=,:1[#6;R1](-[#6;X4])=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1-[#6;X4])-[#6;R0](=O)[C;R0]([H])([H])Cl");
		queries.put("terminal_desaturation_pattern1", "[H][#6;X4R0](-[#6])[C;X4]([H])([H])[H]");
		queries.put("organophosphorothioate", "[#6][#8,#16;A;X2][P;X4]([#8])([#8])=[S;v2X1]");
		queries.put("secondary_heteroalicyclic_carbon_pattern1", "[H][#6;A;X4R]([H])(@-[#6;R])@-[#7,#16]");
		queries.put("guanidine_or_aminohydroazone", "[H][#7;A;X2]=[#6;X3]([#7;A;X3]([#1,#6,#7])[#1,#6,#7])[#7;A;X3]([#1,#6,#7])[#1,#6,#7]");
		queries.put("_5_aryl_14_benzodiazepines", "[!#1!#6]=,:[#6]1[#7;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6](=[#7;R1][C;R1]1([H])[#1,#6])!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");
		queries.put("_4_alkyl_hydroxylation_of_oxoazaphosphorine", "[H][C;R1]1([H])[#6;R1]-[#6;R1]-[#8;R1][P;R1]([#7])(=[O;X1R0])[#7;X3R1]1");
		queries.put("s_oxidation_of_sulfoxide_to_sulfone", "[#6][S;X3]([#6])=[O;X1]");
		queries.put("o_deisopropylation", "[H][#6;A;X4]([#6;A;H3X4])([#6;A;H3X4])[#8;X2R0]-[#6;$([CX4;!$(C([OX2])[O,S,#7,#15])]),c]");
		queries.put("secondary_heteroalicyclic_carbon_pattern2", "[H][#6;A;X4R]([H])(@-[#6,#7,#8,#16])@-[#6;R]@-[#8,#7,#16]");
		queries.put("n_aryl_n_dimethylurea", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("hydroxylation_of_aliphatic_secondary_penultimate_carbon", "[H][C;X4]([H])([H])C([H])([H])[#6]");
		queries.put("methylenedioxy_ring_opening", "[H][C;R1]1([H])[#8]-[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2-[#8]1");
		queries.put("n_dealkylation_of_ureas", "[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3;!$([#7][!#6])][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][H]");
		queries.put("p_substituted_anilide", "[H][*;A;CX4,OX2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#7]([H])-[#6](-[#1,#6])=O");
		queries.put("oxidative_deboronation", "[#6]!@-[#5;X3]([#8;A;H1X2])[#8;A;H1X2]");
		queries.put("alicyclic_tertiary_carbon", "[#6;R]@-[#6;A;H1X4R](@-[#6;A;R])@-[#6;A;R]");
		queries.put("alkyl_formamide", "[H][#7;X3](-[#6;X4])!@-[#6;X3]([H])=[O;X1]");
		queries.put("monosubstituted_benzene_p", "[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*]),F,Cl,Br,I]-[#6;R1]=,:1[#6;H1X3R1]=,:[#6;H1X3][#6;H1X3R1]=,:[#6;H1X3R1][#6;H1X3R1]=,:1");
		queries.put("thioamide", "[#1,$([#6;!$(C=[O,N,S])])][#7;X3]([#1,$([#6;!$(C=[O,N,S])])])[$([#6;X3][#6]),$([#6;X3;H1])]=[S;X1]");
		queries.put("organosulfurous_compound", "[#6]-[#8;X2][S;X3](=[O;X1])[#8;X2]-[#6]");
		queries.put("terminal_aliphatic_trihalide", "[H][C;X4](*)([F,Cl,Br,I])[F,Cl,Br,I]");
		queries.put("vinyl_ether", "[#6]-[#8;X2]-[#6;X3](-[#1,#6])=,:[#6;X3](-[#1,#6])-[#1,#6]");
		queries.put("primary_arylamine", "[H][NX3H2+0,NX4H3+;!$([N][#6]=[!#6;!#1]);!$([N]#[!#6,!#1])]c");
		queries.put("alkylsulfonamide", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4]([#1,#6])(=[O;X1])=[O;X1]");
		queries.put("carbon_alpha_to_conjugated_carbonyl", "[#8,#7,#6,#16,#1,#17,#9,#35,#53]-,:[#6;X3](=[O;X1])[#6]([H])=,:[#6]([*,#1])[*,#1]");
		queries.put("aliphatic_penultimate_carbon_adjacent_to_aromatic_carbon", "[H][C;X4]([H])([H])C([H])([#6;a])[*,#1]");
		queries.put("n_substituted_piperidine_ring", "[H]C1([H])[#6](-[*,#1])C([H])([H])C([H])([H])[#7;A;X3]([#1,c,$([CX4;A])])C1([H])[H]");
		queries.put("aromatic_carbon_para_to_halide_group", "[H]c:*:*:c-[F,Cl,Br,I]");
		queries.put("n_alkyl_diacylurea_pattern2", "[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3](-[#6;X4][H])[#6]([#1,#6])=O");
		queries.put("n_alkyl_diacylurea_pattern1", "[#1,#6][#6](=O)[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3;!$([#7][!#6])]-[#6;X4][H]");
		queries.put("halotrifluoroethane_pattern2", "[H][C;X4]([#17,#9])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("n_in_strained_ring_system", "[#6;R]@-[#7;X3R](@-[#6;R])@-[#6;R]");
		queries.put("o_methylaryl_pattern1", "[H]*[#6](=,:*[H])-[#8;X2R0]C([H])([H])[H]");
		queries.put("acyclic_secondary_amine", "[H][#6;A;X4][#7;A;H1X3R0;+0!$(N*=[#7,#8,#15,#16])][#6]");
		queries.put("carbon_gamma_to_conjugated_carbonyl", "[H][#6;X4]-[#6](-[*,#1])=[#6]([H])-[#6;X3]([#8,#7,#6,#16,#1,#17,#9,#35,#53;A])=[O;X1]");
		queries.put("bisallyl_group", "[H][C;X4](!@-[CX4,#1])([#6;A]([*,#1;#1,#6,Br,Cl,F,I])=[#6;A](/[*,#1;#1,#6,Br,Cl,F,I])[*,#1;#1,#6,Br,Cl,F,I])[#6;A]([*,#1;#1,#6,Br,Cl,F,I])=[#6;A](/[*,#1;#1,#6,Br,Cl,F,I])[*,#1;#1,#6,Br,Cl,F,I]");
		queries.put("n_chloroethyl_group", "[H]C([H])(Cl)C([H])([H])[#7;X3](-[#1,#6])-*");
		queries.put("tertiary_carboxamide_pattern2", "[#6;!$(C=[O,N,S])][#7;X3](!@-[#6;!$(C=[O,N,S])][H])[#6]([#1,#6])!@=[O;R0]");
		queries.put("aromatic_carbon_ortho_to_halide_group", "[H]c(:*):c-[F,Cl,Br,I]");
		queries.put("fused_benzene_ring_pattern2", "[#6]@[#6;R2]1=,:[#6;R2](@[#6,#8,#7,#16])[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1]1[H]");
		queries.put("terminal_methyl", "[H]C([H])([H])[#6;X4R0;!$([CX4H3])]");
		queries.put("_4_substitutted_phenol_pattern2", "[#1,CX4H3]-[#8]!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](!@-[!$([#6]-[#7]=[#7+]);!$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]);!$([#8;A;X2]S(=O)(=O)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F);!$([#8]-[#16](=O)(=O)-[#6]1=,:[#6;R1]-,:[#6;R1]=,:[#6;R1](C([H])([H])[H])-,:[#6;R1]=,:[#6;R1]1);!$([#8]S(=O)(=O)C([H])([H])[H]),$([I,Br,Cl,F]);!$([#8X3+]([H])([H])),$([#8+]([H])[H]);!$([#8X2]-[#7X3+](=[OX1])[OX1-]);!$([#8X2]-[#15](=[O])(=[O])[OX2H1,OX1-]);!$([#8X2]-[#16](=[O])(=[O])[OX2H1,OX1-]);!$([#16X2]-[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]);!$([#16X3+]([H])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]);!$([NX3+0,NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]);!$([#8;X2]-[#6](-[#1,#6])=[#8;X1]);!$([#8]-[#6]1=,:[#6;R1]-,:[#6;R1]=,:[#6;R1]-,:[#6;R1]=,:[#6;R1]1);!$([#8;X1-]);!$([#1])])=,:[#6;R1][#6;R1]=,:1");

		
		
		/**
		 * MACCS 166 Keys from CDK 1.5.13 (maccs.txt) / RDKit (MACCSkeys.py)
		 */
		// 		
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
//		queries.put("maccs_164","[#8]");
//		queries.put("maccs_165","[R]");
//		queries.put("maccs_166","?");	
		
		
		
		/**
		 * MACCS 322 - Yannick's SMARTS implementation of MACCS 322 keys, taken from MayaChemTools
		 * http://www.mayachemtools.org/docs/scripts/html/MACCSKeysFingerprints.html
		 */
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

		
		/**
		 * PubChem keys 264 to 881 (The index xxx in the pubchem_xxx is based on a start from 0)
		 * ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
		 * The first 263 keys cannot be mapped to an atom (e.g. 261	>= 4 aromatic rings)
		 */
		
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

		return queries;
	}
	
	
	
	/**
	 * 
	 * @return : A HashMap with the SMARTS expressions for functional groups and
	 *         structural patterns
	 * @throws Exception
	 */
	
	public static LinkedHashMap<String, String> getCustomBoMFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
	
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic_bond", "*:*");
		queries.put("sulfuric_acid", "[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl_aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
						"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide_bond", "*~[F,Cl,Br,I]");
		queries.put("acyl_chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put("organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#8;R1][#6;R1]=,:1"); // [c;R1]1[c;R1][c;R1][o;R1][c;R1]1
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic_alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put("hydroxylamine",
				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		queries.put("phenol",
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");
		queries.put("primary_carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary_carboxamide",
				"[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])]-,:[#7;H1X3;$([H1][#6;!$(C=[O,N,S])])][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("tertiary_carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl_halide", "[F,Cl,Br,I][c]");
		queries.put("cons_bond", "[#6]~[#7,#8,#16]");
		queries.put("_13_dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary_aliph_ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary_arom_ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o_quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]1=[O;X1]");
		queries.put("p_quinone",
				"[O;X1]=[#6;R1]1[#6;R1]=,:[#6;R1][#6;R1](=[O;X1])[#6;R1]=,:[#6;R1]1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]");
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put("thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put("amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put("n_hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7H1]-[#7H1]-[#6;a]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1_acyloxy_ethyl_ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl_carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl_ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine_monocyclic", "[#6;R1]1=,:[#6;R1][#7;R1]=,:[#6;R1][#7;R1]=,:[#6;R1]1");
		queries.put("pyrimidine", "[#6]1=,:[#6][#7]=,:[#6][#7]=,:[#6]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[SX2][CX2]#[NX1]");
		queries.put("alpha_beta_unsaturated system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic_alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
//		queries.put("CH_acidic",
//				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
//		queries.put("CH_acidic_strong",
//				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N_aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary_aliph_amine",
				"[#6;A][#7;X3H2+0,X4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary_aliph_amine",
				"[#6;A][#7;A;X3H1+0,X4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])][#6;A]");
		queries.put("tertiary_aliph_amine",
				"[#6;A][NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]([#6;A])[#6;A]");
		queries.put("primary_arom_amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary_arom_amine",
				"[#6;a][#7;A;X3H1+0,X4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])][#6;a]");
		queries.put("tertiary_arom_amine",
				"[#6;a][NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]([#6;a])[#6;a]");
		queries.put("peroxo_group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put("thioamide_s_oxide_derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put("thioamide_ss_dioxide_derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary_alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary_alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary_alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n_nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c_nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic_acid_ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic_acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put("steroid_1",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("steroid_2",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("steroid_3", 
				"[H]C1([H])[#6]-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~[#6,#8,#7,#16]-3~[#6,#8,#7,#16;A;R2](~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]-[#6,#7]=,:4[#6,#7]=,:[#6,#7][#6,#7]=,:[#6,#7][#6,#7]-3=,:4)-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6]1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene_monocyclic", "[#6;R1]=,:1[#6;R1]=,:[#6;R1][#16;R1][#6;R1]=,:1");
		queries.put("thiophene", "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1");		
		queries.put("thiophene_s_oxide", "[#8;X1-][S+]-,:1-,:[#6]=,:[#6][#6]=,:[#6]-,:1");
		queries.put("13_thiazole", "[#6]1=,:[#6][#16][#6]=,:[#7]1"); // SMARTCyp - a 2D-method for Prediction of Cytochrome P450 Mediated Drug Metabolism_Supplement
		queries.put("12_thiazole", "[#6]=,:1[#6]=,:[#7][#16][#6]=,:1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "[#6]1=,:[#6][#6]=,:[#7][#7]=,:[#6]1");
		queries.put("13_dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put("biphenyl",
				"[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("_12_aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur_compound", "[#6]~[#16]");
		queries.put("secondary_aliphatic_or_aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("propargyl_type_13_dipolar_organic_compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put("diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put("phenol_ether",
				"[#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("heteroaromatic_bond", "[!#1!#6]:c");
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put("hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put("hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic_ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put("phosphate_monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3_acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n_acyl_aromatic_alpha_amino_acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("n_acyl_aliphatic_alpha_amino_acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl_alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[$([#6]-[#7;X3](-[#1,#6])-[#7;X3+](-[#8;X1-])=[O;X1]),$([#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3])]");
		queries.put("12_oxazole", "[#8]-1-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("13_oxazole", "[#8]-1-[#6]=[#6]-[#7]=[#6]-1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("glycerol_3_phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");
		queries.put("coenzyme_a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty_acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put("2N_linked_ribose_deriv",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic_alcohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("alpha_amino_acid_or_derivative", 
				"[#7;A][#6;X4]-[#6;X3]([!#1!#6])=[O;X1]");  //modified from ClasyFire's alpha-amino-acid-erivative-1 by adding X4 on the C2 carbon.
		queries.put("alpha_amino_acid", 
				"[#7;A;X3,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X4]-[#6;X3]([#8;A;X2H1,X1-])=[O;X1]"); //modified from ClasyFire's alpha-amino-acid-1 by adding X4 on the C2 carbon.
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put("tryptamine", 
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6;R1]1=,:[#6;R1][#7;R1][#6;R2]=,:2-[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]1=,:2");	
		queries.put("flavonol", 
				"[H][#8]-[c;R1]1[#6;R1](=O)c2-,:cc-,:cc-,:c2[#8;R1][c;R1]1-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1]([#8;A;H1X2])[c;R1]-,:1");
		queries.put("flavone", 
				"[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]([#6;R1]=,:1)-[#6;R1]=,:1[#8;R1][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1](=O)[#6;R1]=,:1-[!#8]"); // PMID: 15914008
		queries.put("coumarin", 
				"[$(O=[#6]-1-[#8]-c2ccccc2-[#6]=[#6]-1),$(O=[#6]-1-[#8]-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]=[#6]-1)]");

		queries.put("acetanilide","[#6;X4]-[#6;X3](=[O;X1])-[#7;X3]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("aniline_der","[#7;X3;H1,H2]-[c;R1]1[c;R1][c;R1][c;H1X3R1][c;R1][c;R1]1");
		queries.put("o_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1-[#8;H1X2]");
		queries.put("p_aminophenol","[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][c;R1]1[c;R1][c;R1][c;X3R1](-[#8;H1X2])[c;R1][c;R1]1");
		queries.put("aryl_urea","[#6;a]-[#7;X3]-[#6;X3](-[#7;X3;!$([#7][!#6])])=[O;X1]");
		queries.put("aryl_carbamate","[#6;a;R1]!@-[#7]!@-[#6](=O)-[#8]-[#6;!$(C=[O,N,S])]");
		queries.put("hydroquinone","[H][#8]!@-[c;R1]1[c;R1](-[*,#1;!$([OH])])[c;R1](-[*,#1;!$([OH])])[c;R1](!@-[#8][H])[c;R1](-[*,#1;!$([OH])])[c;R1]1-[*,#1;!$([OH])]");
		queries.put("n_nitrosourea","[#7]!@-[#6;X3](!@=[O;X1])!@-[#7;X3]!@-[#7;X2]!@=O");
		queries.put("n_mustard","[$(Cl[#6]-[#6]-[#7]-[#6]-[#6]Cl),$(F[#6]-[#6]-[#7]-[#6]-[#6]F),$(Br[#6]-[#6]-[#7]-[#6]-[#6]Br),$(I[#6]-[#6]-[#7]-[#6]-[#6]I)]");

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
		queries.put("benzotrifluoride_meta_subs","[F;X1][C;X4]([F;X1])([F;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");
		queries.put("benzotrichloride_meta_subs","[Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1](!@-*)[#6;R1]=,:1");			
		queries.put("benzonitrile_meta_subs","*!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1]1)!@-[C;X2]#[N;X1]");	
		queries.put("benzenesulfate_meta_subs","[#8;A;X1-,X2H1][S;X4](=[O;X1])(=[O;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	
		queries.put("anilinium_meta_subs","[#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1");	
		queries.put("phenyl_1qam_meta_subs","[#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	// benzene C1-substituted with a quaternary ammoniuam salt and C3- substituted any atom.
		queries.put("nitrobenzene_meta_subs","[#8;X1-]-[#7;X3+](=[O;X1])!@-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](!@-*)=,:[#6;R1]1");	

		queries.put("edg","[$([#8;X1-]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([c;X3R1]1[c;R1][c;R1][c;X3R1][c;R1][c;R1]1),$([#7;A;H2X3][c;R1]:[*;R1]!@-*),$([#7;A;H2X3][c;R1]:[*;R1]:[*;R1]:[*;R1]!@-*),$(*!@-[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$(*!@-[*;R1]:[*;R1]:[*;R1]:[c;R1]-[$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])])]),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;R1][c;R1][c;X3R1]1!@-*),$([#8;H1X2]-[c;X3R1]1[c;R1][c;R1][c;X3R1](!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([c,$([#6;X4;!$(C=[N,O,S])])]-[#8;X2]!@-[c;R1]1[c;R1][c;R1]c(!@-*)[c;R1][c;R1]1),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6,#1]-[#6](=[O;X1])[#7;A;X3;$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6]!@-[#6;X3](!@=[O;X1])!@-[#8;X2]!@-[c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$([#7;X3]!@-[#6]([#6,#1;A])!@=O),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1!@-*),$([#6;A;X4R0][c;R1]1[c;R1][c;R1][c;R1](!@-*)[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([c;R1][c;R1]1)-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[#6;A;H1X3R0]=[#6;X3]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1]([#6;A;H1X3R0]=[#6;X3])[c;R1][c;R1]1)]");
		queries.put("ewg","[$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[F,Cl,Br,I]),$(*!@-[c;R1]1[c;R1][c;R1][c;R1](-[F,Cl,Br,I])[c;R1][c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-[#6;A;H1X3]=[O;X1])[c;R1]1),$([#6][#6;A;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#6;!$(C=[O,N,S])]!@-[#8;X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;H1X2]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1]!@-[#6;X3](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([F;X1][C;X4]([F;X1])([F;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([Cl;X1][C;X4]([Cl;X1])([Cl;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$(*!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1]([c;R1]1)!@-[C;X2]#[N;X1]),$([#8;H1X2][S;X4](=[O;X1])(=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H3X4+]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#7;A;H0X4+;!$([N][!#6])!$([N]*~[#7,#8,#15,#16])]!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1),$([#8;X1-]-[#7;X3+](=[O;X1])!@-[c;R1]1[c;R1][c;R1][c;R1][c;R1](!@-*)[c;R1]1)]");	
		

		 // CYP isoform specificity toward drug metabolism: analysis using common feature hypothesis
		 // J Mol Model (2012) 18:709–720: DOI 10.1007/s00894-011-1105-5
	
		queries.put("isopropyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-*");
		queries.put("isobutyl","[#6;A;H3X4]!@-[#6;A;H1X4](!@-[#6;A;H3X4])!@-[#6;A;X4;H2,H3]");
		
		// http://www.ifm.liu.se/compchem/msi/doc/life/catalyst46/tutorials/11_excludeOrTool.doc.html
		// prim. amine OR sec. amine OR ter. amine OR amidine OR guadidino OR amidineH OR (positively charged center except dipole)
//		queries.put("pos_ionizable","[$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([c])[C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3H0+0,X4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]),$([#7;AH1;v3X3,v4X4+][#6;X3]([#7;AH1;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]),$([#7;A;H1X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6]);!$([C][N])]=[#7;A;X2;!$(NC=[O,S])]);!$([+1,+2,+3,+4,+5,+6,+7]~[+0,+1,+2,+3,+4,+5,+6,+7])]");		
		
		queries.put("pos_ionizable_bond","[$([NX3H2+0,NX4H3+]c),$([#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]),$([#7;AH1;v3X3,v4X4+][#6;X3]([#7;AH1;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]),$([#7;A;H1X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6]);!$([C][N])]=[#7;A;X2;!$(NC=[O,S])]);!$([+1,+2,+3,+4,+5,+6,+7]~[+0,+1,+2,+3,+4,+5,+6,+7])]");
		
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

		//  ********** From OpenBabel's FP4		
		// hits chloromethylenethers and other reactive alkylating agents
		queries.put("halogen_acetal_like","[NX3v3,SX2,OX2;!$(*C=[#7,#8,#15,#16])][CX4;!$(C([N,S,O])([N,S,O])[!#6])][FX1,ClX1,BrX1,IX1]");
		
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
				
//		// # 5 or 6-membered ring containing one O and at least one (r5) or two (r6) oxygen-substituents.		
//		queries.put("sugar_pattern_1","[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]");
// 
//		// # 5 or 6-membered ring containing one O and an acetal-like bond at postion 2.
//		queries.put("sugar_pattern_2","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
//				 
//		// # combination of the two above
//		queries.put("sugar_pattern_combi","[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]");
//				
//		// # _5_or_6_membered_cyclic_hemi-acetal
//		queries.put("sugar_pattern_2_reducing","[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]");
//				
//		// # _5_or_6_membered_cyclic_hemi_acetal
//		queries.put("sugar_pattern_2_alpha","[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
//				
//		// # 5 or 6-membered cyclic hemi-acetal
//		queries.put("sugar_pattern_2_beta","[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]");
//				
//		// # pattern1 occours more than once (in same molecule, but moieties don't have to be adjacent!)
//		queries.put("poly_sugar_1","([OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)].[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)])");
//				
//		// # pattern2 occours more than once (in same molecule, but moieties don't have to be adjacent!)
//		queries.put("poly_sugar_2","([OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)].[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)])");
				
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
		
//		// # 1,3 migration of H allowed. Includes keto/enol and amide/enamide. 
//		// # Aromatic rings must stay aromatic - no keto form of phenol 
//		queries.put("_13-tautomerizable","[$([#7X2,OX1,SX1]=*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=*),$([#7X3,OX2,SX2;!H0]*:n)]");
//		queries.put("_15-tautomerizable","[$([#7X2,OX1,SX1]=,:**=,:*[!H0;!$([a;!n])]),$([#7X3,OX2,SX2;!H0]*=**=*),$([#7X3,OX2,SX2;!H0]*=,:**:n)]");
		
//		// # Hits atoms with tetrahedral chirality, if chiral center is specified in the SMILES string
//		// # depictmach does not find oxonium, sulfonium, or sulfoxides!
//		queries.put("chiral_center_specified","[$([*@](~*)(~*)(*)*),$([*@H](*)(*)*),$([*@](~*)(*)*),$([*@H](~*)~*)]");
//		
//		// # Hits atoms with tetrahedral chirality, if chiral center is not specified in the SMILES string
//		// # "@?" (unspecified chirality) is not yet supported in Open Babel Version 2.0 
//		queries.put("chiral_center_unspecified","[$([*@?](~*)(~*)(*)*),$([*@?H](*)(*)*),$([*@?](~*)(*)*),$([*@?H](~*)~*)]");

		queries.put("alkene_epoxidation_pattern2", "[#6]@-[#6;X3R1](-[#1,#6])=[#6;X3R1](@/[#6])-[#1,#6]");
		queries.put("aryldimethylurea", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("_3_hydroxylation_of_coumarins", "[H][#6;R1]1=,:[#6;R1][#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1]1=O");
		queries.put("organophosphorodithioate", "[#6][#16;A;X2][P;X4]([#8])([#8])=[S;v2X1]");
		queries.put("hydroxylation_of_methyl_carbon_adjacent_to_aliphatic_ring", "[H]C([H])([H])[#6;A;R;!$([CX4H3])]");
		queries.put("sp3_adjacent_to_n_in_nitrosamines", "[H][#6;X3]-[#7;X3;!$(N=O)]-[#7;X2]=[O;X1]");
		queries.put("aliphatic_tertiary_amine_pattern1", "[#6;A][#7;A;H0X3;+0!$([N][#6]=[!#6;!#1])!$([N]#[!#6,!#1])]([#6;A])[#6;A]");
		queries.put("hydroxylation_of_benzene_para_to_strongly_edg", "[H][#6;R1]=,:1[#6;R1](-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])])=,:[#6;R1][#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])])])=,:[#6;R1][#6;R1]=,:1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("hydroxylation_of_benzene_meta_to_ewg", "[H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[$(C([Br])([Br])[Br]),$(C([Cl])([Cl])[Cl]),$(C([F])([F])[F]),$(C([I])([I])[I]),$(C#[NX1]),$([NX1+]),$([#16;A](=O)(=O)[OX2H1,OX1-]),$([NX4H3+1]),$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])]),$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]),$([$([CX3H][#6]),$([CX3H2])]=[OX1]),$([#6][CX3](=[OX1])[#6]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])])=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1])!$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])!$([$([CX3H][#6]),$([CX3H2])]=[OX1])!$([#6][CX3](=[OX1])[#6])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])])!$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])]");
		queries.put("halotrifluoroethane", "[H][C;X4]([F,Cl,Br,I])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("carbon_alpha_to_secondary_or_tertiary_alkyl_n", "[H][C;X4]([#6])([#6,#1;A])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("n_deisopropylation", "[#6;A;H3X4][#6;H1X4]([#6;A;H3X4])-[#7;X3]");
		queries.put("hydroxylation_of_aromatic_carbon_meta_to_halide_group", "[H]c:*:c-[F,Cl,Br,I]");
		queries.put("hydroxylation_of_alicyclic_secondary_carbon_pattern1", "[H][#6;A]1([#1,#6])[#6](-[*,#1])-,=[#6]1-[*,#1]");
		queries.put("formation_of_pyridinium_from_4_substituted_piperidine", "[H][#6;R1]-1-[#6;R1]([H])-[#7;X3R1](-[#6;R0])-[#6;R1]([H])-[#6;R1]([H])-[#6;R1]-1-[OX2H1,F]");
		queries.put("_2_hydroxylation_of_14_disubstituted_benzenes_pattern3", "[H][#6]1=,:[#6]([H])[#6](-[$([O][C]([F])([F])[F]),$([O][C]([H])([H])[C]([H])([H])[H]),$([N;X3+](=O)([OX1-])),$(S[C]([H])([H])[H]),$([S](=O)(=O)[N]([H])[H]),$(C(=O)[C]([H])([H])[C]([H])([H])[H]),$(C([H])([H])[H]),$(C([F])([F])[F]),$(C(=O)[OX2H1]),$(C#[N]),$(C([H])([H])C(=O)[O][H]),F,Cl,Br,I])=,:[#6]([H])[#6]([H])=,:[#6]1-[$([*;D1]),$([O][H]),$(O[C]([H])([H])[H]),$(N([H])[H]),$(N([H])[C]([H])([H])[H]),$(N([C]([H])([H])[H])[C]([H])([H])[H]),$(N([H])[C][O][C]([H])([H])[H]),$(N([H])[S](=O)(=O)[C]([H])([H])[H])]");
		queries.put("aromatic_hydroxylation_of_fused_benzene_ring_pattern3", "[#6]@[#6;R2]1=,:[#6;R2](@[#8,#7,#16])[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1]1[H]");
		queries.put("alpha_hydroxylation_of_carbonyl_group", "[H][C;X4]([#1,#6])([#1,#6])[#6;X3]([#8,#7,#6,#16,#1,#17,#9,#35,#53;A])=[O;X1]");
		queries.put("aromatic_hydroxylation_of_fused_benzene_ring_pattern4", "[H][#6;R1]=,:1[#6;R2]=,:2[#6]=,:[#6][#6][#6]3=,:[#6][#6]=,:[#6][#6;R2]([#6]=,:23)=,:[#6;R1](!@-[*,#1])[#6;R1]=,:1[H]");
		queries.put("c_next_to_amide", "[H][#6;X4R0]-[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("acyclic_aliphatic_secondary_carbon", "[H][#6;A;X4R0]([H])([#6;X4](-[*,#1])-[*,#1])[#6;X4](-[*,#1])-[*,#1]");
		queries.put("arylhydroxylamine_pattern1", "[#8;H1X2]-[#7;H1X3R0]-[#6;R1]1=,:[#6][#6]=,:[#6][#8;R1]1");
		queries.put("sterol_7_hydroxylation_pattern1", "[#6;A;H3X4]C12[#6;R1;$([CX4]-[CX3]=O),$([CX3]=O),$([CX4H1]-C(-[CX4H3])CCC(-[#1,CX4H2,CX4H3])C(-[CX4H3])[CX4H3,$([CX4H2]-[OX2H1])])]-[#6;R1]-[#6][#6;A;H1X4]1[#6]-1[#6;A;H2X4][#6;A;X4H2,X3H1]-,=[#6]3[#6;A;H2X4][#6;$([CX4H1]-[OX2H1]),$([CX3]=O)][#6;A;H2X4][#6;A;H2X4]C3([#6;A;H3X4])[#6]-1[#6;A;H2X4][#6;A;H2X4]2");
		queries.put("hydroxylation_of_benzene_ortho_to_strongly_edg", "[H][#6;R1]1=,:[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])])])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("alkyl_methine", "[H][#6;A;R0]([#6;A;X4])([#6;A;X4])[#6;A;X4]");
		queries.put("arene_epoxidation_pattern1", "[H][#6;X3R1]=,:1[#6;X3]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;X3]=,:[#6;X3](-[*,#1;!$([OX2H1])])[#6;X3R1]=,:1[H]");
		queries.put("acyclic_tertiary_amine", "[H][#6;A;X4][#7;A;R0;NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$(N*=[#7,#8,#15,#16])]([#6])[#6]");
		queries.put("pyrrolidine", "[H][#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#7;R1]-1!@-[#6]");
		queries.put("fused_benzene_ring_pattern1", "[H][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1");
		queries.put("n_alkylnitrosamine", "[H][#6;A;X4][#7;!$(N*=O)]-[#7;X2]=[O;X1]");
		queries.put("arylhydroxylamine_pattern4", "[#8;H1X2]-[#7;H1X3]!@-[#6;R1]1=,:[#6][#6][#6][#6]=,:[#6;H1X3R1]1");
		queries.put("formation_of_imminium_ion_from_n_substituted_piperidine", "[H][#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#7;X3R1]-1!@-[#6;R0,R1]");
		queries.put("aliphatic_azaheterocycle_pattern1", "[H][#6;R1]1-[#6;R1]([H])[C;R1]([H])(!@-[#6])[#6;R1]([H])-[#6;R1]([H])[#7;A;X3R1;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]1!@-[#6]");
		queries.put("o_aryl_dealkylation_not_adjacent_to_substituted_carbon", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#8;X2R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1[H]");
		queries.put("alicyclic_tertiary_amine", "[#6][#7;A;X3R1,X3R2;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])](@-[#6;A])@-[#6;A]");
		queries.put("_4_aryl_substituted_14_dihydropyridines", "[H][#7;X3R1]1[#6;R1]=,:[#6;R1]-,:[C;R1]([H])(!@-[#6;a])-,:[#6;R1]=,:[#6;R1]1");
		queries.put("cycloguanide_formation", "[H][#7;X3R0]-[#6;X3R0](-[#7;X2R0]-[#6;X3R0](!@=[#7;X3]-[*,#1])-[#7;X3R0]([H])-[#6;X4R0])!@=[#7;X3]-[*,#1]");
		queries.put("alpha_beta_unsaturated_amide", "[#6]-[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#6;A;X3](=[O;X1])[#6]=[#6]");
		queries.put("alicyclic_secondary_carbon_pattern3", "[H][#6;A;X4]1([H])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("secondary_arylalkyl_amine", "[H][#6;A;X4][#7;A;H1X3;+0!$(N*=[#7,#8,#15,#16])][#6;a]");
		queries.put("secondary_arylamide", "[H][#7;R0](-[*;a;R1])-[#6;R0](-[#1,#6])=[O;R0]");
		queries.put("aliphatic_secondary_antepenultimate_carbon_pattern1", "[H][C;X4R0]([H])([H])[C;R0]([H])([H])[C;R0]([H])([H])[#6]");
		queries.put("halotrifluoroethane_pattern1", "[H][C;X4]([#17,#9])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("prim_sec_aliphatic_amines", "[H][#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("allyl_group", "[*,#1;#1,#6,Br,Cl,F,I][#6](-[*,#1;#1,#6,Br,Cl,F,I])=,:[#6;X3;R0,R1,R2]([*,#1;#1,#6,Br,Cl,F,I])[#6;A;X4;R0,R1]([H])([CX4,#1])[CX4,#1]");
		queries.put("hydrazine", "[H][#6;A;X4][#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])]");
		queries.put("_14_disubstituted_benzenes_pattern1", "[H][#6]1=,:[#6]([H])[#6](-[*;D1])=,:[#6]([H])[#6]([H])=,:[#6]1-[$([O][H]),$(O[C]([H])([H])[H]),$(N([H])[H]),$(N([H])[C]([H])([H])[H]),$(N(C([H])([H])[H])[C]([H])([H])[H]),$(N([H])[C][O][C]([H])([H])[H]),$(N([H])[S](=O)(=O)[C]([H])([H])[H])]");
		queries.put("alicyclic_secondary_carbon_pattern4", "[H][#6;A]1([H])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("mixed_tertiary_amine", "[#6;a][#7;NX3H0+0,NX4H1+;$([N]([c])([C])[#6]);!$(N*=[#7,#8,#15,#16])]([#6])-[#6;X4][H]");
		queries.put("phosphoramide_pattern2", "[#6;A;X4;H1,H2,H3][#7;A;X3][P;X4]([#7;A;X3])([#7;A;X3])=[X1;O,S]");
		queries.put("n_demethoxymethylation", "[H][C;X4]([H])([H])[#8;X2]C([H])([H])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])][#6;a]");
		queries.put("methyl_carbon_adjacent_to_alkyne_group", "[H][#6;A;X4]([H])([H])[#6;A]#[#6;A]");
		queries.put("alicyclic_tertiary_amines_pattern1", "[H][#6;A;X4]!@-[#7;A;X3;+0!$([N]~[!#6])!$(N*=[#7,#8,#15,#16])](@-[#6;A])@-[#6;A]");
		queries.put("aryldimethylformamidine", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])-[#7;X2][#6;A;H1X3]=[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("carbon_adjacent_to_halogen_group", "[H][#6;X4]-[F,Cl,Br,I]");
		queries.put("phosphoramide_pattern1", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]P([#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])(=[O;X1])[#8]-[#1,CX4]");
		queries.put("ns_cleavage", "*-[#7;X3]!@-[#16;X2]-*");
		queries.put("_24_thiazolidinedione", "[H][#7;X3R1]-1-[#6;R1](=O)-[#6;R1]-[#16;R1;SX2,$([SX3]=[OX1]),$([SX3+]-[OX1-])]-[#6;R1]-1=O");
		queries.put("aryldimethylformamidine", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])-[#7;X2][#6;A;H1X3]=[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("arylhydroxylamine_pattern3", "[#8;H1X2]!@-[#7;H1X3]!@-[#6;R1]1=,:[#6][#7;R1]=,:[#6;X3R1][#7;R1]1");
		queries.put("hydroxylation_of_benzene_para_to_edg", "[H][#6;R1]=,:1[#6;R1](-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])])=,:[#6;R1][#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*])])=,:[#6;R1][#6;R1]=,:1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("arylhydroxylamine_pattern2", "[#8;H1X2]-[#7;H1X3]-[#6;X3R1]=,:1[#7;R1][#6]=,:[#6][#7;R1]=,:1");
		queries.put("arene_epoxidation_pattern2", "[H][c;X3R1]1-,:[c;X3R1](-[*,#1;!$([OX2H1])])[c;X3R1](-[*,#1])-,:[c;X3R1](-[*,#1])[c;R1](-[*,#1;!$([OX2H1])])-,:[c;X3R1]1-[F,Cl,Br,I,#1]");
		queries.put("hydroxylation_of_aliphatic_tertiary_penultimate_carbon", "[H][C;X4]([H])([H])[C;X4]([H])([#6])[#6]");
		queries.put("dicarboximide", "[#6,#1;A]-,:[#6;X3](=[O;R0])[#7;X3](!@-[#6;X4][H])[#6;X3]=[O;R0]");
		queries.put("m_hydroxylation_of_monosubstituted_benzene", "[H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[$(C([Br])([Br])[Br]),$(C([Cl])([Cl])[Cl]),$(C([F])([F])[F]),$(C([I])([I])[I]),$(C#[NX1]),$([NX1+]),$([#16;A](=O)(=O)[OX2H1,OX1-]),$([NX4H3+1]),$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])]),$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]),$([$([CX3H][#6]),$([CX3H2])]=[OX1]),$([#6][CX3](=[OX1])[#6]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[ClX1]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]),$([CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])])])=,:[#6;R1]([H])[#6;R1]=,:1[H]");
		queries.put("organophosphosulfur_compounds_pattern2", "[#6;a]-[#8;X2R0][P;X4]([#8])([#8])=[!#1!#6;X1;R0]");
		queries.put("hydroxylation_of_alicyclic_secondary_carbon_pattern5", "[H][#6;A]1([#1,!#6])[#6](-[*,#1])-,=[#6]-,=[#6]-,=[#6]-,=[#6]-,=[#6]1-[*,#1]");
		queries.put("propargylated_sec_and_tert_amines", "[H]C#C[C;R0]([H])([H])[#7;A;X3;+0!$([N]~[!#6])!$([N]*~[#7,#8,#15,#16])]([#6])[#1,#6]");
		queries.put("aryldimethylamine_pattern2", "[#6;A;H3X4][#7;A;X3]([#6;A;H3X4])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#7;R1]=,:[#6;R1][#6;R1]=,:1");
		queries.put("aryldimethylamine_pattern1", "[#6;A;H3X4][#7;A;X3]([#6;A;H3X4])[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6,#7;R1]1");
		queries.put("aliphatic_non_terminal_carbon_adjacent_to_aromatic_ring", "[H][#6;A;!$([CX4H3])](!@-[#6;a])(!@-[#6,#1;A;R0])!@-*");
		queries.put("thiocarbamate", "[#6]-[#16;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#1,#6])-[#1,#6]");
		queries.put("reduction_of_ketone_to_alcohol", "[#6]-[#6;X3](-[#6])=[O;X1]");
		queries.put("sterol_17_alpha_hydroxylation_pattern2", "[H][#6;A;X4]1([#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]2-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1][#6,#8,#7,#16]-[#6]-4=[#6]-[#6](=O)-[#6,#8,#7,#16]-[#6,#8,#7,#16]-[#6,#8,#7,#16]-4([#1,#6;A])~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-[#6,#8,#7,#16]12[#1,#6;A])[#6;R0]([#6;A;H3X4])=O");
		queries.put("sterol_17_alpha_hydroxylation_pattern1", "[H][#8]-[#6]-1-[#6,#8,#7,#16]-[#6,#8,#7,#16]-[#6,#8,#7,#16]-2([#1,#6;A])~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-[#6,#8,#7,#16]4([#1,#6;A])[#6,#8,#7,#16;A;R2](-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16][#6;A;H1X4]4[#6;R0]([#6;A;H3X4])=O)-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1][#6]=[#6]-2-[#6,#8,#7,#16]-1");
		queries.put("sterol_21_hydroxylation", "[#6;A;H3X4R0][#6;R0][#6;A;H1X4]1[#6,#8,#7,#16]-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16;A;R2]~3~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~4~[#6,#8,#7,#16]~3[#6;A;X4;H2,H1][#6,#8,#7,#16]-,=[#6,#8,#7,#16]1-,=2");
		queries.put("dithiocarbamate", "[#1,#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3](-[#1,#6])-[#1,#6]");
		queries.put("o_hydroxylation_of_monosubstituted_benzene", "[H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*]),F,Cl,Br,I])=,:[#6;R1]([H])[#6;R1]=,:1[H]");
		queries.put("aliphatic_hydroxylation_of_acyclic_12_disubstituted_alkene_pattern1", "[H][#6;A;X4]([#6;X3]=*)[#6;A;X3](-[H])=[#6;A;X3](-[H])[#6;X4]");
		queries.put("n_dealkylation_of_hydrazide", "[#1,#6][#7;X3](-[#6;X4][H])[#7;X3](-[#1,#6])[#6;X3]([#1,#6])=[O;X1]");
		queries.put("_4p_hydroxylation_of_oxazaphosphorine", "[#7;X3R0;NX3H1+0,NX3H2+0,$([NX3+0;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])])][P;X4]1(=[O;X1])[#7;NX3H1,$([NX3]-[#6])]-[#6;R1;CX4H2,$([CX4H1]-[#6])]-[#6;R1]-[#6;R1]-[#8]1");
		queries.put("benzene_ortho_to_edg", "[H][#6;R1]1=,:[#6;R1](-[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*])])[#6;R1]=,:[#6;R1][#6;R1](-[*,#1;!$([OX2H1])!$([OX2]-[CX4H3])])=,:[#6;R1]1-[*,#1;!$(C([Br])([Br])[Br])!$(C([Cl])([Cl])[Cl])!$(C([F])([F])[F])!$(C([I])([I])[I])!$(C#[NX1]),$([NX1+])!$([#16;A](=O)(=O)[OX2H1,OX1-])!$([NX4H3+1])!$([NX4H0+;$([N]([#6])([#6])[#6][#6]);!$([N]*~[#7,#8,#15,#16])])!$([$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8])]");
		queries.put("o_alkylated_group", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#8;X2R0]-[#6;$([CX4;!$(C([OX2])[O,S,#7,#15])]),c]");
		queries.put("methyl_carbon_adjacent_to_aromatic_ring", "[H][#6;A;X4]([H])([H])c:*");
		queries.put("s_adjacent_to_aliphatic_nitrogens", "[*,#1][#7;A]([*,#1])[#16;X2][#7;A]([*,#1])[*,#1]");
		queries.put("alkene_epoxidation_pattern1", "[#6][#6;A;X3R0]([#1,#6])=[#6;A;X3R0](/[#1,#6])[#1,#6]");
		queries.put("n_alkoxy_n_aryl_n_chloroacetamides", "[H][#6;X4]-[#8;X2][C;R0]([H])([H])[#7;X3R0](-[#6;R1]=,:1[#6;R1](-[#6;X4])=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1-[#6;X4])-[#6;R0](=O)[C;R0]([H])([H])Cl");
		queries.put("terminal_desaturation_pattern1", "[H][#6;X4R0](-[#6])[C;X4]([H])([H])[H]");
		queries.put("organophosphorothioate", "[#6][#8,#16;A;X2][P;X4]([#8])([#8])=[S;v2X1]");
		queries.put("secondary_heteroalicyclic_carbon_pattern1", "[H][#6;A;X4R]([H])(@-[#6;R])@-[#7,#16]");
		queries.put("guanidine_or_aminohydroazone", "[H][#7;A;X2]=[#6;X3]([#7;A;X3]([#1,#6,#7])[#1,#6,#7])[#7;A;X3]([#1,#6,#7])[#1,#6,#7]");
		queries.put("_5_aryl_14_benzodiazepines", "[!#1!#6]=,:[#6]1[#7;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6](=[#7;R1][C;R1]1([H])[#1,#6])!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1");
		queries.put("_4_alkyl_hydroxylation_of_oxoazaphosphorine", "[H][C;R1]1([H])[#6;R1]-[#6;R1]-[#8;R1][P;R1]([#7])(=[O;X1R0])[#7;X3R1]1");
		queries.put("s_oxidation_of_sulfoxide_to_sulfone", "[#6][S;X3]([#6])=[O;X1]");
		queries.put("o_deisopropylation", "[H][#6;A;X4]([#6;A;H3X4])([#6;A;H3X4])[#8;X2R0]-[#6;$([CX4;!$(C([OX2])[O,S,#7,#15])]),c]");
		queries.put("secondary_heteroalicyclic_carbon_pattern2", "[H][#6;A;X4R]([H])(@-[#6,#7,#8,#16])@-[#6;R]@-[#8,#7,#16]");
		queries.put("n_aryl_n_dimethylurea", "[#6;A;H3X4][#7;X3]([#6;A;H3X4])-[#6;X3](=[O;X1])[#7;A;H1X3][#6;R1]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1");
		queries.put("hydroxylation_of_aliphatic_secondary_penultimate_carbon", "[H][C;X4]([H])([H])C([H])([H])[#6]");
		queries.put("methylenedioxy_ring_opening", "[H][C;R1]1([H])[#8]-[#6;R2]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2-[#8]1");
		queries.put("n_dealkylation_of_ureas", "[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3;!$([#7][!#6])][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][H]");
		queries.put("p_substituted_anilide", "[H][*;A;CX4,OX2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#7]([H])-[#6](-[#1,#6])=O");
		queries.put("oxidative_deboronation", "[#6]!@-[#5;X3]([#8;A;H1X2])[#8;A;H1X2]");
		queries.put("alicyclic_tertiary_carbon", "[#6;R]@-[#6;A;H1X4R](@-[#6;A;R])@-[#6;A;R]");
		queries.put("alkyl_formamide", "[H][#7;X3](-[#6;X4])!@-[#6;X3]([H])=[O;X1]");
		queries.put("monosubstituted_benzene_p", "[$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([OX2H1,OX1-])=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1]1),$([#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1]([OX2H1,OX1-])[#6;R1]=,:[#6;R1]1),$([NX3H0+0,NX4H1+;$([N]([#6])([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H1+0,NX4H2+;$([N]([#6])[#6]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]),$([NX3H2+0,NX4H3+]c),$([OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]),$([#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]-[CX3;$([R0][#6]),$([H1R0])]=[OX1]),$([OX2][CX3;$([R0][#6]),$([H1R0])](=[OX1])[#6;!$(C=[O,N,S])]),$([CX4]),$([CX3]),$([#6;R1]1:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]:,-[#6;R1]:,=[#6;R1]1),$([#6;X3](=,:[#6;X3])-,:[#1,*]),$([#6;A;X3](=[#6;X3])-[#1,*]),F,Cl,Br,I]-[#6;R1]=,:1[#6;H1X3R1]=,:[#6;H1X3][#6;H1X3R1]=,:[#6;H1X3R1][#6;H1X3R1]=,:1");
		queries.put("thioamide", "[#1,$([#6;!$(C=[O,N,S])])][#7;X3]([#1,$([#6;!$(C=[O,N,S])])])[$([#6;X3][#6]),$([#6;X3;H1])]=[S;X1]");
		queries.put("organosulfurous_compound", "[#6]-[#8;X2][S;X3](=[O;X1])[#8;X2]-[#6]");
		queries.put("terminal_aliphatic_trihalide", "[H][C;X4](*)([F,Cl,Br,I])[F,Cl,Br,I]");
		queries.put("vinyl_ether", "[#6]-[#8;X2]-[#6;X3](-[#1,#6])=,:[#6;X3](-[#1,#6])-[#1,#6]");
		queries.put("primary_arylamine", "[H][NX3H2+0,NX4H3+;!$([N][#6]=[!#6;!#1]);!$([N]#[!#6,!#1])]c");
		queries.put("alkylsulfonamide", "[H][#6;A;X4;!$(C([OX2])[O,S,#7,#15])][#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4]([#1,#6])(=[O;X1])=[O;X1]");
		queries.put("carbon_alpha_to_conjugated_carbonyl", "[#8,#7,#6,#16,#1,#17,#9,#35,#53]-,:[#6;X3](=[O;X1])[#6]([H])=,:[#6]([*,#1])[*,#1]");
		queries.put("aliphatic_penultimate_carbon_adjacent_to_aromatic_carbon", "[H][C;X4]([H])([H])C([H])([#6;a])[*,#1]");
		queries.put("n_substituted_piperidine_ring", "[H]C1([H])[#6](-[*,#1])C([H])([H])C([H])([H])[#7;A;X3]([#1,c,$([CX4;A])])C1([H])[H]");
		queries.put("aromatic_carbon_para_to_halide_group", "[H]c:*:*:c-[F,Cl,Br,I]");
		queries.put("n_alkyl_diacylurea_pattern2", "[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3](-[#6;X4][H])[#6]([#1,#6])=O");
		queries.put("n_alkyl_diacylurea_pattern1", "[#1,#6][#6](=O)[#7;X3;!$([#7][!#6])][#6;X3](=[O;X1])[#7;X3;!$([#7][!#6])]-[#6;X4][H]");
		queries.put("halotrifluoroethane_pattern2", "[H][C;X4]([#17,#9])([F,Cl,Br,I])[C;X4](F)(F)F");
		queries.put("n_in_strained_ring_system", "[#6;R]@-[#7;X3R](@-[#6;R])@-[#6;R]");
		queries.put("o_methylaryl_pattern1", "[H]*[#6](=,:*[H])-[#8;X2R0]C([H])([H])[H]");
		queries.put("acyclic_secondary_amine", "[H][#6;A;X4][#7;A;H1X3R0;+0!$(N*=[#7,#8,#15,#16])][#6]");
		queries.put("carbon_gamma_to_conjugated_carbonyl", "[H][#6;X4]-[#6](-[*,#1])=[#6]([H])-[#6;X3]([#8,#7,#6,#16,#1,#17,#9,#35,#53;A])=[O;X1]");
		queries.put("bisallyl_group", "[H][C;X4](!@-[CX4,#1])([#6;A]([*,#1;#1,#6,Br,Cl,F,I])=[#6;A](/[*,#1;#1,#6,Br,Cl,F,I])[*,#1;#1,#6,Br,Cl,F,I])[#6;A]([*,#1;#1,#6,Br,Cl,F,I])=[#6;A](/[*,#1;#1,#6,Br,Cl,F,I])[*,#1;#1,#6,Br,Cl,F,I]");
		queries.put("n_chloroethyl_group", "[H]C([H])(Cl)C([H])([H])[#7;X3](-[#1,#6])-*");
		queries.put("tertiary_carboxamide_pattern2", "[#6;!$(C=[O,N,S])][#7;X3](!@-[#6;!$(C=[O,N,S])][H])[#6]([#1,#6])!@=[O;R0]");
		queries.put("aromatic_carbon_ortho_to_halide_group", "[H]c(:*):c-[F,Cl,Br,I]");
		queries.put("fused_benzene_ring_pattern2", "[#6]@[#6;R2]1=,:[#6;R2](@[#6,#8,#7,#16])[#6;R1]([H])=,:[#6;R1][#6;R1]=,:[#6;R1]1[H]");
		queries.put("terminal_methyl", "[H]C([H])([H])[#6;X4R0;!$([CX4H3])]");
		queries.put("_4_substitutted_phenol_pattern2", "[#1,CX4H3]-[#8]!@-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](!@-[!$([#6]-[#7]=[#7+]);!$([OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]);!$([#8;A;X2]S(=O)(=O)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)F);!$([#8;A;X2]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F);!$([#8]-[#16](=O)(=O)-[#6]1=,:[#6;R1]-,:[#6;R1]=,:[#6;R1](C([H])([H])[H])-,:[#6;R1]=,:[#6;R1]1);!$([#8]S(=O)(=O)C([H])([H])[H]),$([I,Br,Cl,F]);!$([#8X3+]([H])([H])),$([#8+]([H])[H]);!$([#8X2]-[#7X3+](=[OX1])[OX1-]);!$([#8X2]-[#15](=[O])(=[O])[OX2H1,OX1-]);!$([#8X2]-[#16](=[O])(=[O])[OX2H1,OX1-]);!$([#16X2]-[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]);!$([#16X3+]([H])[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]);!$([NX3+0,NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]);!$([#8;X2]-[#6](-[#1,#6])=[#8;X1]);!$([#8]-[#6]1=,:[#6;R1]-,:[#6;R1]=,:[#6;R1]-,:[#6;R1]=,:[#6;R1]1);!$([#8;X1-]);!$([#1])])=,:[#6;R1][#6;R1]=,:1");

		
		
		/**
		 * MACCS 166 Keys from CDK 1.5.13 (maccs.txt) / RDKit (MACCSkeys.py)
		 */
		// 		
//		queries.put("maccs_1","?");
//		queries.put("maccs_2","[#104]");
//		queries.put("maccs_3","[#32,#33,#34,#50,#51,#52,#82,#83,#84]");
//		queries.put("maccs_4","[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]");
//		queries.put("maccs_5","[Sc,Ti,Y,Zr,Hf]");
//		queries.put("maccs_6","[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]");
//		queries.put("maccs_7","[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]");
		queries.put("maccs_8","[!#6;!#1]1~*~*~*~1");
//		queries.put("maccs_9","[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]");
//		queries.put("maccs_10","[Be,Mg,Ca,Sr,Ba,Ra]");
		queries.put("maccs_11","*1~*~*~*~1");
//		queries.put("maccs_12","[Cu,Zn,Ag,Cd,Au,Hg]");
		queries.put("maccs_13","[#8]~[#7](~[#6])~[#6]");
		queries.put("maccs_14","[#16]-[#16]");
		queries.put("maccs_15","[#8]~[#6](~[#8])~[#8]");
		queries.put("maccs_16","[!#6;!#1]1~*~*~1");
		queries.put("maccs_17","[#6]#[#6]");
		queries.put("maccs_18","[#5,#13,#31,#49,#81]");
		queries.put("maccs_19","*1~*~*~*~*~*~*~1");
//		queries.put("maccs_20","[#14]");
		queries.put("maccs_21","[#6]=[#6](~[!#6;!#1])~[!#6;!#1]");
		queries.put("maccs_22","*1~*~*~1");
		queries.put("maccs_23","[#7]~[#6](~[#8])~[#8]");
		queries.put("maccs_24","[#7]-[#8]");
		queries.put("maccs_25","[#7]~[#6](~[#7])~[#7]");
		queries.put("maccs_26","[#6]=;@[#6](@*)@*");
//		queries.put("maccs_27","[I]");
		queries.put("maccs_28","[!#6;!#1]~[CH2]~[!#6;!#1]");
//		queries.put("maccs_29","[#15]");
		queries.put("maccs_30","[#6]~[!#6;!#1](~[#6])(~[#6])~*");
		queries.put("maccs_31","[!#6;!#1]~[F,Cl,Br,I]");
		queries.put("maccs_32","[#6]~[#16]~[#7]");
		queries.put("maccs_33","[#7]~[#16]");
		queries.put("maccs_34","[CH2]=*");
//		queries.put("maccs_35","[Li,Na,K,Rb,Cs,Fr]");
//		queries.put("maccs_36","[#16R]");
		queries.put("maccs_37","[#7]~[#6](~[#8])~[#7]");
		queries.put("maccs_38","[#7]~[#6](~[#6])~[#7]");
		queries.put("maccs_39","[#8]~[#16](~[#8])~[#8]");
		queries.put("maccs_40","[#16]-[#8]");
		queries.put("maccs_41","[#6]#[#7]");
//		queries.put("maccs_42","F");
		queries.put("maccs_43","[!#6;!#1;!H0]~*~[!#6;!#1;!H0]");
		queries.put("maccs_44","[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]");
		queries.put("maccs_45","[#6]=[#6]~[#7]");
//		queries.put("maccs_46","Br");
		queries.put("maccs_47","[#16]~*~[#7]");
		queries.put("maccs_48","[#8]~[!#6;!#1](~[#8])(~[#8])");
//		queries.put("maccs_49","[!+0]");
		queries.put("maccs_50","[#6]=[#6](~[#6])~[#6]");
		queries.put("maccs_51","[#6]~[#16]~[#8]");
		queries.put("maccs_52","[#7]~[#7]");
		queries.put("maccs_53","[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]");
		queries.put("maccs_54","[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]");
		queries.put("maccs_55","[#8]~[#16]~[#8]");
		queries.put("maccs_56","[#8]~[#7](~[#8])~[#6]");
//		queries.put("maccs_57","[#8R]");
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
//		queries.put("maccs_84","[NH2]");
		queries.put("maccs_85","[#6]~[#7](~[#6])~[#6]");
		queries.put("maccs_86","[C;H2,H3][!#6;!#1][C;H2,H3]");
		queries.put("maccs_87","[F,Cl,Br,I]!@*@*");
//		queries.put("maccs_88","[#16]");
		queries.put("maccs_89","[#8]~*~*~*~[#8]");
//		queries.put("maccs_90","[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]");
//		queries.put("maccs_91","[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_92","[#8]~[#6](~[#7])~[#6]");
		queries.put("maccs_93","[!#6;!#1]~[CH3]");
		queries.put("maccs_94","[!#6;!#1]~[#7]");
		queries.put("maccs_95","[#7]~*~*~[#8]");
		queries.put("maccs_96","*1~*~*~*~*~1");
		queries.put("maccs_97","[#7]~*~*~*~[#8]");
		queries.put("maccs_98","[!#6;!#1]1~*~*~*~*~*~1");
		queries.put("maccs_99","[#6]=[#6]");
		queries.put("maccs_100","*~[CH2]~[#7]");
//		queries.put("maccs_101","[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]");
		queries.put("maccs_102","[!#6;!#1]~[#8]");
//		queries.put("maccs_103","Cl");
		queries.put("maccs_104","[!#6;!#1;!H0]~*~[CH2]~*");
//		queries.put("maccs_105","*@*(@*)@*");
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
//		queries.put("maccs_128","[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]");
//		queries.put("maccs_129","[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]");
		queries.put("maccs_130","[!#6;!#1]~[!#6;!#1]");
//		queries.put("maccs_131","[!#6;!#1;!H0]");
		queries.put("maccs_132","[#8]~*~[CH2]~*");
		queries.put("maccs_133","*@*!@[#7]");
//		queries.put("maccs_134","[F,Cl,Br,I]");
		queries.put("maccs_135","[#7]!:*:*");
		queries.put("maccs_136","[#8]=*");
//		queries.put("maccs_137","[!C;!c;R]");
		queries.put("maccs_138","[!#6;!#1]~[CH2]~*");
//		queries.put("maccs_139","[O;!H0]");
//		queries.put("maccs_140","[#8]");
//		queries.put("maccs_141","[CH3]");
//		queries.put("maccs_142","[#7]");
		queries.put("maccs_143","*@*!@[#8]");
		queries.put("maccs_144","*!:*:*!:*");
		queries.put("maccs_145","*1~*~*~*~*~*~1");
//		queries.put("maccs_146","[#8]");
//		queries.put("maccs_147","[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]");
		queries.put("maccs_148","*~[!#6;!#1](~*)~*");
//		queries.put("maccs_149","[C;H3,H4]");
		queries.put("maccs_150","*!@*@*!@*");
//		queries.put("maccs_151","[#7;!H0]");
		queries.put("maccs_152","[#8]~[#6](~[#6])~[#6]");
		queries.put("maccs_153","[!#6;!#1]~[CH2]~*");
		queries.put("maccs_154","[#6]=[#8]");
		queries.put("maccs_155","*!@[CH2]!@*");
		queries.put("maccs_156","[#7]~*(~*)~*");
		queries.put("maccs_157","[#6]-[#8]");
		queries.put("maccs_158","[#6]-[#7]");
//		queries.put("maccs_159","[#8]");
//		queries.put("maccs_160","[C;H3,H4]");
//		queries.put("maccs_161","[#7]");
//		queries.put("maccs_162","a");
		queries.put("maccs_163","*1~*~*~*~*~*~1");
//		queries.put("maccs_164","[#8]");
//		queries.put("maccs_165","[R]");
//		queries.put("maccs_166","?");	
		
		
		
		/**
		 * MACCS 322 - Yannick's SMARTS implementation of MACCS 322 keys, taken from MayaChemTools
		 * http://www.mayachemtools.org/docs/scripts/html/MACCSKeysFingerprints.html
		 */
//		queries.put("maccs_322_001", "*~*~*~*");
//		queries.put("maccs_322_002", "[!#1!#6]");	
	//		queries.put("maccs_322_003", "");	
//		queries.put("maccs_322_004", "*~*~*~*~*");	
//		queries.put("maccs_322_005", "[$(*~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6])]");
//		queries.put("maccs_322_006", "[$(*~[!#1!#6]~[!#1!#6]~[!#1!#6]),$([!#1!#6]~*~[!#1!#6]~[!#1!#6])]");	
//		queries.put("maccs_322_007", "[!#1!#6][H]");	
		queries.put("maccs_322_008", "[$(*~*~[#6]([H])[H]),$(*~C(~*)([H])[H])]");		
		queries.put("maccs_322_009", "*~C([H])([H])[H]");
	//		queries.put("maccs_322_010", "[F,Cl,Br,I]");	//halogen (already here)
//		queries.put("maccs_322_011", "[$(*-*(-*)-*)]");	 // [$(),$(*-*(-*)-*)]
	//		queries.put("maccs_322_012", "");	
//		queries.put("maccs_322_013", "[$(*@-*(@-*)@-*)]"); // [$(),$(*@-*(@-*)@-*)]
		queries.put("maccs_322_014", "*@-*!@-*@-*");	
		queries.put("maccs_322_015", "*!:*:*!:*");	
		queries.put("maccs_322_016", "*!@-*!@-*");
		queries.put("maccs_322_017", "*!@-*@-*!@-*");
		queries.put("maccs_322_018", "*:*!:*:*");	
	//		queries.put("maccs_322_019", "");	// Same as organoheterocycle
		queries.put("maccs_322_020", "[$(*-*(-*)(-*)(-*)-*),$(*@-*(@-*)(@-*)(@-*)@-*),$([!#1!#6!#7!#8!#16!#9!#17!#35!#53])]");	
//		queries.put("maccs_322_021", "[!+0]");
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

		
		/**
		 * PubChem keys 264 to 881 (The index xxx in the pubchem_xxx is based on a start from 0)
		 * ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
		 * The first 263 keys cannot be mapped to an atom (e.g. 261	>= 4 aromatic rings)
		 */
		
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

		return queries;
	}
	

	
	
	
	public LinkedHashMap<String,List<List<Integer>>> generateCustomBoMFingeprint(IAtomContainer atc, 
			LinkedHashMap<String, String> patterns, IChemObjectBuilder builder) throws CDKException{
		LinkedHashMap<String,List<List<Integer>>> atomLists = new LinkedHashMap<String,List<List<Integer>>>();
		
		for (Map.Entry<String, String> item : patterns.entrySet()) {
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);
			if(smartsPattern.matches(atc)){
//				List<List<Integer>> l = new ArrayList<List<Integer>>();
//				for(List<Integer> i : smartsPattern.getMatchingAtoms()){
//					l.add(i);
//				}
				
				atomLists.put(item.getKey(), smartsPattern.getMatchingAtoms());
			}
			else{
				atomLists.put(item.getKey(), null);
			}
		}
		

		return atomLists;
	}
	
	public static void main(String[] args) throws Exception{
		IChemObjectBuilder 	builder 	= SilentChemObjectBuilder.getInstance();
		SmilesGenerator 	smiGen 		= new SmilesGenerator().isomeric();
		SmilesParser		smiParser	= new SmilesParser(builder);
		
		ChemSearcher csearcher = new ChemSearcher();
		
		
		LinkedHashMap<String, String> patterns = getCustomBoMFingerprintPatterns();
		
		IAtomContainer atc =  smiParser.parseSmiles("CCC(C)CCCCC(=O)NC(CCNCS([O-])(=O)=O)C(=O)NC(C(C)O)C(=O)NC(CCNCS([O-])(=O)=O)C(=O)NC1CCNC(=O)C(NC(=O)C(CCNCS([O-])(=O)=O)NC(=O)C(CCNCS([O-])(=O)=O)NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)C(CCNCS([O-])(=O)=O)NC1=O)C(C)O");
		LinkedHashMap<String,List<List<Integer>>> res = csearcher.generateCustomBoMFingeprint(atc, patterns, builder);
		
		for(String s : res.keySet()){
			System.out.println(s + "\t" + res.get(s));
		}
		
		System.out.println(patterns.size());
	}
	
}
