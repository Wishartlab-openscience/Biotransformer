# Important Notice
***************************************************************************************************

BioTransformer's environmental microbial degradation module uses data from the EAWAG's Biodegradation and Biocatalysis Database, which is licensed by [EnviPath](https://envipath.com/license/) under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-nc-sa/4.0/). Therefore users of the environmental microbial degradation module must cite the EnviPath paper below (See Cite section).
To use the environmental microbial module for commercial purposes, users must request an appropriate commercial license from [EnviPath](https://envipath.org/).



# BioTransformer's README

***************************************************************************************************

This is version 3.0.0 of BioTransformer. BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota, as well as the soil/aquatic microbiota. BioTransformer also assists scientists in metabolite identification, based on the metabolism prediction. Please make sure to download the folders database and supportfiles, and save them in the same folder as the .jar file.

BioTransformer is offered to the public as a freely accessible software package under the GNU LGPL, version 3 license. Users are free to copy and redistribute the material in any medium or format. Moreover, they could modify, and build upon the material under the condition that they must give appropriate credit, provide links to the license, and indicate if changes were made. Furthermore, the above copyright notice and this permission notice must be included. Use and re-distribution of the these resources, in whole or in part, for commercial purposes requires explicit permission of the authors. We ask that all users of the BioTransformer software tool, the BioTransformer web server, or BioTransformerDB to cite the BioTransformer reference in any resulting publications, and to acknowledge the authors.


### Cite:

1. Djoumbou Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R, and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small Molecule Metabolism Prediction and Metabolite Identification; Journal of Cheminformatics; (2019) 11:2; DOI: 10.1186/s13321‐018‐0324‐5.


### Contact

To report issues, provide feedback, or ask questions, please send ane-mail the following address: djoumbou@ualberta.ca

***************************************************************************************************

Usage:
java -jar biotransformer-3.0.0.jar -k <BioTransformer Task> -b <BioTransformer Option> [-ismi, -imol, or -isdf <Input file>] 
       [-imol, or -isdf <Output>] [-s <Number of steps>] [-a <Annotate>] [-m <Masses>] [-r <Masses>] [-f <Formulas>] [-q <Sequence>]
       [-cm <Cyp450 Mode>]

This is the version 3.0.0 of BioTransformer. BioTransformer is a software
tool that predicts small molecule metabolism in mammals, their gut
microbiota, as well as the soil/aquatic microbiota. BioTransformer also assists scientists in metabolite identification, based on the metabolism prediction.

 -a,--annotate                       Search PuChem for each product, and
                                     annotate with CID and synonyms

 -b,--btType <BioTransformer Option>   The type of description: Type of
                                     biotransformer - EC-based  (ecbased),
                                     CYP450 (cyp450), Phase II (phaseII),
                                     Human gut microbial (hgut),
                                     human super transformer** (superbio,
                                     or allHuman), Environmental microbial (envimicro)*.
                                     
-cm, cyp450mode						CYP450 prediction Mode here: 
									1) CypReact + BioTransformer rules; 2) CyProduct only;
									3) Combined: CypReact + BioTransformer rules + CyProducts.
									Default mode is 1.  
									
-h,--help                           Prints the usage.

-ismi, --ismiles < Input >             Read the input as a single molecule in 
                                     SMILES format

-imol, --molinput <Input file>       Read the input as a single molecule in
                                     the specified MOL format

-isdf, --sdfinput <Input file>       Read the input from the specified SDF 
                                     format 
                                                                         
-k, --task <BioTransformer Task>     The task to be permed: pred for
                                     prediction, or cid for compound
                                     identification

-m,--mList                           Given the starting compound(s), Find
                                     all their metabolites with the
                                     specified masses, and show a
                                     metabolism pathway. A
                                     whitespace-separated list of masses
                                     is expected.
                                     
-ocsv, --csvoutput <Output file>     Save the results into the specified CSV 
                                     file 

-osdf, --sdfoutput <Output file>     Save the results into the specified SDF 
                                     file 

-f,--formulas <formulas>             Semicolon-separated list of formulas
                                     of compounds to identify                                     
-s,--nsteps <Number of steps>       The number of steps for the
                                     prediction. This option can be set by
                                     the user for the EC-based, CYP450,
                                     Phase II, and Environmental microbial
                                     biotransformers. The default value is 1.

-q, --bsequence <Sequence>			 Ordered Sequence of biotransformation steps.
									  Semi-colon separated pairs of biotransformer
									  types and corresponding number of steps
									  to be simulated. 

-t,--mTolerance                     Mass tolerance for metabolite
                                     identification (default is 0.01).  
                            
                         
                                               
(* ) While the 'superbio' option runs a set number of transformation steps in a
pre-defined order (e.g. deconjugation first, then Oxidation/reduction,
etc.), the 'allHuman' option predicts all possible metabolites from any
applicable reaction(Oxidation, reduction, (de-)conjugation) at each step.


(** ) For the environmental microbial biodegradation, all reactions (aerobic and anaerobic) are reported, and not only the aerobic biotransformations (as per default in the EAWAG BBD/PPS system).


**************
Examples:
**************

1) To predict the biotransformation of a molecule from an SDF input using the human super transformer (option superbio) and annotate the metabolites with names and database IDs (from PubChem), run

      java -jar biotransformer-3.0.0.jar -k pred -b superbio -isdf -i #{input file name} -osdf #{output file} -a.

2) To predict the 2-step biotransformation of Thymol (a monoterpene) using the human super transformer (option allHuman) using the SMILES input, cypMode = 3 (combined mode), and saving to a CSV file, run

      java -jar biotransformer-3.0.0.jar  -k pred -b allHuman -ismi "CC(C)C1=CC=C(C)C=C1O" -ocsv #{replace with output file name} -s 2 -cm 3

Currently, the outputfile is SDF per default.

3) Identify all human metabolites (max depth = 2) of Epicatechin ("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1") with masses 292.0946 Da and 304.0946 Da, with a mass tolerance of 0.01 Da. Provide an annotation (Common name, synonyms, and PubChem CID), when available.

      java -jar biotransformer-3.0.0.jar  -k cid -b allHuman -ismi "O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1" -osdf #{replace with output file name} -s 2 -m "292.0946;304.0946" -t 0.01 -a
    
    - DO NOT forget the quotes around the SMILES string or the list of masses.

4) Simulate an order sequence of metabolism of Atrazine ("CCNC1=NC(=NC(=N1)Cl)NC(C)C"), starting with two steps of Cyp450 oxidation, followed by one step of conjugation.

	java -jar biotransformer-3.0.0.jar -ismi "CCNC1=NC(=NC(=N1)Cl)NC(C)C" -osdf ~/atrazine-sequence.sdf -k pred -q "cyp450:2; phaseII:1"
To report issues, provide feedback, or ask questions, please send an
e-mail the following address: djoumbou@ualberta.ca


**************************************************
Compiling and packaging BioTransformer with Maven:
**************************************************

To compile BioTransformer and package it into an executable .jar file, you need to have maven (v >3.6) installed. Then run:

``` 
1. cd biotransformer
2. ./install-via-maven.sh
```

