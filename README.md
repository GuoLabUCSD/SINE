# SINE:
SINE (Splice Isoform Neoantigen Evaluator) aims to identify potential neoantigens produced as a result of aberrant alternative splicing in tumor tissue and calculate their PHBR binding affinity score through the use of Trinity, netMHCpan, user-provided junctions/regions of interest, and sample BAM files. For a full tutorial and list of arguments to be used with SINE, please see the documentation here: https://github.com/GuoLabUCSD/SINE/tree/main/doc


## Pre-Requisites:
 * BAM and corresponding BAI Index files for samples of interest
 * STAR SJ.out.tab files for samples of interest


## Minimum Software Requirements:
 * bedtools
	* Available at: https://github.com/arq5x/bedtools2
	* Note: bedtools should be available in your PATH variable
 * NetMHCpan
	* Available at: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
 * Python 3 (Tested on Version 3.10.13) with the following modules installed (Recommended to install through pip):
	* pandas (Tested on Version 1.5.3)
	* pysam
   	* argparse
	* biopython
	* numpy
	* matplotlib
	* seaborn
 * SAMtools
 	* Available at: https://www.htslib.org/
 * Trinity 
	* Available at: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#running-trinity-using-singularity


## Installation:
 1. git clone https://github.com/GuoLabUCSD/SINE
 2. From inside the Neoantigen_Pipeline directory users may need to run the below for each function's bash file
	* chmod u+x {file ending in .sh}


  ## Acknowledgements and Citations
-xls_parse.py contains modified code from Rachel Marty's NetMHCpan wrapper PyPresent (https://github.com/Rachelmarty20/pypresent)

	Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PMID: 21572440; PMCID: PMC3571712.

	Marty R, Kaabinejadian S, Rossell D, Slifker MJ, van de Haar J, Engin HB, de Prisco N, Ideker T, Hildebrand WH, Font-Burgada J, Carter H. MHC-I Genotype Restricts the Oncogenic Mutational Landscape. Cell. 2017 Nov 30;171(6):1272-1283.e15. doi: 10.1016/j.cell.2017.09.050. Epub 2017 Oct 26. PMID: 29107334; PMCID: PMC5711564.

	Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res. 2020 Jul 2;48(W1):W449-W454. doi: 10.1093/nar/gkaa379. PMID: 32406916; PMCID: PMC7319546.

