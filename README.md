# Neoantigen_Pipeline:
This tool aims to identify potential neoantigens produced as a result of aberrant alternative splicing in tumor tissue and calculate their PHBR binding affinity score through the use of Trinity, netMHCpan, user-provided junctions of interest, and sample BAM files.


## Pre-Requisites:
 * BAM and corresponding BAI Index files for samples of interest
 * STAR SJ.out.tab files for samples of interest


## Minimum Software Requirements (Recommended to install via pip):
 * bedtools
	* Available at: https://github.com/arq5x/bedtools2
	* Note: bedtools should be available in your PATH variable
 * NetMHCpan
	* Available at: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
 * Python 3 (Tested on Version 3.6.8) with the following modules installed:
	* pandas (Tested on Version 1.1.5)
	* pysam
	* biopython
	* numpy
	* matplotlib
	* seaborn
 * SAMtools
 * Singularity
 * Trinity 
	* Available at: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#running-trinity-using-singularity


## Installation:
 1. git clone https://github.com/GuoLabUCSD/Neoantigen_Pipeline
 2. From inside the Neoantigen_Pipeline directory users may need to run the below for each function's bash file
	* chmod u+x {file ending in .sh}


## Example Workflow:
 ## 1. Get wild-type junctions
    
	./get_wt.sh -i junctions_list.txt -b BAM_File_Directory/ -n my_experiment -o wt_output_directory -p ../Neoantigen_Pipeline

### Arguments:

		-i      
			New line separated list of junctions of interest to evaluate
		-b      
			Path to a Directory containing BAM files for the samples of interest and STAR's SJ.tab.out files
		-n      
			Prefix to what the output files should be named
		-o      
			Output Directory to store files
		-p
			Path to the Neoantigen_Pipeline Directory

### Output:
		{Sample_Name}_junctions.tsv:
			New-line separated list of junctions of interest and their corresponding wild-type junction per individual sample.
	
		{Filename_Prefix}_junctionsOfInterest_plusWT.tsv:
			Table of junctions of interest and their corresponding wild-type junction for all samples with relevant data extracted from STAR's SJ.out.tab files.

 ## 2. Extract junction reads, perform Trinity De novo assembly, and translate all reading frames of the present isoforms. 
 -It is recommended to run this in parallel on a cluster for all samples in the experiment before proceeding to step 3.
    
	./main.sh -b Sample_Name.bam -j Sample_Name_junctions.tsv -g Mus.musculus.GRCm39.107.gtf -f Mus_musculus.GRCm39.cds.all.fa -t trinityrnaseq_latest.sif -o isoforms_output_directory -p ../Neoantigen_Pipeline

### Arguments:
		
		-b
			Path to Tumor RNA BAM file aligned to the genome for a single sample
		-j
			Path to Junctions of interest file
		-g
			Path to GTF annotation file
		-f
			Path to CDS fasta file
		-t
			Path to Trinity Singularity Image
		-o
			Path to ouput directory where results should be stored
		-p
			Path to the directory where this pipeline is

### Output:
		..{isoforms_output_directory}/{Sample_Name/} Directory containing:
			{Sample_Name}.identify_junctions.sh: 
				Bash file with the Trinity execution
			intermediate_results: 
				{chromosome}:{start-end} directory with Trinity output
			neopeptides.tsv:
				Tab separated table containing the assembled isoform, junction position (python-indexed), peptides (kmer length = 15), and the transcript ID

## 3. Prep and Run NetMHCPan to get peptide binding to MHC class I. 
-It is recommended to run this in parallel on a cluster for all samples in the experiment before proceeding to step 4.

-See ./examples to see how the allele table should be formatted.
   
	./get_netmhcpan.sh -i isoforms_output_directory -a mhc_allele_table.txt -f output_directory_fastas -m output_directory_mapping_files -o output_directory_netmhcpan_results -s Sample_Name -p ../Neoantigen_Pipeline

### Arguments:
		
		-i
			Path to a directory of samples containing isoform results from step 2. Ie) Output Directory from step 2 is the input directory here
		-f
			Path to an output directory to store fasta files of isoforms. This directory should have subdirectories labelled "mhc_i" and "mhc_ii"
		-m
			Path to an output directory to store mapping files.
		-o
			Path to an output directory to store NetMHCpan results
		-a
			Path to a tab-separated text file of MHC alleles to check peptide binding against for the sample
		-s
			Name of one of the sample directories in the isoform directory input (-i input) to run through NetMHCpan
		-p
			Path to the directory where this pipeline is

### Output:
		../{output_directory_fastas}/mhc_i/{Sample_Name}.fasta
			Fasta file containing the netmhcpan id # for isoform and the isoform sequence
		../{output_directory_fastas}/mhc_i/{Sample_Name}.peptides
			Fasta file containing the netmhcpan id # of the isoform and all possible k-mers of length 8-11 spanning the junction.
		../{output_directory_fastas}/mhc_ii/{Sample_Name}.fasta
			Fasta file containing the netmhcpan id # for isoform and the isoform sequence
		../{output_directory_fastas}/mhc_ii/{Sample_Name}.peptides
			Fasta file containing the netmhcpan id # of the isoform and all possible k-mers of length 15 spanning the junction.
		../{output_directory_mapping_files}/{Sample_Name}.tsv
			Tab separated file containing the same results from step 2, but with an additional NetMHCpan id # to match peptides to their corresponding isoform and junction event after running NetMHCpan.
		../{output_directory_netmhcpan_results}/{Sample_Name}.xlsoutput
			Tab separated output table from NetMHCPan containing binding affinity scores for every possible peptide at every provided junction for the sample
			
## 4. Parse the NetMHCpan .xls results files
-See ./examples to see how the allele table should be formatted.
   
	./xls_parser.sh -i output_directory_netmhcpan_results -f output_directory_fastas/mhc_i -o output_directory_parsed_affinity -a mhc_allele_table -p ../Neoantigen_Pipeline

### Arguments:

		-i
			Path to a directory of NetMHCpan .xlsoutput results files
		-f
			Path to a directory of Fasta files containing the isoforms. Ie) The mhc_i subdirectory in The fasta output directory from step 3 is the input directory here.
		-o
			Path to a directory to store newly parsed affinity files
		-a
			Path to a tab-separated text file of MHC alleles to check peptide binding against for the sample
		-p
			Path to the directory where this pipeline is

### Output:
		{Sample_Name}.output:
			Tab separated output table containing the affinity scores for each allele for the highest binding peptide per isoform

## 5. Calculate the PHBR score for each peptide, output results table, and provide a plot of ratios of event expression
   
	./get_phbr.sh -j Filename_Prefix_junctionsOfInterest_plusWT.tsv -a output_directory_parsed_affinity -m output_directory_mapping_files -o output_directory_final_results -p ../Neoantigen_Pipeline

### Arguments:
		
		-j
			Path to the full junctions table generated from Step 1.
		-a
			Path to a directory of allele affinity scores. The output directory from step 4 is the input directory here
		-m
			Path to a directory containing sample files with mapping IDs
		-o
			Path to a directory to store the file results tables and plots
		-p
			Path to the directory where this pipeline is

### Output:
		fin_prioritization.high_express_junctions.tsv:
			Tab separated results table of PHBR scores and corresponding peptide results for the highest binders
		
		fin_sorted_prioritization.high_express_junctions.tsv:
			Tab separated results table of PHBR scores, corresponding peptide results, and the percent of samples containing the potential neoantigen/peptide. Sorted numerically by junction of interest.
		
		fin_lhbr_filtered_prioritization.high_express_junctions.tsv:
			Tab separated results table filtered where the PHBR of the alternative junction of interest is less than the PHBR of the wild-type junction (Lower PHBR = Higher Binding Affinity)
	
		fin_exp_and_phbr_filtered_prioritization.high_express_junctions.tsv
			Tab separated results table of filtered both by PHBR and where the expression level of the junction of interest is greater than the expression level of the wild-type junction
		
		ase_vs_other.png:
			Boxplot comparing overall PHBR scores between the alternatively spliced junctions of interest and the wild-type junctions

  ## Acknowledgements and Citations
-xls_parse.py contains modified code from Rachel Marty's NetMHCpan wrapper PyPresent (https://github.com/Rachelmarty20/pypresent)

	Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PMID: 21572440; PMCID: PMC3571712.

	Marty R, Kaabinejadian S, Rossell D, Slifker MJ, van de Haar J, Engin HB, de Prisco N, Ideker T, Hildebrand WH, Font-Burgada J, Carter H. MHC-I Genotype Restricts the Oncogenic Mutational Landscape. Cell. 2017 Nov 30;171(6):1272-1283.e15. doi: 10.1016/j.cell.2017.09.050. Epub 2017 Oct 26. PMID: 29107334; PMCID: PMC5711564.

	Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res. 2020 Jul 2;48(W1):W449-W454. doi: 10.1093/nar/gkaa379. PMID: 32406916; PMCID: PMC7319546.

