#! /bin/bash

# to suppress perl warnings
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

usage() {
  cat <<EOF 
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script runs the default pipeline for identifying alt splicing neopeptides from RNA 
to creating a custom database for searching in mass-spec proteomics data.

- Minimum required input: tumor RNA bam, ROI file (see wiki for format)


Available options:

-b    RNA bam aligned to the genome [REQUIRED]
-j    region of interest file [REQUIRED]
-g    path to gtf file [REQUIRED]
-f    path to cds fasta file [REQUIRED]
-t    path to trinity singularity image [REQUIRED]
-o    output directory [REQUIRED]
-p    pipeline directory [REQUIRED]
-s    sample_name [REQUIRED]
-l    integer specifying the minimum contig length for Trinity Fasta Output [OPTIONAL]
-r    integer specifying the factor with which to duplicate the number of reads by [OPTIONAL]
-d    regions of interest file also contains direction info [OPTIONAL]
-m    data is from a mouse genome [OPTIONAL]
EOF
  exit
}
contig_length=100
while getopts "b:j:g:f:t:o:p:s:l::r::dmh" flag ; do

	case "${flag}" in
		b) tumor_rna_bam=${OPTARG};;
		j) insertion_file=${OPTARG};;
		g) gtf_file=${OPTARG};;
		f) cds_fasta_file=${OPTARG};;
		t) trinity_image=${OPTARG};;
		p) pipeline_dir=${OPTARG};;
		o) output_dir=${OPTARG};;
		s) sample_name=${OPTARG};;
		l) contig_length=${OPTARG};;
		r) read_boost=${OPTARG};;
		m) mouse='This_is_a_Mouse_Genome';;
		d) direct='Data_contains_direction_info';;
		h) usage;;
	esac
done

# check required args are given
if [ -z $tumor_rna_bam ] || [ -z $insertion_file ] || [ -z $output_dir ] || [ -z $pipeline_dir ] || [ -z $gtf_file ] || [ -z $trinity_image ] || [ -z $cds_fasta_file ] || [ -z $sample_name ]; then
	echo "the default pipeline requires a tumor RNAseq bam pair, a ROI file, a designated output directory, a gtf file, the CDS fasta file, the trinity image, and the pipeline directory. Exiting..."
	usage
	exit 1
fi

# check input files exist
if [ ! -f $tumor_rna_bam ] || [ ! -f $insertion_file ]; then
	echo "Either $tumor_rna_bam or $insertion_file does not exist. Exiting...
	"
	exit 1
fi

# make output directory if necessary
if [ ! -d $output_dir ]; then
	mkdir -p $output_dir
fi


echo "tumor_rna_bam: $tumor_rna_bam"
echo "ROI_file: $insertion_file"
echo "output_dir: $output_dir
";



echo "Running default pipeline
"

temp=$(basename $tumor_rna_bam)
rna_genome_basename=${temp%.*}

echo "
####################################
##### 0. Checking file formats #####
####################################
"
# make sure ROI file has required columns 
if [ -z $direct ]; then
	expected_num=1
else
	expected_num=3
fi

obs_num1=$(grep -Eo 'ROI' $insertion_file | uniq | wc -l)
obs_num2=$(grep -Eo 'direction' $insertion_file | uniq | wc -l)
obs_num3=$(grep -Eo 'junction' $insertion_file | uniq | wc -l)
obs_num=$(($obs_num1 + $obs_num2 + $obs_num3))

if [ $expected_num == $obs_num ]; then
	echo "	- Insertions file appears to have all expected columns, continuing..."
else
	echo "	- Insertions file ($insertion_file) does not have all columns. See GitHub Wiki for more info. Exiting...
	"
	exit 1
fi

# warn if RNA bam seems to be aligned to transcriptoms
if [[ "$sample_name" == *"Aligned.toTranscriptome"* ]]; then
	echo "	- WARNING - BAM SEEMS TO BE ALIGNED TO TRANSCRIPTOME... ATTEMPTING TO PROCEED..."
fi	

# make sure RNA bam has index
if [ ! -f $tumor_rna_bam.bai ]; then
	echo "	- $tumor_rna_bam.bai does not exist, please create and retry. Exiting...
	"
	exit 1
fi

echo "
############################################
##### 1. Create script to run pipeline #####
############################################
"

# make intermediate results dir
if [ -z $read_boost ]; then
	python $pipeline_dir/alt_splice_neoantigens/scripts/prepare_trinity_ROI_scripts.py --bam_filepath $tumor_rna_bam --insertion_file $insertion_file --output_script_path $output_dir/$rna_genome_basename.identify_ROI.sh --work_dir $output_dir/intermediate_results --trinity_sif $trinity_image --contig_length $contig_length --pipeline_dir $pipeline_dir/alt_splice_neoantigens

else
	python $pipeline_dir/alt_splice_neoantigens/scripts/prepare_trinity_ROI_scripts.py --bam_filepath $tumor_rna_bam --insertion_file $insertion_file --output_script_path $output_dir/$rna_genome_basename.identify_ROI.sh --work_dir $output_dir/intermediate_results --trinity_sif $trinity_image --contig_length $contig_length --pipeline_dir $pipeline_dir/alt_splice_neoantigens --read_boost $read_boost
fi

echo "
##################################
##### 2. Quietly run Trinity #####
##################################
"

bash $output_dir/$rna_genome_basename.identify_ROI.sh

echo "
###################################
##### 3. Generate neopeptides #####
###################################
"
if [ -z $mouse ]; then
	if [ -z $direct ]; then
		python $pipeline_dir/alt_splice_neoantigens/scripts/human_generate_ROI_neopeptides.py -i $output_dir/intermediate_results/ -G $gtf_file -F $cds_fasta_file -s $sample_name -b $tumor_rna_bam
	else
		python $pipeline_dir/alt_splice_neoantigens/scripts/human_generate_ROI_neopeptides.py -i $output_dir/intermediate_results/ -G $gtf_file -F $cds_fasta_file -s $sample_name -b $tumor_rna_bam -j $insertion_file
	fi

else
	if [ -z $direct ]; then
		python $pipeline_dir/alt_splice_neoantigens/scripts/mouse_generate_ROI_neopeptides.py -i $output_dir/intermediate_results/ -G $gtf_file -F $cds_fasta_file -s $sample_name -b $tumor_rna_bam
	else
		python $pipeline_dir/alt_splice_neoantigens/scripts/mouse_generate_ROI_neopeptides.py -i $output_dir/intermediate_results/ -G $gtf_file -F $cds_fasta_file -s $sample_name -b $tumor_rna_bam -j $insertion_file
	fi

fi
