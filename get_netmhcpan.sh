#! /bin/bash

#Get required arguments
usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script will prep and run netmhcpan for one sample in the isoforms directory provided by main.sh.

Arguments:

-i      Path to a Directory of isoforms generated from Main.sh
-f      Location to store output fasta file
-m      Location to store mapping files
-o      Output Directory to store files
-a      Filepath to a tab-separated text file of MHC alleles to check peptide binding against for each sample
-s      Patient name to run netmhcpan on
-p	Pipeline Directory
EOF
  exit
}

while getopts "i:f:m:o:a:s:p:h" flag ; do

        case "${flag}" in
                i) isoforms_dir=${OPTARG};;
                f) output_directory_fastas=${OPTARG};;
                m) output_directory_mapping=${OPTARG};;
                o) netmhcpan_output_directory=${OPTARG};;
                a) allele_table=${OPTARG};;
                s) patient_sample=${OPTARG};;
		p) pipeline_directory=${OPTARG};;
                h) usage;;
        esac
done

if [ -z $isoforms_dir ] || [ -z $output_directory_fastas ] || [ -z $output_directory_mapping ] || [ -z $netmhcpan_output_directory ] || [ -z $allele_table ] || [ -z $patient_sample ] || [ -z $pipeline_directory ]; then
        echo "Missing required arguments"
        usage
        exit 1
fi

# make mhc_i and mhc_ii directories if necessary
if [ ! -d $output_directory_fastas/mhc_i ]; then
        mkdir -p $output_directory_fastas/mhc_i
fi

# make output directory if necessary
if [ ! -d $output_directory_fastas/mhc_ii ]; then
        mkdir -p $output_directory_fastas/mhc_ii
fi


#Run Python Script
echo "Running"


python $pipeline_directory/supplemental/prepare_netmhcpan.py \
--isoforms_dir $isoforms_dir \
--output_directory_fastas $output_directory_fastas \
--output_directory_mapping $output_directory_mapping \
--netmhcpan_output_directory $netmhcpan_output_directory \
--allele_table $allele_table \
--patient_sample $patient_sample \
--pipeline_directory $pipeline_directory

#Run Bash Script

bash $pipeline_directory/supplemental/temp_files/$patient_sample\_netmhcpan.sh

#Remove Bash Script

rm $pipeline_directory/supplemental/temp_files/$patient_sample\_netmhcpan.sh

