#! /bin/bash

#Get required arguments
usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script will run a modified version of pypresent and parse the affinity results from netmhcpan.

Arguments:

-i      Path to a file of junctions
-f      Directory where peptide fasta files are stored per sample
-o      Output Directory to store parsed affinity results
-p	Pipeline Directory
EOF
  exit
}

while getopts "i:f:o:p:h" flag ; do

        case "${flag}" in
                i) netmhcpan_results_directory=${OPTARG};;
                f) peptide_directory=${OPTARG};;
                o) output_directory=${OPTARG};;
		p) pipeline_directory=${OPTARG};;
                h) usage;;
        esac
done

if [ -z $netmhcpan_results_directory ] || [ -z $peptide_directory ] || [ -z $output_directory ] || [ -z $pipeline_directory ]; then
        echo "Missing required arguments"
        usage
        exit 1
fi

if [ ! -d $output_directory ]; then
        echo "The specified output directory does not exist"
        exit 1
fi

#Run Python Script
echo "Running"

python $pipeline_directory/supplemental/xls_parse.py \
--netmhcpan_results_directory $netmhcpan_results_directory \
--peptide_directory $peptide_directory \
--output_directory $output_directory
