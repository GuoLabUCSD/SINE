#! /bin/bash

#Get required arguments
usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script will calculate a PHBR score for each relevant event and output a list of potential neoantigens.

Arguments:
-i      Path to a Directory of isoforms generated from main_ROI.sh
-a      Path to a directory of allele affinity scores generated from xls_parser.sh
-m      Path to a directory containing patient files with mapping ID number to netmhcpan generated from get_netmhcpan.sh
-o      Location to store phbr score results and plots
-p	Pipeline Directory
EOF
  exit
}

while getopts "i:a:m:o:p:h" flag ; do

        case "${flag}" in
		i) isoforms_dir=${OPTARG};;
                a) affinity_scores=${OPTARG};;
                m) mapping_directory=${OPTARG};;
                o) output_directory=${OPTARG};;
                p) pipeline_directory=${OPTARG};;
                h) usage;;
        esac
done

if [ -z $affinity_scores ] || [ -z $mapping_directory ] || [ -z $output_directory ] || [ -z $pipeline_directory ] || [ -z $isoforms_dir ]; then
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

python $pipeline_directory/supplemental/calc_ROI_phbr.py \
--isoforms_dir $isoforms_dir \
--affinity_scores $affinity_scores \
--mapping_directory $mapping_directory \
--output_directory $output_directory
