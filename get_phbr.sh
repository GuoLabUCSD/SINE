#! /bin/bash

#Get required arguments
usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <args>

This script will calculate a PHBR score for each relevant junction event, compare the junctions of interest to the wild type, and provide relevant plots.

Arguments:

-j      Path to the full junctions table generated from get_wt.sh
-a      Path to a directory of allele affinity scores generated from xls_parser.sh
-m      Path to a directory containing patient files with mapping ID number to netmhcpan generated from get_netmhcpan.sh
-o      Location to store phbr score results and plots
-p	Pipeline Directory
EOF
  exit
}

while getopts "j:a:m:o:p:h" flag ; do

        case "${flag}" in
                j) junctions_table=${OPTARG};;
                a) affinity_scores=${OPTARG};;
                m) mapping_directory=${OPTARG};;
                o) output_directory=${OPTARG};;
                p) pipeline_directory=${OPTARG};;
                h) usage;;
        esac
done

if [ -z $junctions_table ] || [ -z $affinity_scores ] || [ -z $mapping_directory ] || [ -z $output_directory ] || [ -z $pipeline_directory ]; then
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

python $pipeline_directory/supplemental/calc_phbr.py \
--junctions_table $junctions_table \
--affinity_scores $affinity_scores \
--mapping_directory $mapping_directory \
--output_directory $output_directory
