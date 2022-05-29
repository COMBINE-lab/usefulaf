#!/bin/bash

display_usage() { 
  echo "This script returns the gene id to gene name mapping as a TSV file using gffread"
  echo -e "\nUsage: $0 [options]" 
  echo -e "\toptions:" 
  echo -e "\t -g, --gtf REQUIRED path to a GTF file"
  echo -e "\t -o, --output REQUIRED path to the output TSV file (will be created if it doesn't exist)"
  echo -e "\t --gffread, path to the gffread binary"
  echo -e "\t -h, --help display this help message"
} 

v3_plist="https://umd.box.com/shared/static/eo0qlkfqf2v24ws6dfnxty6gqk1otf2h"
v2_plist="https://umd.box.com/shared/static/jbs2wszgbj7k4ic2hass9ts6nhqkwq1p"

while [[ "$#" -gt 0 ]]; do 
  case $1 in 
    -g|--gtf) gtf="$2"; shift;;
    -o|--output) output="$2"; shift ;;
    --gffread) gffread="$2"; shift ;;
    -h|--help) help=1; shift;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac 
  shift 
done

if [[ -n "$help" ]]; then
  display_usage
  exit 0
fi

if [[ -z "$output" || -z "$gtf" ]]; then
  display_usage
  exit 1
fi

if [[ -z "$gffread" ]]; then 
gffread="gffread"
fi

# make temp dir
WORK_DIR=`mktemp -d`

# check if tmp dir was created
if [[ ! "$WORK_DIR" || ! -d "$WORK_DIR" ]]; then
    echo "Could not create temp dir"
    exit 1
fi

# deletes the temp directory
function cleanup {      
    rm -rf "$WORK_DIR"
    # echo "  - Deleted temp working directory $WORK_DIR"
}

# generate gff from gtf
gff_cmd="$gffread ${gtf} -o $WORK_DIR/genes.gff"
eval $gff_cmd

# make the file
grep "gene_name" $WORK_DIR/genes.gff | cut -f9 | cut -d';' -f2,3 | sed 's/=/ /g' | sed 's/;/ /g' | cut -d' ' -f2,4 | sort | uniq > ${output}

trap cleanup EXIT
