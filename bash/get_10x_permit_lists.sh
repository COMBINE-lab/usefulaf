#!/bin/bash

display_usage() { 
  echo "This script downloads a 10x chromium v2 or v3 permit list to the specified output directory"
  echo -e "\nUsage: $0 [options]" 
  echo -e "\toptions:" 
  echo -e "\t -o, --output REQUIRED path to output directory (will be created if it doesn't exist)"
  echo -e "\t -l, --list REQUIRED permit list to download, one of {v2, v3}"
  echo -e "\t -h, --help display this help message"
} 

v3_plist="https://umd.box.com/shared/static/eo0qlkfqf2v24ws6dfnxty6gqk1otf2h"
v2_plist="https://umd.box.com/shared/static/jbs2wszgbj7k4ic2hass9ts6nhqkwq1p"

while [[ "$#" -gt 0 ]]; do 
  case $1 in 
    -o|--output) output="$2"; shift ;;
    -l|--list) list="$2"; shift;;
    -h|--help) help=1; shift;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac 
  shift 
done

if [[ -n "$help" ]]; then
  display_usage
  exit 0
fi

if [[ -z "$output" || -z "$list" ]]; then
  display_usage
  exit 1
fi

if [ "$list" == "v2" ]; then
  mkdir -p $output
  wget -v -O $output/10x_v2_permit.txt -L $v2_plist
fi

if [ "$list" == "v3" ]; then
  mkdir -p $output
  wget -v -O $output/10x_v3_permit.txt -L $v3_plist
fi

echo -e "\n\noutput written to ${output}/10x_${list}_permit.txt\n\n"
