#!/bin/bash

v3_plist="https://umd.box.com/shared/static/eo0qlkfqf2v24ws6dfnxty6gqk1otf2h"
v2_plist="https://umd.box.com/shared/static/jbs2wszgbj7k4ic2hass9ts6nhqkwq1p"

while [[ "$#" -gt 0 ]]; do 
  case $1 in 
    -o|--output) output="$2"; shift ;;
    -l|--list) list="$2"; shift;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac 
  shift 
done

if [ "$list" == "v2" ]; then
  mkdir -p $output
  wget -v -O $output/10x_v2_permit.txt -L $v2_plist
fi

if [ "$list" == "v3" ]; then
  mkdir -p $output
  wget -v -O $output/10x_v3_permit.txt -L $v3_plist
fi
