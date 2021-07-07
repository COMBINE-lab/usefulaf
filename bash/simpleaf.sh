#!/bin/bash

display_usage() { 
  echo "This script runs alevin-fry to quantify a single-cell RNA-seq experiment"
  echo -e "\nUsage: $0 [options]" 
  echo -e "\toptions:" 
  echo -e "\t -1, --r1 REQUIRED comma separated list of left reads"
  echo -e "\t -2, --r2 REQUIRED comma separated list of right reads"
  echo -e "\t -i, --index REQUIRED path to a (sparse or dense) salmon splici index"
  echo -e "\t -o, --output REQUIRED path to output directory (will be created if it doesn't exist)"
  echo -e "\t -f, --fmode REQUIRED permit list filter mode, one of {knee, k, unfilt, u}"
  echo -e "\t -c, --chem REQUIRED chemistry of experiment, one of {v2, v3}"
  echo -e "\t -r, --res REQUIRED resolution strategy for alevin-fry, one of {cr-like, cr-like-em}"
  echo -e "\t -m, --t2g REQUIRED three-column txp-to-gene file to pass to alevin-fry quant command"
  echo -e "\t -t, --threads OPTIONAL number of threads to use when running [default: min(16, num cores)]"
  echo -e "\t -h, --help display this help message"
} 


salmon="${SALMON_BIN:-salmon}"
fry="${FRY_BIN:-alevin-fry}"

time="/usr/bin/time -v -o"
threads=16
np=`nproc`

if [ "$threads" -gt "$np" ]; then
  threads="$np"
fi

while [[ "$#" -gt 0 ]]; do 
  case $1 in 
    -1|--r1) read1="$2"; shift ;;
    -2|--r2) read2="$2"; shift ;;
    -i|--index) index="$2"; shift ;;
    -t|--threads) threads="$2"; shift ;;
    -o|--output) output="$2"; shift ;;
    -f|--fmode) fmode="$2"; shift ;;
    -r|--res) res="$2"; shift ;;
    -c|--chem) chem="$2"; shift ;;
    -m|--t2g) t2g="$2"; shift ;;
    -h|--help) help=1; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac 
  shift 
done

set -o errexit -o pipefail

if [[ -n "$help" ]]; then
  display_usage
  exit 0
fi


if [[ -z "$read1" || -z "$read2" || -z "$index" || -z "$threads" || -z "$output" || -z "$fmode" || -z "$res" || -z "$chem" || -z "$t2g" ]]; then
  display_usage
  exit 1
fi

## Check that the fitler mode is one of knee,k,unfilt,u
if [[ "$fmode" == "k" || "$fmode" == "knee" || "$fmode" == "unfilt" || "$fmode" == "u" ]]; then
  echo "filter mode is : $fmode"
else
  echo "filter mode must be one of {knee, k, unfilt, u}"
  exit 1
fi

## Check that the chemistry is either v2 or v3
if [[ "$chem" == "v2" || "$chem" == "v3" ]]; then
  echo "chemistry is : 10x chromimum $chem"
else
  echo "chemistry mode must be one of {v2, v3}"
  exit 1
fi

## If the chemistry is v2, set the chemflag and the fitler flags
if [[ $chem == "v2" ]]; then
  if [[ "$fmode" == "unfilt" || "$fmode" == "u" ]]; then
    v2file=$output/plist/10x_v2_permit.txt
    if [ ! -f "$v2file" ]; then
          echo "10x v2 permit list does not exist, downloading now"
          bash get_10x_permit_lists.sh -o $output/plist -l v2 
    fi
    permitmode="-u $output/plist/10x_v2_permit.txt"
  else
    permitmode="-k"
  fi
  chemflag="--chromium"
fi

## If the chemistry is v3, set the chemflag and the fitler flags
if [[ $chem == "v3" ]]; then
  if [[ "$fmode" == "unfilt" || "$fmode" == "u" ]]; then
    v3file=$output/plist/10x_v3_permit.txt
    if [ ! -f "$v3file" ]; then
          echo "10x v3 permit list does not exist, downloading now"
          bash get_10x_permit_lists.sh -o $output/plist -l v3 
    fi
    permitmode="-u $output/plist/10x_v3_permit.txt"
  else
    permitmode="-k"
  fi
  chemflag="--chromiumV3"
fi

mkdir -p $output/logdir
logdir="$output/logdir"

## turn comma separated list into space separated list
read1=`echo $read1 | tr ',' ' '`
read2=`echo $read2 | tr ',' ' '`

## map
cmd="$time $logdir/map.time $salmon alevin -l ISR -i $index -1 $read1 -2 $read2 -p $threads $chemflag -o $output/alevin_map --sketch"
echo "mapping:"
echo "command: $cmd"
echo "============="
eval $cmd

### generate permit list
cmd="$time $logdir/gpl.time $fry generate-permit-list $permitmode -d fw -i $output/alevin_map -o $output/gpl/ |& stdbuf -oL tr '\r' '\n' > $logdir/gpl.log"
echo "gpl:"
echo "command: "$cmd
echo "============="
eval $cmd

### collate
cmd="$time $logdir/collate.time $fry collate -i $output/gpl/ -r $output/alevin_map -t $threads |& stdbuf -oL tr '\r' '\n' > $logdir/collate.log"
echo "collate:"
echo "command: "$cmd
echo "============="
eval $cmd

### quant
cmd="$time $logdir/quant.time $fry quant -r $res --use-mtx -m $t2g -i $output/gpl/ -o $output/quant -t $threads |& stdbuf -oL tr '\r' '\n' > $logdir/quant.log"
echo "quant:"
echo "command: "$cmd
echo "============="
eval $cmd

echo "Finished! Quantification results are at $output/quant"
