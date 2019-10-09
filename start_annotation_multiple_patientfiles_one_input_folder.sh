#!/bin/bash

# this script calls input tool and queues annotation of identidied mps for each patient separately

# $1 = directory with icam input files
# $2 = output folder
# $3 = allele file
# $4 = cohort name

for icamfile in $1/*.txt;
do
  # extract patient id from path
  echo $icamfile
  patfil=${icamfile##*/}
  #patfil="$(cut -d'/' -f12 <<<$icamfile)"
  echo $patfil
  pat="$(cut -d'.' -f1 <<<$patfil)"
  echo $pat
    #file=$dir/scratch/$pat
  file2=$icamfile
  outfile="$2"/"$4".input."$pat".txt
  output="$2"/"$pat".log
  error="$2"/"$pat".err
  echo $file2
  echo $outfile
  echo $output
  echo $error
  module load anaconda/2/2018
  module load bioinf/netMHCpan/4.0
  module load R/3.6.0
  sbatch -p Compute --job-name INPuT_"$4"_$pat -c 1 --time=30-00:00:00 --mem 2G --output=$output --error=$error --wrap="python /projects/SUMMIT/WP1.2/input/development/predict_all_epitopes.py -i $file2 -a $3 > $outfile"
done
