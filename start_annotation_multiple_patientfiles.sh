#!/bin/bash

# this script calls input tool and queues annotation of identidied mps for each patient separately

# $1 = directory of icam output folders
# $2 = output folder
# $3 = allele file

for dir in $1/Pt*;
do
  if [ -f $dir/scratch/*_mut_set.txt.transcript.squish.somatic.freq ]
  then
    pat="$(cut -d'/' -f10 <<<$dir)"
    file=$dir/scratch/$pat
    file2="$file"_mut_set.txt.transcript.squish.somatic.freq
    outfile="$2"/riaz.input."$pat".txt
    output="$2"/"$pat".log
    error="$2"/"$pat".err
    echo $file2
    echo $outfile
    echo $output
    echo $error
    sbatch -p CentOS --job-name INPuT_riaz_$pat -c 1 --time=30-00:00:00 --mem 2G --output=$output --error=$error --wrap="python /projects/SUMMIT/WP1.2/input/productive/addannot1.0/predict_all_epitopes.py -i $file2 -a $3 > $outfile"
  fi;
done
