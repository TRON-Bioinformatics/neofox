#!/bin/bash

# this script calls input tool and queues annotation of identidied mps for each patient separately

# $1 = directory of icam output folders
# $2 = output folder
# $3 = allele file
# $4 = cohort name

for dir in $1/Pt*;
do
  if [ -f $dir/scratch/*_mut_set.txt.transcript.squish.somatic.freq ]
  then
    w="$dir"
    # extract patient id from path
    num="$(echo $w | tr -cd "/" | wc -c)"
    #echo $num
    sum=$(($num+1))
    echo $dir
    #echo $sum
    pat="$(cut -d'/' -f$sum <<<$dir)"
    echo $pat
    file=$dir/scratch/$pat
    file2="$file"_mut_set.txt.transcript.squish.somatic.freq
    outfile="$2"/"$4".input."$pat".txt
    output="$2"/"$pat".log
    error="$2"/"$pat".err
    echo $file2
    echo $outfile
    echo $output
    echo $error
    module load anaconda/2/2018
    module load R/3.6.0
    module load bioinf/netMHCpan/4.0
    sbatch -p Compute --job-name INPuT_"$4"_$pat -c 1 --time=30-00:00:00 --mem 2G --output=$output --error=$error --wrap="python /projects/SUMMIT/WP1.2/input/development/predict_all_epitopes.py -i $file2 -a $3 > $outfile"
  fi;
done
