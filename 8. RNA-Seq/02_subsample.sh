#!/bin/bash

input_dir='/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/1_trimmed/subsample_mappingtest'


R1=($input_dir/*.pair1.truncated)
R2=($input_dir/*.pair2.truncated)
for ((i=0;i<=${#R1[@]};i++))
do
  seqkit sample -p 0.1 -s 11 -o "${R1[i]}.subsample" ${R1[i]}
  seqkit sample -p 0.1 -s 11 -o "${R2[i]}.subsample" ${R2[i]}
done
