#!/usr/bin/bash

file=orthogroups_Tfas_Tlei_Acom.counts.no_TEs.fromN0.txt
# All orthogroups
aco_all=`awk -F'\t' '{sum+=$2;}END{print sum;}' $file`
tfas_all=`awk -F'\t' '{sum+=$3;}END{print sum;}' $file`
tlei_all=`awk -F'\t' '{sum+=$4;}END{print sum;}' $file`
# Singlecopy
aco_one=`awk -F'\t' '$2 == 1 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_one=`awk -F'\t' '$3 == 1 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_one=`awk -F'\t' '$4 == 1 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
aco_til_one=`awk -F'\t' '$3 == 1 && $4 == 1 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_til_one=`awk -F'\t' '$3 == 1 && $4 == 1 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_til_one=`awk -F'\t' '$3 == 1 && $4 == 1 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
bro_one=`awk -F'\t' '$2 == 1 && $3 == 1 && $4 == 1 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
# Multi copy
aco_m=`awk -F'\t' '$2 > 1 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_m=`awk -F'\t' '$3 > 1 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_m=`awk -F'\t' '$4 > 1 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
aco_mT=`awk -F'\t' '$3 > $4 && $4 > 0 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_mT=`awk -F'\t' '$3 > $4 && $4 > 0 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_mT=`awk -F'\t' '$3 > $4 && $4 > 0 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
aco_mL=`awk -F'\t' '$3 < $4 && $3 > 0 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_mL=`awk -F'\t' '$3 < $4 && $3 > 0 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_mL=`awk -F'\t' '$3 < $4 && $3 > 0 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
aco_mE=`awk -F'\t' '$3 == $4 && $3 > 1 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_mE=`awk -F'\t' '$3 == $4 && $3 > 1 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_mE=`awk -F'\t' '$3 == $4 && $3 > 1 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`
# Unique Genes
aco_u=`awk -F'\t' '$2 != 0 && $3 == 0 && $4 == 0 {print $0}' $file | awk -F'\t' '{sum+=$2;}END{print sum;}'`
tfas_u=`awk -F'\t' '$2 == 0 && $3 != 0 && $4 == 0 {print $0}' $file | awk -F'\t' '{sum+=$3;}END{print sum;}'`
tlei_u=`awk -F'\t' '$2 == 0 && $3 == 0 && $4 != 0 {print $0}' $file | awk -F'\t' '{sum+=$4;}END{print sum;}'`

# print
echo '	Acom	Tfas	Tlei'
echo "total_orthogroups	${aco_all}	$tfas_all	$tlei_all"
echo "single_copy_species	$aco_one	$tfas_one	$tlei_one"
echo "single_copy_bro	$bro_one"
echo "single_copy_til	$aco_til_one	$tfas_til_one	$tlei_til_one"
echo "multi_copy	$aco_m	$tfas_m	$tlei_m"
echo "multi_copy_T	$aco_mT	$tfas_mT	$tlei_mT"
echo "multi_copy_L	$aco_mL	$tfas_mL	$tlei_mL"
echo "multi_copy_E	$aco_mE	$tfas_mE	$tlei_mE"
echo "uniq	$aco_u	$tfas_u	$tlei_u"
