#!/bin/bash
bam_file=${1}
list_bam_files=${2}
te_gff=${3}
output_prefix=${4}    
thread=${5}
overlap=${6}
expectation=${7}
max_iterations_outer=${8}
max_iterations_inner=${9}
abe_step=${10}
SCRIPT_DIR=${11}
PARENT_DIR=${12}
reserve_tmp_dir=${13}

#
# echo "Executing: batch_mode.sh"
current_path=$(pwd) # && echo ${current_path}

# print the parameters
# echo "Parameters: "
# echo ${bam_file} ${list_bam_files} ${te_gff} ${output_prefix} ${thread} ${overlap} ${expectation} ${max_iterations_outer} ${max_iterations_inner} ${abe_step} ${SCRIPT_DIR} ${PARENT_DIR} ${reserve_tmp_dir} | awk '{for(i=1;i<=NF;i++) printf "\t%s\n", $i}'

# invoking base model
# echo -e ${list_bam_files} 
echo -e ${list_bam_files} | xargs -I {} -P ${thread} ${PARENT_DIR}/scripts/base_mode.sh {} "BatchMode" ${te_gff} ${output_prefix} ${thread} ${overlap} ${expectation} ${max_iterations_outer} ${max_iterations_inner} ${abe_step} ${SCRIPT_DIR} ${PARENT_DIR} ${reserve_tmp_dir}

# union RowName
ls ${output_prefix}*_LATTE_TE_Expression | head -n 1 | xargs -i awk -v FS="\t" 'NR>1 {print $1}' {} > ${output_prefix}_LATTE_ReadCount
ls ${output_prefix}*_LATTE_TE_Expression | head -n 1 | xargs -i awk -v FS="\t" 'NR>1 {print $1}' {} > ${output_prefix}_LATTE_TPM
[[ -f ${output_prefix}_LATTE_ColName ]] && rm ${output_prefix}_LATTE_ColName

# build Matrix
echo -e ${list_bam_files} | while read AbsPath_bam_file   
do
    # 
    sample_name=$(basename $AbsPath_bam_file)
    echo $sample_name >> ${output_prefix}_LATTE_ColName
    # 纵向合并
    awk -v FS="\t" 'NR>1 {print $2}' ${output_prefix}${sample_name}_LATTE_TE_Expression > ${output_prefix}tmpLATTE_ReadCount
    paste -d "\t" ${output_prefix}_LATTE_ReadCount ${output_prefix}tmpLATTE_ReadCount > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_ReadCount
    awk -v FS="\t" 'NR>1 {print $3}' ${output_prefix}${sample_name}_LATTE_TE_Expression > ${output_prefix}tmpLATTE_TPM
    paste -d "\t" ${output_prefix}_LATTE_TPM ${output_prefix}tmpLATTE_TPM > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_TPM
    # rm
    rm ${output_prefix}${sample_name}_LATTE_TE_Expression
done

# 删除中间文件
rm ${output_prefix}tmpLATTE_ReadCount ${output_prefix}tmpLATTE_TPM

# process the output file
awk '{printf "%s\t", $0}' ${output_prefix}_LATTE_ColName | awk '{gsub(/\t$/, "\n", $0); printf "\t%s", $0}'  > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_ColName
awk -v OFS="\t" '{for(i=1;i<=NF;i++){ if(i==1){printf "`%s`\t", $i}else{printf "%s%s", $i, (i == NF ? "\n" : OFS) }} }' ${output_prefix}_LATTE_TPM > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_TPM
awk -v OFS="\t" '{for(i=1;i<=NF;i++){ if(i==1){printf "`%s`\t", $i}else{printf "%s%s", $i, (i == NF ? "\n" : OFS) }} }' ${output_prefix}_LATTE_ReadCount > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_ReadCount

# format for eQTL mapping
cat ${output_prefix}_LATTE_ColName ${output_prefix}_LATTE_TPM > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_TPM
cat ${output_prefix}_LATTE_ColName ${output_prefix}_LATTE_ReadCount > ${output_prefix}tmp_LATTE && mv ${output_prefix}tmp_LATTE ${output_prefix}_LATTE_ReadCount

