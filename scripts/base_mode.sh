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

if [[ ${2} == "BatchMode" ]]; then
    file_name=$(basename ${bam_file})
    BacthMode_output_prefix="${output_prefix}${file_name}"
    output_prefix=${BacthMode_output_prefix}
fi

# print the parameters
echo "[being processed] ${bam_file}"
# echo "Parameters: "
# echo ${bam_file} ${list_bam_files} ${te_gff} ${output_prefix} ${thread} ${overlap} ${expectation} ${max_iterations_outer} ${max_iterations_inner} ${abe_step} ${SCRIPT_DIR} ${PARENT_DIR} ${reserve_tmp_dir}

# making temporary directory
tmp_dir=$(date +%s%N | sha256sum | awk '{print $1}')
mkdir ${tmp_dir}

# change directory
cd ${tmp_dir}

# preparing data
samtools view -e '[AS] && [NH]==1' ${bam_file} -b -o single.bam
samtools view -e '[AS] && [NH]!=1' ${bam_file} -b -o multi.bam
## bamtobed
bedtools bamtobed -split -i single.bam > single.bed
bedtools bamtobed -split -i multi.bam > multi.bed
## intersect
bedtools intersect -wa -wb -f ${overlap} -a single.bed -b ${te_gff} | awk '!seen[$0]++' > single_te.bed
bedtools intersect -wa -wb -f ${overlap} -a multi.bed -b ${te_gff} | awk '!seen[$0]++' > multi_te.bed

# outer loop cycle
out_cycle_time=1
EndingCycleCount=$(awk '{print $4}' multi_te.bed | sort -u | wc -l | awk '{print $0 / 10}' | cut -f 1 -d ".")
while [ $out_cycle_time -lt $max_iterations_outer ]; 
do
    # echo "outside cycle time " ${out_cycle_time}
    # echo "current expectation " ${expectation}
    [ ${out_cycle_time} -eq 1 ] && LastSettleReadCount=0
    out_cycle_time=$(expr ${out_cycle_time} + 1)
    
    # ByTESubfamily
    [ -e "ByTESubfamily" ] || mkdir ByTESubfamily && cd ByTESubfamily # && echo "ByTESubfamily"
    ## pre
    awk -v OFS="|" ' !seen[$0]++ {print$16 "\t" $1,$2,$3 "\t" $4}' ../multi_te.bed > MultiTEReads
    ## executing
    ${PARENT_DIR}/scripts/ByTESubfamily.sh ../single_te.bed MultiTEReads ${expectation} ${max_iterations_inner}
    ## connecting
    ${PARENT_DIR}/scripts/insideConnection.sh $out_cycle_time "ByTESubfamily" ../single_te.bed ../multi_te.bed ../RecordSettledMultiReads
    cd ../

    # ByTEConsensus_SubFamily
    [ -e "ByTEConsensus_SubFamily" ] || mkdir ByTEConsensus_SubFamily && cd ByTEConsensus_SubFamily # && echo "ByTEConsensus_SubFamily"
    ## reformatting
    awk -v OFS="|" ' !seen[$0]++ {print$7,$10,$11,$16,$17,$18 "\t" $1,$2,$3 "\t" $4}' ../multi_te.bed > MultiTEReads
    awk -v OFS="|" '{print$7,$10,$11,$16,$17,$18 "\t" $1,$2,$3 "\t" $4}' ../single_te.bed | awk -v OFS="\t" '{gsub(/\|/,"\t",$2); print $2,$1,$3}'  > SingleTEReads.bed
    awk -v OFS="\t" '{gsub(/\|/,"\t",$2); print $2,$1,$3}' MultiTEReads | sort -k5,5 > MultiTEReads.bed
    ## executing
    ${PARENT_DIR}/scripts/ByTEConsensus_SubFamily.sh ./SingleTEReads.bed ./MultiTEReads.bed ${expectation} ${max_iterations_inner}
    ## connecting
    ${PARENT_DIR}/scripts/insideConnection.sh $out_cycle_time "ByTEConsensus" ../single_te.bed ../multi_te.bed ../RecordSettledMultiReads
    cd ..

    # ByTEAnnotation
    # echo "ByTEAnnotation"
    [ -e "ByTEAnnotation" ] || mkdir ByTEAnnotation && cd ByTEAnnotation
    ## reformatting
    awk -v OFS="|" ' !seen[$0]++ {print$7,$10,$11,$16,$17,$18 "\t" $1,$2,$3 "\t" $4}' ../multi_te.bed > MultiTEReads
    ## executing
    ${PARENT_DIR}/scripts/ByTEAnnotation.sh ../single_te.bed ./MultiTEReads ${expectation} ${max_iterations_inner}
    ## connecting
    ${PARENT_DIR}/scripts/insideConnection.sh $out_cycle_time "ByTEAnnotation" ../single_te.bed ../multi_te.bed ../RecordSettledMultiReads
    cd ../

    # Ending the useless cycle
    TotalSettleReadCount=$(wc -l RecordSettledMultiReads | awk '{print $1}')
    CycleSettleReadCount=$(expr ${TotalSettleReadCount} - ${LastSettleReadCount})
    # echo "CurrentSettledReads" ${CycleSettleReadCount}
    if [[ ${CycleSettleReadCount} -lt ${EndingCycleCount} ]]; then
        break
    fi
    mid_value=$(expr ${CycleSettleReadCount} + ${LastSettleReadCount})
    LastSettleReadCount=${mid_value}
done

# ByJointProbability
[ -e "ByJointProbability" ] || mkdir ByJointProbability && cd ByJointProbability # && echo "ByJointProbability"
${PARENT_DIR}/scripts/ByJointProbability.sh ../multi_te.bed
## connecting
${PARENT_DIR}/scripts/insideConnection.sh $out_cycle_time "ByJointProbability" ../single_te.bed ../multi_te.bed ../RecordSettledMultiReads
cd ../

# AbnormalERVDetection.sh
if [ ${abe_step} == "True" ];then
        cp ${PARENT_DIR}/data/ERVFamily ERVFamily
        ${PARENT_DIR}/scripts/AbnormalERVDetection.sh ${te_gff} ./single_te.bed ${PARENT_DIR}
        [ -s "abnormal_ERVs" ] && cat abnormal_ERVs > ${output_prefix}_AbnormalERVs
fi

# Output
awk -v OFS="|" '{ gsub(/"Motif:|"/,"",$16); print$7,$10,$11,$16,$17,$18,$6 "\t" $1,$2,$3 "\t" $4}' single_te.bed > SingleTEReads
#
awk -v OFS="\t" '{ gsub(/"Motif:|"/,"",$16); te_length[$7"|"$10"|"$11"|"$16"|"$17"|"$18"|"$6]=$11-$10 ; sum[$7"|"$10"|"$11"|"$16"|"$17"|"$18"|"$6]+=1 } END { for (category in sum) printf "%s\t%d\t%.2f\n", category, sum[category], (sum[category]/te_length[category])/(NR/1000000) } ' single_te.bed > ReadCountAndTPM
# RowName
awk -v OFS="|" '$1!~/#/{ gsub(/"Motif:|"/,"",$10); print $1,$4,$5,$10,$11,$12,$7 }' $te_gff > RowName
# merge
awk -v OFS="\t" '
  NR==FNR { k[$1]=$0; next }
    { if(k[$1]!="") {
      print k[$1]
      } else {
      print $1,0,0} 
    }'  ReadCountAndTPM RowName | sed "1i TE\tCount\tTPM" > ${output_prefix}_LATTE_TE_Expression

cat ./SingleTEReads > ${output_prefix}_LATTE_TE_Reads
# echo "."

# delect the temporary directory
cd ../
[ ${reserve_tmp_dir} == "True" ] || rm -rf ${tmp_dir}

