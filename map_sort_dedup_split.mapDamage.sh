#!/bin/bash

#USAGE: Requires trimmed_fastq files in trimData directory and 2 text files: expandedGenomeList.txt & chrID.txt
#       Specify variable for location of Ref Seq genomes for downloading and alignment

## Input: Trimmed FASTQ file
## Output:
##  - Deduplicated bam (DELETE after split and count)
##  - Split bam by 30 genomes
##  - mapDamage on each split bam

# Batch directory variables
# Just to make sure we dont overwrite everything on our phoenix directory
#batchDir=$1

#Specify variables
ROOTDIR=/home/a1698312
TRIMDIR=$ROOTDIR/testPhoenix/trimData
MAPDIR=$ROOTDIR/testPhoenix/mapData
MAP_DAMAGE_DIR=$ROOTDIR/testPhoenix/mapDamageData
LOGFILE=$ROOTDIR/testPhoenix/mapDamageLog.txt
COUNTFILE=$ROOTDIR/testPhoenix/mapData/map_count.txt

##### create log file for script ####

if [ ! -f mapDamageLog.txt ]
then
echo -e 'Creating file mapDamageLog.txt'
echo -e "task\tfile\ttime(s)" > ${LOGFILE}
else
  echo -e 'mapDamageLog.txt already exists'
fi

#### create directories - mapData; mapDamageData ####

#Create mapDamageData directory
if [ ! -d ${MAP_DAMAGE_DIR} ]
then
echo "Creating ${MAP_DAMAGE_DIR}"
mkdir -p ~/simData/mapDamageData
else
  echo "${MAP_DAMAGE_DIR} already exists"
fi

##### bwa_build alignment index ####

#Create mapData directory
if [ ! -d ${MAPDIR} ]
then
echo "Creating ${MAPDIR}"
mkdir -p mapData
echo "Changing into ${MAPDIR}"
cd ${MAPDIR}
else
  echo "${MAPDIR} already exists. Changing into ${MAPDIR}"
cd ${MAPDIR}
fi

#log time to build index (including downloading of fasta files)
STARTTIME=$(date +%s)

#Download fasta file for each genome in genomeList.txt file
while read -r line
do
link=$(echo "${line}" | cut -f8)
ref=$(echo "${line}" | cut -f1)
wget -c "${link}" -O "${ref}".fna.gz
done < ${ROOTDIR}/expandedGenomeList.txt

#unzip fasta files
#gunzip *fna.gz

#concatenate fasta files
cat *fna.gz > combined.fna.gz

unpigz combined.fna.gz

#build-index for alignment
bwa index -p bwaidx combined.fna

## Jimmy testing and suggesting (maybe test - but dont worry to much)
# bwa index -p bwaidx <(zcat combined.fna.gz)

#Add time stamp to logfile
ENDTIME=$(date +%s)
echo -e "build_index\tNA\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}

################ BWA Alignment #################

#Change into directory where trimmed_fastq files located
if [ -d ${TRIMDIR} ]
then
echo "Changing to trimData directory"
cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

#bwa alignment of collapsed reads
for fastq_file in *fastq.gz
do
STARTTIME=$(date +%s)
echo "Aligning ${fastq_file}"
bwa aln -n 0.01 -o 2 -l 1024 -t 4  $MAPDIR/bwaidx $fastq_file > ${fastq_file/%.fastq.gz/_MAPPED.sai}
ENDTIME=$(date +%s)
echo -e "bwa_alignment\t${fastq_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done

#Convert .sai alignment file to bam format with the sam header. Use -q 25. Exclude unmapped reads.
for map_file in *_MAPPED.sai
do
echo "Converting ${map_file} to bam format"
PREFIX=${map_file%%_MAPPED.sai}
bwa samse $MAPDIR/bwaidx \
${PREFIX}_MAPPED.sai \
${PREFIX}.fastq.gz | \
samtools view -q 25 -bSh -F0x4 -> $MAPDIR/${PREFIX}_bwa.bam
done


#Remove .sai files as no longer needed
rm *_MAPPED.sai

################### sambamba sort and rmdup #####################

#Change into directory where mapped files located
if [ -d ${MAPDIR} ]
then
echo "Changing to mapData directory"
cd ${MAPDIR}
else
  echo "Cannot find ${MAPDIR}"
exit1
fi

#log time to sort and deduplicate .bam files
STARTTIME=$(date +%s)

#for bam_file in *_bwa.bam
do
PREFIX2=${bam_file%%_bwa.bam}
echo "Sorting bam file for ${bam_file}"
sambamba sort -o ${PREFIX2}_sorted.bam ${bam_file}
done

for sort_file in *_sorted.bam
do
PREFIX3=${sort_file%%_sorted.bam}
echo "Removing duplicates ${sort_file}"
sambamba markdup -r ${sort_file} ${PREFIX3}_rmdup.bam
done

#Remove _sorted.bam.bai files as no longer needed
rm *_sorted.bam.bai
#Remove _sorted.bam files as no longer needed
rm *_sorted.bam

#Add time stamp to logfile
ENDTIME=$(date +%s)
echo -e "sort_dedupe\tNA\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}

################# split bam file ####################

#use samtools view & chromosome ID to split into separate bam files
for rmdup_file in *_rmdup.bam
do
STARTTIME=$(date +%s)
PREFIX4=${rmdup_file%%_rmdup.bam}
echo -e "Splitting ${rmdup_file}"
# index combined.fna files
samtools faidx combined.fna
samtools view -q 25 ${rmdup_file} | awk '$3 ~ /NZ_LWMU/ { print $0 }' | samtools view -bt combined.fna.fai > ${PREFIX4}_M.oralis_split.bam
while read -r line; do
chrID=$(echo "${line}" | cut -f2)
ref=$(echo "${line}" | cut -f1)
samtools view -q 25 -bSh ${rmdup_file} ${chrID} > ${PREFIX4}_${ref}_split.bam
done < ${ROOTDIR}/chrID.txt
ENDTIME=$(date +%s)
echo -e "split_file\t${rmdup_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done

#remove *split files with fewer than 1000 lines
for split_file in *_split.bam
do
value=$( samtools view -c $split_file)
if [ $value -lt 500 ]
then
rm $split_file
fi
done

################# count reads in each alignment file #################

#Generate text file for storing alignment count data
if [ ! -f ${COUNTFILE} ]
then
echo -e "Creating file ${COUNTFILE}"
echo -e "count\tMAPQ" > ${COUNTFILE}
else
  echo  'map_count file already exists'
fi

#Timing variable
STARTTIME=$(date +%s)

#Count total number of reads in each bwa.bam, rmdup.bam, and split.bam file
for bam_file in *.bam
do
echo "Counting reads in ${bam_file}"
MAPCOUNT=$(samtools view ${bam_file} | cut -f5 | sort | uniq -c)
echo -e "${bam_file}" >> ${COUNTFILE}
echo -e "${MAPCOUNT}" >> ${COUNTFILE}
done

#calculate time to perform map_count
ENDTIME=$(date +%s)
echo -e "map_count\tNA\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}

################# mapDamage #########################

#Run mapDamage
for split_file in *_split.bam
do
STARTTIME=$(date +%s)
echo "Running mapDamage on ${split_file}"
mapDamage -d ${MAP_DAMAGE_DIR}/${split_file%%_split.bam} -i ${split_file} -r ${MAPDIR}/combined.fna
ENDTIME=$(date +%s)
echo -e "map_Damage\t${split_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done
