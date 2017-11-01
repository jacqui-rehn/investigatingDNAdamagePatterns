#!/bin/bash

#Collate the number of TP and FP for each genome in each .bam file

#Specify variables
ROOTDIR=/home/a1698312/simData
BWA_MAPDIR=$ROOTDIR/bwaMapData
BWAMEM_MAPDIR=$ROOTDIR/bwa_memMapData
bt2_MAPDIR=$ROOTDIR/bt2MapData
COUNT_FILE=$ROOTDIR/map_TP_FP_counts.txt

############## Create count file for storing data ###############

#Generate text file for storing alignment counts
if [ ! -f ${COUNT_FILE} ]
then
  echo -e "Creating ${COUNT_FILE}"
  echo -e "count\trefID\tMAPQ" > ${COUNT_FILE}
else
  echo  "${COUNT_FILE} already exists"
fi

######### Counts for bwa.bam files ############

#Change into directory where *bwa.bam files located
if [ -d ${BWA_MAPDIR} ]
then
  echo "Changing to bwaMapData directory"
  cd ${BWA_MAPDIR}
else
  echo "Cannot find ${BWA_MAPDIR}"
exit1
fi

#Count TP and FP for each *_bwa.bam file
for bwa_file in *_bwa.bam
do
  echo -e "Collating TP and FP for $bwa_file"
  #Collate number of TP (aligned to correct genome) without quality filtering
  TPcount=$(samtools view ${bwa_file} | awk '$1 ~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bwa_file}_TP" >> ${COUNT_FILE}
  echo -e "${TPcount}" >> ${COUNT_FILE}
  #count number of FP (aligned to incorrect genome) without quality filtering
  FPcount=$(samtools view ${bwa_file} | awk '$1 !~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bwa_file}_FP" >> ${COUNT_FILE}
  echo -e "${FPcount}" >> ${COUNT_FILE}
done

######### Counts for bwa_mem.bam files ############

#Change into directory where *bwa.bam files located
if [ -d ${BWAMEM_MAPDIR} ]
then
  echo "Changing to bwa_memMapData directory"
  cd ${BWAMEM_MAPDIR}
else
  echo "Cannot find ${BWAMEM_MAPDIR}"
exit1
fi

#Count TP and FP for each *_bwamem.bam file
for bwamem_file in *_bwamem.bam
do
  echo -e "Collating TP and FP for $bwamem_file"
  #Collate number of TP (aligned to correct genome) without quality filtering
  TPcount=$(samtools view ${bwamem_file} | awk '$1 ~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bwamem_file}_TP" >> ${COUNT_FILE}
  echo -e "${TPcount}" >> ${COUNT_FILE}
  #count number of FP (aligned to incorrect genome) without quality filtering
  FPcount=$(samtools view ${bwamem_file} | awk '$1 !~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bwamem_file}_FP" >> ${COUNT_FILE}
  echo -e "${FPcount}" >> ${COUNT_FILE}
done

######### Counts for bt2.bam files ############

#Change into directory where *bwa.bam files located
if [ -d ${bt2_MAPDIR} ]
then
  echo "Changing to bt2MapData directory"
  cd ${bt2_MAPDIR}
else
  echo "Cannot find ${bt2_MAPDIR}"
exit1
fi

#Count TP and FP for each *_bwamem.bam file
for bt2_file in *_bt2.bam
do
  echo -e "Collating TP and FP for $bt2_file"
  #Collate number of TP (aligned to correct genome) without quality filtering
  TPcount=$(samtools view ${bt2_file} | awk '$1 ~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bt2_file}_TP" >> ${COUNT_FILE}
  echo -e "${TPcount}" >> ${COUNT_FILE}
  #count number of FP (aligned to incorrect genome) without quality filtering
  FPcount=$(samtools view ${bt2_file} | awk '$1 !~ $3 {print $3, $5}' | sort | uniq -c)
  echo -e "${bt2_file}_FP" >> ${COUNT_FILE}
  echo -e "${FPcount}" >> ${COUNT_FILE}
done


