#!/bin/bash

#USAGE: Requires trimmed_fastq files in specified directory
#       Specify variable for location of Ref Seq genomes for downloading and alignment


#Specify variables
ROOTDIR=/home/a1698312/simData
TRIMDIR=$ROOTDIR/data
MAPDIR=$ROOTDIR/bwaMapData
BWAMEM_MAPDIR=$ROOTDIR/bwa_memMapData
bt2_MAPDIR=$ROOTDIR/bt2MapData
LOGFILE=$ROOTDIR/mapDamageLog.txt

##### create log file for script ####

if [ ! -f mapDamageLog.txt ]
then
echo -e 'Creating file mapDamageLog.txt'
echo -e "task\tfile\ttime(s)" > ${LOGFILE}
else
  echo -e 'mapDamageLog.txt already exists'
fi

##### bwa_build alignment index ####

#Change into bwaMapData directory
if [ ! -d ${MAPDIR} ]
then
  echo "Creating ${MAPDIR}"
  mkdir -p bwaMapData
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
  wget -c "$link" -O "${ref}".fna.gz
done < ${ROOTDIR}/simData/expandedGenomeList.txt

#unzip fasta files
gunzip *fna.gz

#concatenate fasta files
cat *fna > combined.fna

#build-index for alignment
bwa index -p bwaidx combined.fna

#Add time stamp to logfile
ENDTIME=$(date +%s)
echo -e "build_index\tNA\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}

################ BWA Alignment #################

#Change into directory where trimmed_fastq files located
if [ -d ${TRIMDIR} ]
then
  echo "Changing to data directory"
  cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

#bwa alignment of collapsed reads
for sim_file in *_endo*fa.gz
do
  STARTTIME=$(date +%s)
  echo "Aligning ${sim_file}"
  bwa aln -n 0.01 -o 2 -l 1024 $MAPDIR/bwaidx $sim_file -0 -t 4 > ${sim_file/%.fa.gz/_MAPPED.sai}
  ENDTIME=$(date +%s)
  echo -e "bwa_alignment\t${sim_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done

#Convert .sai alignment file to bam format with the sam header. Exclude unmapped reads.
for map_file in *_MAPPED.sai
  do 
    echo "Converting ${map_file} to bam format"
    PREFIX=${map_file%%_MAPPED.sai}
    bwa samse $MAPDIR/bwaidx \
              ${PREFIX}_MAPPED.sai \
              ${PREFIX}.fa.gz | \
                samtools view -bSh -F0x4 -> $MAPDIR/${PREFIX}_bwa.bam
  done


#Remove .sai files as no longer needed
rm *_MAPPED.sai

#Change to root directory
cd ${ROOTDIR}

################### BWA-MEM Alignment ###########################

#Create directory for storing bwa_memMapData
if [ ! -d ${BWAMEM_MAPDIR} ]
then
  echo "Creating ${BWAMEM_MAPDIR}"
  mkdir -p ~/simData/bwa_memMapData
else
  echo "${BWAMEM_MAPDIR} already exists"
fi

#Change into directory where trimmed_fastq files located
if [ -d ${TRIMDIR} ]
then
  echo "Changing to data directory"
  cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

#perform bwa-mem alignment on all simulated files
for sim_file in *fa.gz
do
  STARTTIME=$(date +%s)
  echo "Aligning ${sim_file} with bwa-mem"
  bwa mem -t 4 $MAPDIR/bwaidx $sim_file | samtools view -bSh -F0x4 -> $BWAMEM_MAPDIR/${sim_file/%.fa.gz/_bwamem.bam}
  ENDTIME=$(date +%s)
  echo -e "bwa_mem_alignment\t${sim_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done

#Change back to root directory
cd ${ROOTDIR}

##### build bt2 alignment index ####

#Change into bt2MapData directory
if [ ! -d ${bt2_MAPDIR} ]
then
  echo "Creating ${bt2_MAPDIR}"
  mkdir -p bt2MapData
  echo "Changing into ${bt2_MAPDIR}"
  cd ${bt2_MAPDIR}
else
  echo "${bt2_MAPDIR} already exists. Changing into ${bt2_MAPDIR}"
  cd ${bt2_MAPDIR}
fi

#log time to build index (including downloading of fasta files)
STARTTIME=$(date +%s)

#build-index for bowtie2 alignment
bowtie2-build -f ~/simData/bwaMapData/combined.fna ~/simData/bt2MapData/bt2idx

#Add time stamp to logfile
ENDTIME=$(date +%s)
echo -e "build_bt2_index\tNA\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}


################ bowtie2 Alignment - default parameters #################

#Change into directory where trimmed_fastq files located
if [ -d ${TRIMDIR} ]
then
  echo "Changing to data directory"
  cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

#bt2 alignment of collapsed reads
gunzip *fa.gz

#for sim_file in *fa
do
  STARTTIME=$(date +%s)
  echo "Aligning ${sim_file}"
  bowtie2 -f -p 4 -x $bt2_MAPDIR/bt2idx $sim_file | samtools view -bS -F0x4 - > ${bt2_MAPDIR}/${sim_file/%.fa/_bt2.bam}
  gzip $sim_file
  ENDTIME=$(date +%s)
  echo -e "bt2_alignment\t${sim_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}
done

