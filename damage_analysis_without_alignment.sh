#!/bin/bash

#Assess fragment length of reads and calculate proportion of each base at both ends of fastq reads

#Specify variables
ROOTDIR=/home/a1698312
TRIMDIR=$ROOTDIR/weyrich/trimData
LOGFILE=$TRIMDIR/damageAnalysisLog.txt

################ Count read lengths ######################

#Change into directory where trimmed_fastq files located
if [ -d ${TRIMDIR} ]
then
echo "Changing to trimData directory"
cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

##### create log file for script ####

if [ ! -f damageAnalysisLog.txt ]
then
  echo -e 'Creating file damageAnalysisLog.txt'
  echo -e "file\ttime" > ${LOGFILE}
else
  echo -e 'damageAnalysisLog.txt already exists'
fi

#Generate text file for storing length count data
if [ ! -f fastq_length.txt ]
then
  echo -e 'Creating file fastq_length.txt'
  echo -e 'occ\tlength' > fastq_length.txt
else
  echo  'fastq_length.txt already exists'
fi

#Generate text file for storing base proportions
if [ ! -f base_proportions.txt ]
then
echo -e 'Creating file base_proportions.txt'
echo -e "fileName\tend\tpos\tA\tC\tG\tT" > base_proportions.txt
else
  echo  'base_proportions.txt already exists'
fi

#Count fastq lengths and proportion of each base for each file

for fastq_file in *fastq.gz
  do
    #variable to start timing loop
    STARTTIME=$(date +%s)
    echo -e "Collating read lengths in ${fastq_file}"

#extract DNA sequence from file, count length of each sequence, sort lengths and count number of unique
    LENGTH=$(zcat ${fastq_file} | sed -n '2~4p' | awk '{print length($1)}' | sort -n | uniq -c)

    #print fileName and collated lengths to bottom of fastq_length.txt
    echo -e "${fastq_file}"  >> fastq_length.txt
    echo -e "${LENGTH}"  >> fastq_length.txt

    #variable to count number of reads in file
#    FASTQcount=$(zcat ${fastq_file} | sed -n '1~4p' | wc -l)

    #start counter
    n=1
    echo -e "Collating base proportions in ${fastq_file}"

    #loop which increases value of n each time up to 25
    while [ $n -le 25 ]
    do

      #extract DNA sequence from file, print base at position n in sequence, count number of each base
      COUNT5p=$(zcat ${fastq_file} | sed -n '2~4p' |  awk 'length($0) > '"$(($n*2))"'' | \
        gawk -F '' '{print $'"$n"'}' | sort | uniq -c | awk '{print}' ORS='\t ')

      echo -e "${fastq_file}\t5p\t${n}\t${COUNT5p}" >> base_proportions.txt

      #extract DNA sequences, reverse them (so now reading from 3' end of read), print base at position n, count number of each base
      COUNT3p=$(zcat ${fastq_file} | sed -n '2~4p' | awk 'length($0) > '"$(($n*2))"'' | rev | \
        gawk -F '' '{print $'"$n"'}' | sort | uniq -c | awk '{print}' ORS='\t ')

      echo -e "${fastq_file}\t3p\t${n}\t${COUNT3p}" >> base_proportions.txt

      #add 1 to the counter
      let n=n+1

      #end inner/while loop
    done

    #calculate time taken to complete analysis on this file
    ENDTIME=$(date +%s)
    echo -e "${fastq_file}\t$(($ENDTIME - $STARTTIME))" >> ${LOGFILE}

  #end outer/for loop
  done