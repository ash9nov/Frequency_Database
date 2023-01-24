#!/bin/bash
#
# input args = 1)INPUTFOLDER where 2)normalized vcf file (not zipped) is located

INFOLDER=$1
COHORT=$2

DATE=`date +'%y-%m-%d'`
mkdir -p $INFOLDER/inhousedb_$DATE

python indbcreator_stolav_normalized.py $INFOLDER $COHORT $DATE

#bgzip
bgzip -c $INFOLDER/inhousedb_$DATE/inhousedb.vcf > $INFOLDER/inhousedb_$DATE/inhousedb.vcf.gz

#tabix of bgzipped file
tabix $INFOLDER/inhousedb_$DATE/inhousedb.vcf.gz

## Upload to TSD through tacl p874 --upload $INFOLDER/inhousedb_$DATE/