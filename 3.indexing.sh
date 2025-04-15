#!/bin/sh

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference. If you have a large genome, this will req$
##              First need to use gffread to convert annotation file from .gff3 to .gft format.
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing   InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files
## Samtools  Convert .sam to .bam and sort          Input: Alignment files,  .sam
##                                                  Output: Sorted  .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]",
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
## queue: bigmem
## core: 6
## time limit (HH:MM:SS): 360:00:00  
## Memory: 250gb 
###############################################

#### Load all the programs you are going to use in this script.
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load gffread



#  Set the stack size to unlimited
ulimit -s unlimited
# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## Replace the [#] with paths to define these variable
MyID=aubclsd0335          ## Example: MyID=aubtss

WD=/scratch/$MyID/RNAProject_dog/                      ## Example:/scratch/$MyID/PracticeRNAseq
CD=/scratch/$MyID/RNAProject_dog/CleanData             ## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files a$
REFD=/scratch/$MyID/RNAProject_dog/DogRefGenome    ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed re$
MAPD=/scratch/$MyID/RNAProject_dog/Map_HiSat2           ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/scratch/$MyID/RNAProject_dog/Counts_StringTie       ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/RNAProject_dog/Counts_H_S          ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S

REF=genomic                   ## This is what the "easy name" will be for the genome reference
## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
cp /home/${MyID}/RNAProject_dog/ncbi_dataset/data/GCF_000002285.5/${REF}.fna .
cp /home/${MyID}/RNAProject_dog/ncbi_dataset/data/GCF_000002285.5/${REF}.gff .

###  Identify exons and splice sites on the reference genome
gffread ${REF}.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Tasha_index


