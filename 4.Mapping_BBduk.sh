#!/bin/sh


#### Load all the programs you are going to use in this script.
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread
#module load gffcompare


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
CD=/scratch/$MyID/RNAProject_dog/CleanBBduk             ## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files a$
REFD=/scratch/$MyID/RNAProject_dog/DogRefGenome    ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed re$
MAPD=/scratch/$MyID/RNAProject_dog/BBdukMap_HiSat2           ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/scratch/$MyID/RNAProject_dog/BBdukCounts_StringTie       ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/RNAProject_dog/BBdukCounts_H_S          ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S

REF=genomic                   ## This is what the "easy name" will be for the genome reference
## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD
##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
#cp /home/${MyID}/RNAProject_dog/ncbi_dataset/data/GCF_000002285.5/${REF}.fna .
#cp /home/${MyID}/RNAProject_dog/ncbi_dataset/data/GCF_000002285.5/${REF}.gff .

###  Identify exons and splice sites on the reference genome
#gffread ${REF}.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
#hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
#hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
#hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna TashaDog_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd ${CD}  #### This is where our clean paired reads are located.
## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${CD}/list  .

## process the samples in the list, one by one using a while loop
while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "${REFD}"/Tasha_index       \
    -1 "${CD}"/"$i"_1_paired.fastq  -2 "${CD}"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam
    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model.
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons

mkdir "${COUNTSD}"/"$i"
stringtie -p 6 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"_sorted.bam

done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt ${RESULTSD}

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix.
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory.
cp *.csv ${RESULTSD}
## move these results files to your personal computer for downstream statistical analyses in R studio.



