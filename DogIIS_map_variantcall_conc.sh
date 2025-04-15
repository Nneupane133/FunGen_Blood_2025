#! /bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to map RNAseq data to a set of CDS, call variants an take concensus 
##          sequence with reference replaced with new variants
##      Start with cleaned reads
## For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##  You may need to increase these for bigger datasets
## 		queue: medium
##		core: 6
##		time limit (HH:MM:SS): 6:00:00 
##		Memory: 6gb
##		
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1
module load multiqc
module load trimmomatic/0.39
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread


##########  Define variables and make directories
## !!!!!!!Replace with Your specific information!!!!!!!!

WD=/scratch/aubclsd0335/RNAProject_dog        	#Working directory  MAKE
#DD=$WD/RawData					# Raw Data downloaded  IF YOUR DATA IS ALREADY DOWNLOADED AND TRIMMED YOU DO NOT NEED THIS!
CD=$WD/CleanData            	## Cleaned data IF YOUR DATA ARE ALREADY CLEANNED THEN JUST POINT TO THAT DIRECTORY
REFD=$WD/Ref          			 # this directory contains the indexed IIS transcripts fasta file and where the HiSat2 index will be
REF=IIS_CDS                 ## This is what the "easy name" will be for the genome reference
MAPD=$WD/Map_HiSat2           	## Map_HiSat2      bams
REFDVAR=$WD/REF_SNP        # This is where the Samtools index of the reference file will be for calling SNPs
VAR=$WD/Variants				#  will hold the results files. the .vcf and .concensus fasta files and a copy of the sorted_MapOnly.bams.  
                        ## This folder will be tarballed and you can bring back to yoru computer.

## Make the directories and all subdirectories defined by the variables above
#mkdir -p $DD
mkdir -p $CD
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $REFDVAR
mkdir -p $VAR

vdb-config --interactive

############################***********  Mapping and Calling SNPs ************##########################
##################  Prepare the Reference Index of transcripts for mapping with HiSat2   #############################
cd ${REFD}
### Copy the reference set of cds (.fasta) to this REFD directory  **** You will need to replace my ASC ID (aubtss) for YOUR ASC ID
cp /home/aubclsd0335/class_shared/Dog_Tasha_GCF_000002285.5/data/GCF_000002285.5/${REF}.fasta   .

#### Create a HISAT2 index for the reference. NOTE every mapping program will need to build a its own index.
hisat2-build ${REF}.fasta IIS_CDS_index

###############  Prepare the index of the refence for variant calling with BCFtools.
cd ${REFDVAR}
### Copy the reference set of cds (.fasta) to this  directory  **** You will need to replace my ASC ID (aubtss) for YOUR ASC ID
cp /home/aubclsd0335/class_shared/Dog_Tasha_GCF_000002285.5/data/GCF_000002285.5/${REF}.fasta .

samtools faidx ${REF}.fasta

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
    -x "${REFD}"/IIS_CDS_index       \
    -l "${CD}"/"$i"_1_paired.fastq  -2 "${CD}"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  	## remove unmapped reads. Since most reads won't map to our selected transcripts, this makes a very small bam file that we can work with on our laptops
  samtools view -F 0x04 -b "$i"_sorted.bam > "$i"_sorted_mapOnly.bam


##############  Calling SNPs  #######################
cd ${MAPD}
###Call SNPs  and output a .vcf file
bcftools mpileup -f $REFDVAR/${REF}.fasta "$i"_sorted_mapOnly.bam | bcftools call -mv -Ov -o ${VAR}/"$i"_variants.vcf

### Call snps and generate the consensus sequence
bcftools mpileup -Ou -f $REFDVAR/${REF}.fasta "$i"_sorted_mapOnly.bam | bcftools call -c | /apps/x86-64/apps/bcftools_1.13/bin/vcfutils.pl vcf2fq > ${VAR}/"$i"_consensus.fastq

done<list


#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your SMALL BAM files - copy to the Sariants folder
cp ${MAPD}/*_sorted_mapOnly.bam ${VAR}

cd ${WD}
tar cvzf Variants.tar.gz Variants  ### This tarball you can bring back to your computer


