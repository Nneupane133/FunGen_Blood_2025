#! /bin/bash

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsd0335         ## Example: MyID=aubtss

  ## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/${MyID}/RNAProject_dog/RawData/                     ## Example: DD=/scratch/${MyID}/PracticeRNAseq/RawData
WD=/scratch/${MyID}/RNAProject_dog                          ## Example: WD=/scratch/${MyID}/PracticeRNAseq
RDQdog=RawDataQualitydog

 
##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p ${DD}
## move to the Data Directory
cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
	## -F 	Defline contains only original sequence name.
	## -I 	Append read id after spot id as 'accession.spot.readid' on defline.
	## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are from Bioproject PRJNA629466, PRJNA803741, PRJNA823683, PRJNA1011250, PRJNA1028815. 
## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
        ## then run the script by using "run_script [script name]"
        ## suggested paramenters are below to submit this script.
                ## queue: bigmem
                ## core: 6
                ## time limit (HH:MM:SS): 360:00:00  
                ## Memory: 200gb 
                ## run on asax
###############################################

vdb-config --interactive

# List of SRA accession numbers
sra_list=(SRR11650122
SRR11650127
SRR11650126
SRR11650131
SRR11650134
SRR11650138
SRR11650142
SRR18645116
SRR18645125
SRR17893388
SRR17893389
SRR18645117
SRR11650128
SRR11650129
SRR11650140
SRR11650141
SRR25832772
SRR25832774
SRR26403891
SRR26403890
SRR26403899
SRR17893386
SRR26403894
SRR26403914
SRR26403927
)

# Loop through each accession and download
for sra in "${sra_list[@]}"; do
    echo "Downloading $sra..."
    fastq-dump -F --split-files "$sra"
done

echo "All downloads complete."

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
mkdir ${WD}/${RDQdog}
fastqc *.fastq --outdir=${WD}/${RDQdog}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd ${WD}/${RDQdog}
tar cvzf ${RDQdog}.tar.gz  ${WD}/${RDQdog}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.

