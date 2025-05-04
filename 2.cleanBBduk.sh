# Load required modules
source /apps/profiles/modules_asax.sh.dyn
module load bbmap
module load fastqc/0.10.1

# Define user-specific variables
MyID=aubclsd0335  # Replace with your actual ASC ID

# Define directories
WD=/scratch/$MyID/RNAProject_dog  # Working directory
DD=/scratch/$MyID/RNAProject_dog/RawData  # Raw data directory
CD=/scratch/$MyID/RNAProject_dog/CleanBBduk  # Cleaned data directory
PCQ=PostCleanQualityBBduk  # Quality control results directory
adapters=AdaptersToTrim_All.fa  # Adapter sequences file

# Create necessary directories
mkdir -p ${CD}
mkdir -p ${WD}/${PCQ}

# Move to raw data directory
cd ${DD}

# Generate a list of unique sample names
ls | grep ".fastq" | cut -d "_" -f 1 | sort | uniq > list

# Copy adapter file from shared directory (Modify path if needed)
cp /home/${MyID}/class_shared/AdaptersToTrim_All.fa .

# Loop through each sample and run BBDuk for adapter and quality trimming
while read i
 do
    bbduk.sh in1="$i"_1.fastq in2="$i"_2.fastq \
    out1=${CD}/"$i"_1_paired.fastq out2=${CD}/"$i"_2_paired.fastq \
    ref=AdaptersToTrim_All.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=rl trimq=30 minlen=36

    # Run FastQC on the cleaned reads
    fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
    fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}

done < list  # End of loop

# Move to FastQC results directory
cd ${WD}/${PCQ}

# Compress FastQC results for easy download
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*

# Instructions to download the results:
# Use 'scp' or 'rsync' to transfer ${PCQ}.tar.gz to your local machine for inspection.

