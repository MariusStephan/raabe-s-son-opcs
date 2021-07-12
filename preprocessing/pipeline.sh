#!/bin/bash

# Convert BCL files to FASTQ and demultiplex
bcl2fastq --create-fastq-for-index-reads --minimum-trimmed-read-length=8 --no-lane-splitting --mask-short-adapter-reads=8 --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -r 6 -w 6 -R ~/Project/BCL/201002_NS500378_0362_AHTKYTBGXG --sample-sheet=/media/ss-admin/Filmdatenbank/SampleSheet.csv --output-dir=~/Project/FASTQ

#cd ~/Project/FASTQ
for sample in NPC_S1 OPC1_S_S5 OPC2_SON_S6
do
cat ${sample}_L001_I1_001.fastq.gz ${sample}_L002_I1_001.fastq.gz ${sample}_L003_I1_001.fastq.gz ${sample}_L004_I1_001.fastq.gz >  ${sample}_I1_001.fastq.gz
cat ${sample}_L001_R1_001.fastq.gz ${sample}_L002_R1_001.fastq.gz ${sample}_L003_R1_001.fastq.gz ${sample}_L004_R1_001.fastq.gz >  ${sample}_R1_001.fastq.gz
cat ${sample}_L001_R2_001.fastq.gz ${sample}_L002_R2_001.fastq.gz ${sample}_L003_R2_001.fastq.gz ${sample}_L004_R2_001.fastq.gz >  ${sample}_R2_001.fastq.gz
done

# Generate Meta Reference bundles for DropSeq
## GRCh38.p13 and Ensembl Chr. Annotation Release 100 were used with the addition of the 3xFLAG and LTR sequence regions from lentiviral constructs used.
## Pipeline works with unmodified files perfectly well

mkdir ~/DropSeq_MetaBundles/Metabundle_Human_S_SON
cd ~/DropSeq_MetaBundles/Metabundle_Human_S_SON
~/Drop-seq_tools-2.3.0/create_Drop-seq_reference_metadata.sh -d ~/Drop-seq_tools-2.3.0/ -n hg38_S_SON_DropSeq_MetaData -r ~/Genomes/Ensembl/GRCh38_S_SON_Constructs.fa -s human -g ~/Genomes/Ensembl/GRCh38_S_SON_Constructs.gtf

cd ~/Project/

mkdir FASTQC
mkdir BAM
mkdir ddSeeker_Summary
mkdir TAGGED_FASTQ
mkdir STAR
mkdir MERGED_BAM
mkdir ANNOTATED_BAM
mkdir DropSeq_Reports
mkdir FINAL_BAM
mkdir DGE

for sample in NPC_S1 OPC1_S_S5 OPC2_SON_S6
do

### Sequence quality control
fastqc -t 4 -o FASTQC FASTQ/${sample}_R1_001.fastq.gz
fastqc -t 4 -o FASTQC FASTQ/${sample}_R2_001.fastq.gz

### Cell barcode tagging
~/ddSeeker/code/ddSeeker.py --input FASTQ/${sample}_R1_001.fastq.gz FASTQ/${sample}_R2_001.fastq.gz --output BAM/${sample}_tagged.bam --cores 4 -s ddSeeker_Summary/${sample}

### Cell barcode quality control/Knee plot
~/ddSeeker/code/make_graphs.R ddSeeker_Summary/${sample} 2000

### Sort reads by query name
cd ~/Project/BAM/
java -jar ~/Drop-seq_tools-2.3.0/3rdParty/picard/picard.jar SortSam SORT_ORDER=queryname I=${sample}_tagged.bam O=${sample}_tagged_qsorted.bam

### Filter reads for rejected cell barcodes
cd ~/Project/BAM/
~/Drop-seq_tools-2.3.0/FilterBam TAG_REJECT=XE I=${sample}_tagged_qsorted.bam O=${sample}_tagged_qsorted_filtered.bam
rm ${sample}_tagged_qsorted.bam


### Convert filtered BAM back to FASTQ for mapping
cd ~/Project/
java -Xmx500m -jar ~/Drop-seq_tools-2.3.0/3rdParty/picard/picard.jar SamToFastq I=BAM/${sample}_tagged_qsorted_filtered.bam FASTQ= TAGGED_FASTQ/${sample}_tagged_qsorted_filtered.fastq

### Map Read2 to genome with STAR aligner
cd ~/Project/
STAR --genomeDir ~/DropSeq_MetaBundles/Metabundle_Human_S_SON/STAR/ --runThreadN 4 --outFileNamePrefix STAR/${sample}_ --readFilesIn TAGGED_FASTQ/${sample}_tagged_qsorted_filtered.fastq
rm TAGGED_FASTQ/${sample}_tagged_qsorted_filtered.fastq

### Sort STAR output to be on the safe side
cd ~/Project/STAR
java -jar ~/Drop-seq_tools-2.3.0/3rdParty/picard/picard.jar SortSam SORT_ORDER=queryname I=${sample}_Aligned.out.sam O=${sample}_aligned_qsorted.bam
rm ${sample}_Aligned.out.sam
		
### Merge aligned with unaligned BAM
cd ~/Project/
java -Xmx4000m -jar ~/Drop-seq_tools-2.3.0/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=~/DropSeq_MetaBundles/Metabundle_Human_S_SON/hg38_S_SON_DropSeq_MetaData.fasta.gz UNMAPPED_BAM=BAM/${sample}_tagged_qsorted_filtered.bam ALIGNED_BAM=STAR/${sample}_aligned_qsorted.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false O=MERGED_BAM/${sample}_merged.bam
rm ${sample}_aligned_qsorted.bam

### Tag merged BAM with functional annotation
cd ~/Project/
~/Drop-seq_tools-2.3.0/TagReadWithGeneFunction I=MERGED_BAM/${sample}_merged.bam O=ANNOTATED_BAM/${sample}_merged_annotated.bam ANNOTATIONS_FILE=~/DropSeq_MetaBundles/Metabundle_Human_S_SON/hg38_S_SON_DropSeq_MetaData.refFlat

### Repair substitutions in cell barcodes
cd ~/Project/ANNOTATED_BAM
~/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors INPUT=${sample}_merged_annotated.bam OUTPUT=${sample}_merged_annotated_repaired.bam MIN_UMIS_PER_CELL=20 OUTPUT_REPORT=~/Project/DropSeq_Reports/${sample}_substitution_error_report.txt TMP_DIR=~/FASTQ
rm ${sample}_merged_annotated.bam

### Repair synthesis errors in cell barcodes
cd ~/Project/
~/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors INPUT=ANNOTATED_BAM/${sample}_merged_annotated_repaired.bam MIN_UMIS_PER_CELL=20 OUTPUT_STATS=DropSeq_Reports/${sample}_synthesis_error_stats.txt SUMMARY=DropSeq_Reports/${sample}_synthesis_error_summary.txt REPORT=DropSeq_Reports/${sample}_synthesis_error_report.txt CREATE_INDEX=true TMP_DIR=~/FASTQ OUTPUT=FINAL_BAM/${sample}_final.bam

echo "---------------------------"
echo "###########################"
echo ""
echo "$sample was repaired!"
echo ""
echo "###########################"
echo "---------------------------"

### UMI count and DGE on the 2000 most abundant cell barcodes
cd ~/Project/
~/Drop-seq_tools-2.3.0/DigitalExpression I=FINAL_BAM/${sample}_final.bam TMP_DIR=./ O=DGE/${sample}_DGE_2000cells.txt.gz SUMMARY=DropSeq_Reports/${sample}_DGE_2000cells.summary.txt NUM_CORE_BARCODES=2000

	echo "---------------------------"
	echo "###########################"
	echo ""
	echo "$sample is finished!"
	echo ""
	echo "###########################"
	echo "---------------------------"
done

cd ~/Project/
echo "---------------------------"
echo ""
echo "Done!"
echo ""
echo "---------------------------"

### DropEst 

mkdir DropEst
mkdir Velocyto

## DropEst DGE
~/dropEst/bin/dropest -m -V -b -f -g ~/DropSeq_MetaBundles/Metabundle_Human_S_SON/hg38_S_SON_DropSeq_MetaData.gtf -o DropEst/${sample} -L eiEIBA -c ~/dropEst/configs/drop_seq_velocyto.xml FINAL_BAM/${sample}_final.bam >> dropest.log 2>&1

## Correct mutated Cell Barcode (missed by Drop-Seq Pipeline?)
velocyto tools dropest-bc-correct DropEst/${sample}_final.tagged.bam DropEst/${sample}.rds

## Run velocyto
velocyto run -o Velocyto -e ${sample} DropEst/correct_${sample}_final.tagged.bam ~/DropSeq_MetaBundles/Metabundle_Human_S_SON/hg38_S_SON_DropSeq_MetaData.gtf >> velocyto.log 2>&1
