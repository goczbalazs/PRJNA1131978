#!/bin/bash
#####FastQC analysis#####
export PATH=$PATH:/mnt/c/Users/gocz.balazs/Documents/tools/FastQC

# Creating the folder for the QC results
mkdir qc

# Running fastqc
for f in *fastq 
do 
  fastqc -t 12 $f -o qc
done 

######First round of adapter trimming with trimmomatic####
TRIM='/mnt/c/Users/gocz.balazs/Documents/tools/Trimmomatic-0.39/trimmomatic-0.39.jar'

for f in *.fastq
do 
  o=${f/'.fastq'/'_trimmed.fastq'}
  log=${f/'.fastq'/'_trim.out.log'}
  java -jar $TRIM SE -threads 12 $f $o LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> $log
  echo complete ${f}
done

#####Second round of adapter trimming with Cutadapt#####
for f in *_trimmed.fastq
do
log=${f/'.fastq'/'_trim.out.log'}
    cutadapt -j 12 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --minimum-length 36 -o ${f/'_trimmed.fastq'/'_trimmed_final.fastq'} $f > $log
	echo complete ${f}
done

#####Aligning with STAR######
for f in *_trimmed_final.fastq
do
  root=${f/'_final.fastq'/''} 
STAR -- genomeDir /mnt/c/Users/gocz.balazs/Documents/tools/ref_mouse_111 -- readFilesIn $f -- outFileNamePrefix 'FIRST_GRCm39_111_'$root -- outReadsUnmapped Fastx -- outSAMtype BAM SortedByCoordinate -- outMultimapperOrder Random -- genomeLoad LoadAndKeep --limitBAMsortRAM 40000000000 -- runThreadN 12
done

for f in *_trimmed_final.fastq
do
  root=${f/'_final.fastq'/''} 
STAR -- genomeDir /mnt/c/Users/gocz.balazs/Documents/tools/ref_mouse_111 -- readFilesIn $f -- outFileNamePrefix 'Remove_'$root -- outReadsUnmapped Fastx -- outSAMtype BAM SortedByCoordinate -- outMultimapperOrder Random -- runThreadN 12 -- genomeLoad Remove --limitBAMsortRAM 40000000000
done

for f in *_trimmed_final.fastq
do
  root=${f/'.fastq'/''}
  STAR -- genomeDir /mnt/c/Users/gocz.balazs/Documents/tools/ref_mouse_111 -- sjdbFileChrStartEnd *SJ.out.tab -- readFilesIn $f -- outFileNamePrefix 'final_bams/Final_GRCm39_111_'$root -- outReadsUnmapped Fastx -- outSAMtype BAM SortedByCoordinate -- outMultimapperOrder Random -- runThreadN 12
done

cd final_bams

####Assigned with FeatureCounts#####
mkdir feature_counts_2_strand_no_multimappers

featureCounts -s 2 -T 12 -a /mnt/c/Users/gocz.balazs/Documents/tools/Mus_musculus.GRCm39.111.gtf -o feature_counts_2_strand_no_multimappers/featureCounts_GRCm39.111_Mouse_GnRH.txt Final_GRCm39_111*bam

annotationGTF='/mnt/c/Users/gocz.balazs/Documents/tools/Mus_musculus.GRCm39.111_COLLAPSED.gtf'

#####Quality checking with RSeQC#####
for b in *_trimmed_finalAligned.sortedByCoord.out.bam
do
  s=${b/'_trimmed_finalAligned.sortedByCoord.out.bam'/''}
  s=${s/'GRCm39_111_'/''}
 /mnt/c/Users/gocz.balazs/Documents/tools/rnaseqc $annotationGTF $b rnaseqcRes --sample=$s --unpaired
done