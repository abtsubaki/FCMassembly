#!/bin/bash
#PBS -N JOBNAME
#PBS -e jobname.err
#PBS -o jobname.out
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -l mem=5gb
#PBS -m abe

cd $PBS_O_WORKDIR
#####################################################################################################
###A bash script using several different modules for the hybrid assembly of the false codling moth genome 
###Assembly is done using PacBio long reads and Illumina short reads
###Assembly is based on https://www.melbournebioinformatics.org.au/tutorials/tutorials/pacbio/
#####################################################################################################

#set directories
DIR=$HOME/

#####################################################################################################
#QC of PacBio and Illumina reads
#fastQC
module load app/fastqc

echo "starting fastqc on Illumina reads\n"
gunzip -c $DIR/forwardreads.fastq.gz | fastqc -o $DIR
gunzip -c $DIR/reversereads.fastq.gz | fastqc -o $DIR
echo "finished Illumina reads\n"

echo "starting fastqc on PacBio reads"\n
gunzip -c $DIR/pacbioreads.fastq.gz | fastqc -o $DIR
echo "completed fastqc on PacBio reads"\n
#####################################################################################################
#Trim low quality bases
#FastXToolkit
module load app/fastx_toolkit

echo "starting trimreads on Illumina forward reads\n"
gunzip -c $DIR/forwardreads.fastq.gz | fastx_trimmer -f 1 -l 104 -o $DIR/forwardreads.trim.fastq
echo "finsihed trimreads on Illumnia forward reads\n"

echo "starting trimreads on Illumina reverse reads\n"
gunzip -c $DIR/reversereads.fastq.gz | fastx_trimmer -f 1 -l 89 -o $DIR/reversereads.trim.fastq
echo "finished trimreads on Illumina reverse reads\n"
#####################################################################################################
#Trim adaptors from Illumina reads
#Trimmomatic
module load app/Trimmomatic

echo "starting trim adaptors on R1\n"
trimmomatic PE -phred33 $DIR/forwardreads.trim.fastq $DIR/reversereads.trim.fastq \
$DIR/forwardreads.trimadapt.FP.fastq $DIR/forwardreads.trimadapt.FU.fastq \
$DIR/reversereads.trimadapt.RP.fastq $DIR/reversereads.trimadapt.RU.fastq \
ILLUMINACLIP:$DIR/adaptors.fa:2:30:10 MINLEN:50
echo "finished trimadaptors on R2\n"
#####################################################################################################
#Filter out contamination ~ remove potential human contaminants by mapping to human genome v38
module load app/bwa
module load app/samtools

#index human reference genome
echo "starting indexing\n"
bwa index -p href38 -a bwtsw $DIR/GRCh38.fasta
echo "completed indexing\n"

#Map 
echo "now running BWA mem\n"
bwa mem -t 8 $DIR/href38 $DIR/forwardreads.trimadapt.FP.fastq $DIR/reversereads.trimadapt.RP.fastq \
| samtools view -bS -h -o $DIR/fcmreads_trim_href_map.bam
echo "complete BWA mem\n"

#flagstat
echo "starting samtools flagstat\n"
samtools flagstat $DIR/fcmreads_trim_href_map.bam > $DIR/fcmreads_trim_href_map.flagstat
echo "completed flagstat\n"

#####################################################################################################
#ASSEMBLY
#kmer count with kmergenie
module load app/kmergenie

echo "starting kmergenie\n"
kmergenie $DIR/forwardreads.trimadapt.FP.fastq
kmergenie $DIR/reversereads.trimadapt.RP.fastq
echo "kmergenie complete\n"

echo "starting kmergenie\n"
gunzip $DIR/pacbioreads.fastq.gz
kmergenie $DIR/pacbioreads.fastq
echo "kmergenie complete\n"

#####################################################################################################
#assemble PacBioreads w/ Canu
module load app/canu/2.1.1

echo "starting assembly with Canu\n"
canu -p fcm -d $DIR/Canu1 genomeSize=0.7g gridOptions="-l walltime=300:00:00" -pacbio-hifi $DIR/pacbioreads.fastq
echo "completed assembly with Canu\n"
#####################################################################################################
#check assembly w/ Quast
module load app/QUAST/5.0

quast.py -o $DIR/Canu1_PacBio_Quast -e --large $DIR/Canu1/fcmpacbio.contigs.fasta
#####################################################################################################
#Align illumina reads to PacBio assembly
#index reference genome
echo "starting bwa index\n"
bwa index -p canu1 -a bwtsw $DIR/Canu1/fcmpacbio.contigs.fasta
echo "completed bwa index\n"

#Map
echo "now running BWA mem\n"
bwa mem -t 8 $DIR/canu1 $DIR/forwardreads.trimadapt.FP.fastq $DIR/reversereads.trimadapt.RP.fastq \
| samtools view -bS -h -o $DIR/fcmreads_canu1_map.bam
echo "complete BWA mem\n"

#sort bam file
echo "sorting bam file\n"
samtools sort $DIR/fcmreads_canu1_map.bam -o $DIR/fcmreads_canu1_map_sorted.bam
echo "done sorting bam file\n"

#index bam file
echo "indexing bam file\n"
samtools index $DIR/fcmreads_canu1_map_sorted.bam
echo "done indexing bam file\n"

#flagstat
echo "starting samtools flagstat\n"
samtools flagstat $DIR/fcmreads_canu1_map_sorted.bam > $DIR/fcmreads_canu1_map_sorted.flagstat
echo "completed flagstat\n"

#Check mapping quality for potential read discarding
echo 'this is bwa header\n'
samtools view -H $DIR/fcmreads_canu1_map_sorted.bam | tail -n 10

echo 'this is bwa alignments\n'
samtools view $DIR/fcmreads_canu1_map_sorted.bam | head -n 300

#Filter out mapped reads
echo "extract unmapped reads\n"
samtools view -h -f 4 -b $DIR/fcmreads_canu1_map_sorted.bam > $DIR/fcmreads_canu1_map_sorted.unmappedonly.bam
echo "done extracting unmapped reads\n"

#Check mapping quality for potential read discarding
echo 'this is bwa header\n'
samtools view -H $DIR/fcmreads_canu1_map_sorted.unmappedonly.bam | tail -n 10

echo 'this is bwa alignments\n'
samtools view $DIR/fcmreads_canu1_map_sorted.unmappedonly.bam | head -n 300

#convert bam file to fasta
echo "convert unmapped reads to fasta\n"
samtools fasta $DIR/fcmreads_canu1_map_sorted.unmappedonly.bam > $DIR/fcmreads_canu1_map_sorted.unmappedonly.fasta
echo "done converting to fasta\n"

#####################################################################################################
#assemble unmapped Illumina reads w/ Spades
module load app/SPAdes/3.13.0

echo "starting assembly with SPAdes\n"
spades.py --only-assembler --iontorrent -s $DIR/fcmreads_canu1_map_sorted.unmappedonly.bam -k 19,21 -o $DIR/Spades_Ill_Unmapd
echo "completed assembly with SPAdes\n"

#####################################################################################################
#Check assembly w/ Quast
echo "checking contigs with quast\n"
quast.py -o $DIR/Spades_Ill_Unmapd_contigs_Quast/ -e $DIR/Spades_Ill_Unmapd/contigs.fasta
echo "checking scaffolds with Quast\n"
quast.py -o $DIR/Spades_Ill_Unmapd_scaff_Quast/ -e $DIR/Spades_Ill_Unmapd/scaffolds.fasta
echo "done checking with quast\n"

#####################################################################################################
#Concatenate assembly 1 and 2
echo "start concatenation\n"
cat $DIR/Canu1/fcmpacbio.contigs.fasta $DIR/Spades_Ill_Unmapd/contigs.fasta > $DIR/fcmgenome_v1.fasta
echo "done concatenating\n"
#####################################################################################################
#Correct the assembly w/ Pilon
module load app/pilon

#MAPPING Illumina reads to fcmgenome_v1
#index reference genome
echo "starting bwa index\n"
bwa index -p fcmgenome1 -a bwtsw $DIR/fcmgenome_v1.fasta
echo "completed bwa index\n"

echo "starting fasta index\n"
samtools faidx $INDIR/fcmgenome_v1.fasta
echo "completed fasta index\n"

#Map
echo "now running BWA mem\n"
bwa mem -t 8 $DIR/fcmgenome1 $DIR/forwardreads.trimadapt.FP.fastq $DIR/reversereads.trimadapt.RP.fastq \
| samtools view -bS -h -o $DIR/fcmreads_fcmgenome1_map.bam
echo "complete BWA mem\n"

#sort bam file
echo "sorting bam file\n"
samtools sort $DIR/fcmreads_fcmgenome1_map.bam -o $DIR/fcmreads_fcmgenome1_map_sorted.bam
echo "done sorting bam file\n"

#index bam file
echo "indexing bam file\n"
samtools index $DIR/fcmreads_fcmgenome1_map_sorted.bam
echo "done indexing bam file\n"

#flagstat
echo "starting samtools flagstat\n"
samtools flagstat $DIR/fcmreads_fcmgenome1_map_sorted.bam > $DIR/fcmreads_fcmgenome1_map_sorted.flagstat
echo "completed flagstat\n"

#Check mapping quality for potential read discarding
echo 'this is bwa header\n'
samtools view -H $DIR/fcmreads_fcmgenome1_map_sorted.bam | tail -n 10

echo 'this is bwa alignments\n'
samtools view $DIR/fcmreads_fcmgenome1_map_sorted.bam | head -n 300

echo "starting pilon\n"
java -Xmx200g -jar /apps/pilon/pilon-1.23.jar --genome $DIR/fcmgenome_v1.fasta --frags $DIR/fcmreads_fcmgenome1_map.bam \
--output pilon1 --fix all --mindepth 0.5 --changes --verbose --outdir $DIR/Pilon --vcf
echo "pilon completed\n"
#####################################################################################################
#Deduplicate assembly with Purge Haplotigs
module load app/purge_haplotigs/current

#"First align your genome against the raw reads to produce a bam file, and use bam file below"

#ref="Reference Assembly"
#read_file="long read sequencing data"

##STEP-1
purge_haplotigs hist -b $DIR/fcmgenomev1_PacBioreads_bwa_sorted.bam -g $DIR/fcmgenome_v1.fasta -t 16


## STEP-2 Only run step two after checking the histogram.. Refer to paper or github for example
purge_haplotigs cov -i $DIR/fcmgenomev1_PacBioreads_bwa_sorted.bam.gencov -l 1 -m 6 -h 110

## STEP-3
purge_haplotigs purge -b $DIR/fcmgenomev1_PacBioreads_bwa_sorted.bam -g $DIR/fcmgenome_v1.fasta -c coverage_stats.csv -o "fcmgenomev2_purge" -d -a 55

######################################################################################################

#Evaluate assembly w/ BUSCO
module load app/BUSCO/5.2


echo "starting BUSCO\n"
busco -i $DIR/fcmgenomev2_purge.fasta -m genome -l $DIR/lepidoptera_odb10 -o fcmgenome_Purged_Busco
echo "completed BUSCO\n"

echo "plotting BUSCO results"\n
generate_plot.py -wd $DIR
echo "done plotting"\n
#####################################################################################################


