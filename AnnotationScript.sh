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
###A bash script using several different modules for the preliminary annotation of the false codling moth genome 
###Annotation is done using the Codling Moth genome, Cydia pomonella, GCA_003425675.2 Cpom.V2, Wan et al. 2019
###Annotation is performed in a comparative manner using the CDS sequences from Wan et al., 2019's paper and BLASTn
###Annotation is performed in a comparative manner using the OR and IR genes and associated SNPs from Wan et al., 2019's paper
#####################################################################################################

#set directories
DIR=$HOME/

#####################################################################################################################################################
#BLAST - first, use BLASTn to identify regions in the FCM assembly that match CDS regions from the Codling moth genome (Wan et al. 2019)
module load app/NCBI/2.4.0+

echo "starting to make custom BLAST database of Cpom CDS\n"
makeblastdb -in $DIR/CodlingMoth_CDS.fasta -dbtype nucl -out CodlingMoth_CDS_BLASTdb
echo "completed custom database\n"

echo "starting blastn\n"
blastn -db $DIR/CodlingMoth_CDS_BLASTdb -outfmt "7 qacc sacc evalue pident mismatch gapopen qstart qend sstart send length" -query $INDIR/fcmgenome_v2.fasta \
-out $OUTDIR/CpomCDSBLASTfcmgenome2 \
-num_threads $NP
echo "done blastn\n"

#Filter BLAST hits to Codling moth CDS for sequence identities > 90%, >85%, >80% and >75% respectively
grep -v '#' CpomCDSBLASTfcmgenome2 | sort -n -k 4 | awk '$4>"75"' | awk '$4<"80"' > CpomCDSBLASTfcmgenomev2_sim75
grep -v '#' CpomCDSBLASTfcmgenome2 | sort -n -k 4 | awk '$4>"80"' | awk '$4<"85"' > CpomCDSBLASTfcmgenomev2_sim80
grep -v '#' CpomCDSBLASTfcmgenome2 | sort -n -k 4 | awk '$4>"85"' | awk '$4<"90"' > CpomCDSBLASTfcmgenomev2_sim85
grep -v '#' CpomCDSBLASTfcmgenome2 | sort -n -k 4 | awk '$4>"90"' > CpomCDSBLASTfcmgenomev2_sim90

#Determine unique hits between the four similarity brackets of > 90%, >85%, >80% and >75%
cut -f 2 CpomCDSBLASTfcmgenomev2_sim75 | sort -u | wc -l)
cut -f 2 CpomCDSBLASTfcmgenomev2_sim80 | sort -u | wc -l)
cut -f 2 CpomCDSBLASTfcmgenomev2_sim85 | sort -u | wc -l)
cut -f 2 CpomCDSBLASTfcmgenomev2_sim90 | sort -u | wc -l)

#Make list of only query names for BLAST hits
cut -f 2 CpomCDSBLASTfcmgenomev2_sim75 > querylist_sim75.txt
cut -f 2 CpomCDSBLASTfcmgenomev2_sim80 > querylist_sim80.txt
cut -f 2 CpomCDSBLASTfcmgenomev2_sim85 > querylist_sim85.txt
cut -f 2 CpomCDSBLASTfcmgenomev2_sim90 > querylist_sim90.txt
 
#Obtain query sequences based on the query names
#use python script gistfile.py
module load python/3.5.1

python3 gistfile.py CodlingMoth_CDS.fasta querylist_sim75.txt > queryseqs_sim75.fa
python3 gistfile.py CodlingMoth_CDS.fasta querylist_sim80.txt > queryseqs_sim80.fa
python3 gistfile.py CodlingMoth_CDS.fasta querylist_sim85.txt > queryseqs_sim85.fa
python3 gistfile.py CodlingMoth_CDS.fasta querylist_sim90.txt > queryseqs_sim90.fa
  
#Obtain peptide sequences of queries
#use python script gistfile.py
python3 gistfile.py CodlingMoth_PEPTIDE.fasta querylist_sim75.txt > querypep_sim75.fa
python3 gistfile.py CodlingMoth_PEPTIDE.fasta querylist_sim80.txt > querypep_sim80.fa
python3 gistfile.py CodlingMoth_PEPTIDE.fasta querylist_sim85.txt > querypep_sim85.fa
python3 gistfile.py CodlingMoth_PEPTIDE.fasta querylist_sim90.txt > querypep_sim90.fa

#Extract only the unique identifiers across the filtered hits
#First, concatenate all the query names files for each filter bracket of > 90%, >85%, >80% and >75%
cat querylist_sim90.txt querylist_sim85.txt querylist_sim80.txt querylist_sim75.txt > querylistall.txt
#Extract unique lines only
sort -u querylistall.txt > querylist_uniq.txt
          
#get the Codling Moth CDS sequences (CodlingMoth_CDS.fasta) based on the unique querynames/unique BLAST hits only (querylist_uniq.txt)
#use python script mergefasta.py 
python3 mergefasta.py #outputs queryseqs_uniq.fa
#get the Codling moth peptide sequences (CodlingMoth_PEPTIDE.fasta) based on the querynames/unique BLAST hits only (querylist_uniq.txt)
python3 mergefasta.py #outputs querypep_uniq.fa
 
###########output two files containing Codling moth CDS nucleotide and peptide sequences that match our FCM assembly with >75% similarity############

#####################################################################################################################################################
#As the Codling Moth CDS BLAST hits do not have sequence identifiers yet, we need to do a general BLAST search to identify these sequences.
#To narrow the search down to only lepidopteran sequences, use -gilist lepidoptera.gi in your BLAST script

#!/bin/bash
#PBS -N Annot
#PBS -e annot.err
#PBS -o annot.out
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=300:00:00
#PBS -l mem=100gb
#PBS -m abe

NP=$(( $(cat $PBS_NODEFILE | wc -l) - 1 ))
cd $PBS_O_WORKDIR

module load app/NCBI/2.4.0+

DIR=$HOME/

echo "starting blastp\n"
blastp -db /apps/NCBI/BLAST/db/nr -outfmt "7 qacc sacc evalue pident mismatch gapopen qstart qend sstart send length" -query $INDIR/querypep_uniq.fa \
-gilist $REFDIR/lepidoptera.gi \
-out $OUTDIR/querypep_uniq_BLASTp_lepi \
-num_threads $NP
echo "done blastp\n"

echo "starting blastn\n"
blastn -db /apps/NCBI/BLAST/db/nt -outfmt "7 qacc sacc evalue pident mismatch gapopen qstart qend sstart send length" -query $INDIR/queryseqs_uniq.fa \
-gilist $REFDIR/lepidoptera.gi \
-out $OUTDIR/queryseqs_uniq_BLASTn_lepi \
-num_threads $NP
echo "done blastn\n"

####output two files containing Codling moth CDS nucleotide and peptide sequences with sequence identifiers that match our FCM assembly with >75% similarity######

#Each query has multiple BLAST hits to the NCBI database. We need to extract the best hit for each query based on evalue, %identity and length.
#First, count how many sequences we have: 
$grep -v '#' queryseqs_uniq_BLASTn_lepi | cut -f 1 | sort | uniq | wc -l
#next, extract all hits with a sequence similarity of 100% and count how many we have: 
$grep -v '#' queryseqs_uniq_BLASTn | awk '$4 > "99"' | wc -l
#869 sequences ~ that’s not many so let’s drop down our percentage criteria...At 96.9 we get 8247 sequences. Let’s work with that and ask if all 8090 Cpom identifiers are represented: 
$grep -v '#' queryseqs_uniq_BLASTn | awk '$4 > "96.9"' | cut -f 1 | sort | uniq | wc -l
#Let’s divide the datasets into similarity brackets of >99, 90, 80, 70: 
$grep -v '#' querypep_uniq_BLASTp_lepi | awk '$4 > "99"' > Filter/querypep_uniq_BLASTp_lepi_sim99
#Now extract the BLAST hit names and identify them.
cut -f 2 queryseqs_uniq_BLASTn_Sim99 > queryseqs_uniq_BLASTn_Sim99_hits 

#Use Entrez Batch search to retrieve sequence identifiers for the list. Save ‘complete record’ to a ‘file’ as a ‘summary:
#Retrieve records for 280 UID(s) (BLASTn)
#Retrieve records for 812 UID(s) (BLASTp)

#######################output file with the names/IDs of highly similar sequences that match the Codling moth CDS that match the FCM genome assembly###########################

#####################################################################################################################################################
#Wan et al. 2019 lists specific groups of chemosensory and insecticide resistance genes indetified in their Codling Moth genome. We will endeavour to find those same genes in our FCM assembly
#*OR: Chemosensory genes: Odorant receptors: 
#Cydia pomonella odorant receptors (44 available on GenBank)
#Download multifasta of all 44 sequences. = CpomORgenes.fasta

echo "starting to make custom BLAST database of fcm genome assembly\n"
makeblastdb -in $DIR/fcmgenomev2.fasta -dbtype nucl -out fcmv2
echo "completed custom database\n"

echo "starting blastn\n"
blastn -db $DIR/fcmv2 -outfmt "7 qacc sacc evalue pident mismatch gapopen qstart qend sstart send length" -query $INDIR/CpomORgenes.fasta \
-out $OUTDIR/fcmgenome2BLASTCpomORgenes \
-num_threads $NP
echo "done blastn\n"

#23 of the 44 OR genes had BLAST hits to the fcmgenomev2 assembly. Percentage identity ranged from 73% to 100% and e-values ranging from 0.0 to 9.92e-09.
##########################Now we have a dataset, of chemosensory genes in our assembly->CpomORgenes_fcmgenomev2_BLASTn###############################
      
#*IR: Insecticide resistance:
#Wan2019_SupplData.xlsx lists 667 genes the authors identifed from previous work as having a potential role in insecticide resistance.
#The authors then investigated differential SNPs in these genes using GWAS between resistant and susceptible insects

#Wan 667 list    
#Can we find the 667 genes in our assembly as well?
#Take the chromosomal region listed in Wan2019_SupplData.xlsx, extract the sequence from the Codling Moth genome as fasta sequence and BLAST against FCM assembly.
cd $PBS_O_WORKDIR

module load app/bedtools
module load app/samtools

#index fasta file
echo "starting fasta index\n"
samtools faidx $DIR/CodlingMothGeome.fasta
echo "completed fasta index\n"

#run bedtools to extract fasta sequences from Codling moth genome based on chromosomal regions given in Wan2019_SupplData.xlsx
echo "starting bedtools\n"
bedtools getfasta -fi $DIR/CodlingMothGeome.fasta -bed $DIR/Wan2019_SupplData.csv -fo $DIR/Wan2019_SupplData.fasta -name
echo "finsihed bedtools\n"

#BLAST Wan2019_SupplData.fasta sequences against fcmgenomev2 assembly to obtain matches in our assembly
echo "starting blastn\n"
blastn -db $DIR/fcmv2 -outfmt "7 qacc sacc evalue pident mismatch gapopen qstart qend sstart send length" -query $DIR/Wan2019_SupplData.fasta \
-out $DIR/Wan2019_SupplData_BLASTn_fcmgenomev2 \
-num_threads $NP
echo "done blastn\n"

#Find lengths of fasta sequences in Wan2019_SupplData.fasta by subtracting start and end sequence position in Wan2019_SupplData.xlsx. Add lengths to that file.
#Filter each BLAST hit in Wan2019_SupplData_BLASTn_fcmgenomev2 to the best hits that are (a) longest; (b) greater than 70% similarity and place in ‘filtered’ folder.
grep -v '#' Wan2019_SupplData_BLASTn_fcmgenomev2 | awk '$4 > "70"' | sort -n -k 11 | LC_NUMERIC=C awk '$11 > 100' | head

##########################Now we have a dataset, of insecticide resistance genes in our assembly->Wan2019_SupplData_BLASTn_fcmgenomev2###############################

#Wan SNPs           
#Wan SuplTable 25 lists Insecticide Resistance genes showing significant differences between resistant and susceptible moths. Gene IDs are given
#Find gene IDs in SuplTable 25 in the filtered BLAST hits generated above -> CpomSeqIDList
#Obtain sequences for gene IDs from the Codling Moth CDS sequences using gistfile.py
python3 gistfile.py CodlingMoth_CDS.fasta CpomSeqIDList > CpomIRSNPSeqs.fa

#Find the matching FCM scaffold number for these Seq IDs. -> FCMSeqIDList
#Obtain matching FCM sequences for gene IDs from the fcm genome assembly using gistfile.py
python3 gistfile.py fcmgenome_v2.fasta FCMSeqIDlist > FCMIRSNPSeqs.fa

#Align the CPOM gene and FCM sequences and see if the Nucleotide variation in the Wan Supplementary table is present
#Align sequences using MAFFT to obtain visual representation of differences
mafft --reorder --adjustdirection --auto --clustalout mafftinput.txt > out.txt #where mafftinput.txt consists of sequences to align from the files CpomIRSNPSeqs.fa and FCMIRSNPSeqs.fa

#***Note, the fasta files output that I want to use in MAFFT, have line breaks which makes Mauve unable to read them properly!
#Use this script to remove line breaks prior to MAFFT:
cat input.fasta | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > joinedlineoutput.fasta

#################Now we have a visual representation of Codling Moth sequences associated with insecticide resistance aligned to FCM sequences showing sequence variants between the two######################






























