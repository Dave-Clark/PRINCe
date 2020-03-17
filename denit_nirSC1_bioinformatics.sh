# nirS

# starting directory
cd /media/dave/storage/PRINCE/data/sequenceData/run_2/190529_M04364_0052_000000000-C36NK/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in *R1*.fastq.gz; do
cutadapt -g "ATCGTCAACGTCAARGARACVGG;min_overlap=22" -G "TTCGGGTGCGTCTTSABGAASAG;min_overlap=22" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/denit_nirSC1/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/denit_nirSC1/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/denit_nirSC1
rm -f Undetermined*

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# quality filter and merge sequences
# remove seqs if > 20% is under q20
# below 200 bases after trimming
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 150 -y -A -5 -c --overlap_len_require 15 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
done

cd mergedSeqs

### Error correct all seqs ###
for f in *.fastq;
	do ~/bioinformatics/SPAdes-3.13.0-Linux/bin/spades.py -o ${f%q_trimmed_merged.fastq}_errorCorrected --only-error-correction -s $f -t 8 -m 32 --disable-gzip-output;
done

# copy all error corrected sequences and convert to fastq and length filter
mkdir allErrorCorrected
mv *errorCorrected/corrected/*.fastq allErrorCorrected

cd allErrorCorrected

# convert to fasta format
for f in *.fastq; do
  seqtk seq -a $f > ${f%_q_trimmed_merged.00.0_0.cor.fastq}.fna
done

# get read lengths for one file
grep "^[^>]" C8x6x5CMxW_S64.fna | awk "{print length}" > seqLengths.txt

# take all reads between 360 - 380 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 360 { print $0 "\n" seq }' $f | awk '!/^>/ { next } { getline seq } length(seq) <= 380 { print $0 "\n" seq }' > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "_" '{print $1}')-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# discard samples with fewer than 1000 seqs
rm -f {S1x1x5CMxS_S32_len_filtered.fna,S4x6x5CMxS_S52_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 1340619

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/nirS -N -i 0.5 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/nirS/fungene_nirS.fa uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' nirS_corr_prot.fasta > uc_nirS_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_nirS_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_nirS_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout nirSC1_AAV.txt -qmask none -dbmask none
