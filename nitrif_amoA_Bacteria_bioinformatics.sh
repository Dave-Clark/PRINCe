# Bacterial amoA
cd /media/dave/storage/PRINCE/data/sequenceData/run_3/190612_M04364_0055_000000000-CJ9R5/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer

# target samples shared with AOA first
for f in AOA*-*R1*.fastq.gz; do
cutadapt -g "GGGGTTTCTACTGGTGGT;min_overlap=17" -G "CCCCTCKGSAAAGCCTTCTTC;min_overlap=20" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOB/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOB/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

# target samples that are AOB/comammox only
SAMPS=$(ls *R1* | grep -v "AOA*" | grep -v "Undetermined*")
for f in ${SAMPS}; do
cutadapt -g "GGGGTTTCTACTGGTGGT;min_overlap=17" -G "CCCCTCKGSAAAGCCTTCTTC;min_overlap=20" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOB/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOB/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOB/

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# try pair end alignment. seems to work for ~60% of seqs
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 150 -y -A -5 -c --overlap_len_require 10 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
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
grep "^[^>]" C8x4x5CMxS_S188.fna | awk "{print length}" > seqLengths.txt

# take all reads between > 448
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 448 { print $0 "\n" seq }' $f > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in AOA*_len_filtered.fna; do
	samp=">$(echo $f | awk -F "[-._]" '{print $2}')-"
	sed -i "s/>/$samp/" $f
done

SAMPS=$(ls *_len_filtered.fna | grep -v "AOA*")
for f in ${SAMPS}; do
	samp=">$(echo $f | awk -F "[-._]" '{print $1}')-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# discard samples with fewer than 200 seqs
rm -f {AOAxC1x1x5CMxWPL-C3x2x5CMxS_S79_len_filtered.fna,AOAxC1x4x5CMxW-C10x5x5CMxW_S95_len_filtered.fna,AOAxC3x1x5CMxS-C6x2x5CMxS_S91_len_filtered.fna,AOAxC8x2x5CMxWPL-C8x5x5CMxW_S75_len_filtered.fna,AOAxS1x4x5CMxS-NEG3_S19_len_filtered.fna,AOAxS3x2x5CMxS-C10x2x5CMxS_S29_len_filtered.fna,AOAxS3x6x5CMxS-C1x4x5CMxS_S56_len_filtered.fna,C2x2x5CMxWPL_S118_len_filtered.fna,NEG1_S166_len_filtered.fna,NEG2_S136_len_filtered.fna,S1x6x5CMxW_S116_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 1022541

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/AOB_amoA -N -i 0.5 -l 50 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOB_amoA/derep_AOB_amoA.fna uniq_DNA_seqs.fna


cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' AOB_amoA_corr_prot.fasta > uc_AOB_amoA_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_AOB_amoA_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_AOB_amoA_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout AOB_AAV.txt -qmask none -dbmask none
