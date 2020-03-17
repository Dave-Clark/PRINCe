# anammox hzo

# starting directory
cd /media/dave/storage/PRINCE/data/sequenceData/run_2/190529_M04364_0052_000000000-C36NK/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didnt contain both the forward and reverse primer
for f in *R1*.fastq.gz; do
cutadapt -g "AAGACNTGYCAYTGGGGWAAA;min_overlap=20" -G "GACATACCCATACTKGTRTANACNGT;min_overlap=24" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/anammox_hzo/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/anammox_hzo/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/anammox_hzo
rm -f Undetermined*

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# quality filter and merge sequences
# remove seqs if > 20% is under q20
# below 200 bases after trimming
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 150 -y -A -5 -c --overlap_len_require 50 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
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

# take all reads between 175 - 178 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 175 { print $0 "\n" seq }' $f | awk '!/^>/ { next } { getline seq } length(seq) <= 178 { print $0 "\n" seq }' > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "_" '{print $1}')-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# discard samples with fewer than 1000 seqs
rm -f {C5x2x5CMxS_S34_len_filtered.fna,C6x2x5CMxS_S43_len_filtered.fna,C9x2x5CMxS_S35_len_filtered.fna,S1x5x5CMxS_S97_len_filtered.fna,S3x3x5CMxS_S49_len_filtered.fna,S3x6x5CMxS_S36_len_filtered.fna,S4x1x5CMxS_S18_len_filtered.fna,S4x1x5CMxWPL_S114_len_filtered.fna,S4x2x5CMxS_S2_len_filtered.fna,S4x2x5CMxWPL_S68_len_filtered.fna,S4x3x5CMxS_S67_len_filtered.fna,S4x3x5CMxWPL_S99_len_filtered.fna,S4x4x5CMxS_S94_len_filtered.fna,S4x5x5CMxS_S111_len_filtered.fna,S4x5x5CMxW_S110_len_filtered.fna,S4x6x5CMxS_S52_len_filtered.fna,S4x6x5CMxW_S12_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 1897257

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

mkdir framebotOut
/media/dave/storage/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut
/hzo -N -i 0.5 -l 40 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/ana
mmox_hzo/derep_hzo_seqs.fna uniq_DNA_seqs.fna


# up to here
# need to run through framebot once site is working
# need hzo refset or use online tool
