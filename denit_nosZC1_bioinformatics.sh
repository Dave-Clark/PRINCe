# nosZ Clade I
cd /media/dave/storage/PRINCE/data/sequenceData/run_4/190819_M04364_0067_000000000-CF725/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in nosZ*R1*.fastq.gz; do
cutadapt -g "WCSYTGTTCMTCGACAGCCAG;min_overlap=20" -G "ATGTCGATCARCTGVKCRTTYTC;min_overlap=22" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/denit_nosZC1/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/denit_nosZC1/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/denit_nosZC1/

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# quality filter and merge sequences
# remove seqs if > 20% is under q20
# below 200 bases after trimming
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 130 -y -A -5 -c --overlap_len_require 30 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
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
grep "^[^>]" nosZxS4x5x5CMxW_S206.fna | awk "{print length}" > seqLengths.txt

# take all reads between 200 - 210 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print $0 "\n" seq }' $f | awk '!/^>/ { next } { getline seq } length(seq) <= 210 { print $0 "\n" seq }' > ${f%.fna}_len_filtered.fna
done

# strip gene name label
for f in *_len_filtered.fna; do
    [ -f "$f" ] || continue
    mv "$f" "${f//nosZx/}"
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "_" '{print $1}')-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# discard samples with fewer than 1000 seqs
rm -f {NEG1_S176_len_filtered.fna,NEG2_S161_len_filtered.fna,NEG3-nirKxS1x4x5CMxS_S103_len_filtered.fna,C6x1x5CMxWPL_S169_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 5840370

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/nosZ_C1 -N -i 0.5 -l 60 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/nosZ_CladeI/derep_nosZ_C1.fna uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' nosZ_C1_corr_prot.fasta > uc_nosZ_C1_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_nosZ_C1_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_nosZ_C1_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout nosZC1_AAV.txt -qmask none -dbmask none
