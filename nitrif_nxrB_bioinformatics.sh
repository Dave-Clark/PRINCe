# Nitrospira nxrB

cd /media/dave/storage/PRINCE/data/sequenceData/run_4/190819_M04364_0067_000000000-CF725/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in nosZ*R1*.fastq.gz; do
cutadapt -g "TACATGTGGTGGAACA;min_overlap=15" -G "CGGTTCTGGTCRATCA;min_overlap=15" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/nitrif_nxrB/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/nitrif_nxrB/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/nitrif_nxrB/

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# quality filter and merge sequences
# remove seqs if > 20% is under q20
# below 200 bases after trimming
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 200 -y -A -5 -c --overlap_len_require 15 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
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

# sort sequences by length to check for artefacts
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' nosZxS4x5x5CMxW_S206.fna |\
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
    sort -k1,1n | cut -f 2- | tr "\t" "\n" > lenSortedSeqs.fa

# take all reads between 450 - 455 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 450 { print $0 "\n" seq }' $f | awk '!/^>/ { next } { getline seq } length(seq) <= 455 { print $0 "\n" seq }' > ${f%.fna}_len_filtered.fna
done

# strip gene name label
for f in *.fna; do
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

# discard samples with fewer than 1511 seqs (n = 2)
rm -f {C6x1x5CMxWPL_S169_len_filtered.fna,S3x3x5CMxS-nirKxC8x1x5CMxWPL_S145_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 2198680

# dereplicate downloaded fungene database seqs
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/nitrospira_nxrB/fungene_9.7_nxrB_nitrospira_290_unaligned_protein_seqs.fa --fastaout /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/nitrospira_nxrB/derep_nxrB.fna --minuniquesize 1

# dereplicate reads by sample to preserve abundance info
# don't discard singletons
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

mkdir framebotOut

# translate to frameshift corrected amino acid sequences with per sample
# abundance information preserved.
java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/nitrospira_nxrB -N -i 0.5 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/nitrospira_nxrB/derep_nxrB.fna uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' nitrospira_nxrB_corr_prot.fasta > uc_nitrospira_nxrB_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_nitrospira_nxrB_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_nitrospira_nxrB_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout nxrB_AAV.txt -qmask none -dbmask none

# remove singleton AAVs in R
