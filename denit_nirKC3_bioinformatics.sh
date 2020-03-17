# nirK Clade 3

cd /media/dave/storage/PRINCE/data/sequenceData/run_4/190819_M04364_0067_000000000-CF725/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in *nirK*R1*.fastq.gz; do
cutadapt -g "CATCGGCAACGGCATGYAYGGNGC;min_overlap=23" -G "CGACCATGGCCGTGGSWNACRAANGG;min_overlap=25" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/denit_nirKC3/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/denit_nirKC3/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/denit_nirKC3/
rm -f Undetermined*

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
grep "^[^>]" nirKxC8x3x5CMxS_S4.fna | awk "{print length}" > seqLengths.txt

# NOT LOOKING GOOD!
# try taking all and use framebot to screen out crap.
# take all reads between 360 - 380 bp
# add sample labels to fasta headers
for f in *.fna; do
	samp=">$(echo $f | cut -d "-" -f2 | cut -d "_" -f1)-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *.fna > finalLibSizes.txt

# discard samples with fewer than 1000 seqs
rm -f {S1x1x5CMxS_S32_len_filtered.fna,S4x6x5CMxS_S52_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 1340619
