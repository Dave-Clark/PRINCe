mkdir primerTrimmed

# cut primer sequences
# remove any sequences that didnt contain both the forward and reverse primer
for f in *R1*.fastq.gz; do
cutadapt -g "CCTACGGGNGGCWGCAG;min_overlap=14" -G "GACTACHVGGGTATCTAATCC;min_overlap=18" -j 0 --discard-untrimmed -o primerTrimmed/${f%001.fastq.gz}primer_trimmed.fastq -p primerTrimmed/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# quality filter and merge sequences
# remove seqs if > 20% is under q20
# below 200 bases after trimming
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 200 -y -A -5 -c --overlap_len_require 15 --merge --merged_out mergedSeqs/${f%L001_R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
done

# check sample sizes after quality trimming and pair end alignment
grep -c "^@M0" *.fastq > library_sizes_postQC
# look ok!

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

# take all reads between 350 - 450 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 350 { print $0 "\n" seq }' $f | awk '!/^>/ { next } { getline seq } length(seq) <= 450 { print $0 "\n" seq }' > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "_" '{print $1}')_"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# remove two small samples
# S3x6x5CMxWPL_S131_len_filtered.fna
# C3x2x5CMxWPL_S112_len_filtered.fna
# keep negatives in
rm -f {S3x6x5CMxWPL_S131_len_filtered.fna,C3x2x5CMxWPL_S112_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 3786559

# dereplicate sequences for OTU clustering, discard singletons
vsearch --derep_fulllength allSequences.fna --output derepSeqs.fna --minuniquesize 2 --relabel OTU --fasta_width 0 --sizeout

# ref based chimera check
vsearch --uchime3_denovo derepSeqs.fna --chimeras bacChimeras.fna --nonchimeras chimFreeSeqs.fna

# sort by cluster size
vsearch --sortbysize chimFreeSeqs.fna --output derepSorted.fna

# pick otu centroids
vsearch --cluster_size derepSorted.fna --id 0.97 --centroids bacCentroids.fna

# create map reads to OTUs
vsearch --usearch_global allSequences.fna --db bacCentroids.fna --id 0.97 --log bacteriaOtus.log --otutabout bacOtuTab.txt

# assign taxonomy to OTU centroid sequences
java -Xmx4g -jar ~/bioinformatics/RDPTools/classifier.jar -f filterbyconf -o bacteriaTaxonomy.txt bacCentroids.fna
