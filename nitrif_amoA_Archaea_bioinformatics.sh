# Archaeal amoA
cd /media/dave/storage/PRINCE/data/sequenceData/run_3/190612_M04364_0055_000000000-CJ9R5/Data/Intensities/BaseCalls

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in AOA*R1*.fastq.gz; do
cutadapt -g "ATGGTCTGGCTWAGACG;min_overlap=16" -G "GCCATCCATCTGTATGTCCA;min_overlap=19" -j 0 --discard-untrimmed -o /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOA/${f%001.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOA/${f%1_001.fastq.gz}2_primer_trimmed.fastq $f ${f%R1_001.fastq.gz}R2_001.fastq.gz
done

cd /media/dave/storage/PRINCE/data/sequenceData/nitrif_AOA/

# make directories for merged sequence files, html output and discarded sequences
mkdir {qualTrimmedSeqs,htmlReports}

# quality filter
# remove seqs if > 20% is under q20
# below 200 bases after trimming
# paired reads won't merge - too long an amplicon
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -q 20 -u 20 -l 200 -y -A -5 -o qualTrimmedSeqs/${f%L001_R1_primer_trimmed.fastq}qualTrimmed.fastq
  -h htmlReports/${f%L001_R1_primer_trimmed.fastq}quality_report.html
done

cd qualTrimmedSeqs

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
  seqtk seq -a $f > ${f%_qualTrimmed.00.0_0.cor.fastq}.fna
done

### NEED TO LOOK AT LENGTH FILTERS
# get read lengths for one file
grep "^[^>]" AOAxC5x4x5CMxW-S4x6x5CMxW_S7.fna | awk "{print length}" > seqLengths.txt

# take all reads between > 280 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 280 {print $0 "\n" seq}' $f > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "[-._]" '{print $1}')-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
grep -c "^>" *_len_filtered.fna > finalLibSizes.txt

# discard samples with fewer than 3000 seqs
rm -f {AOAxNEG1_S165_len_filtered.fna,AOAxNEG2_S135_len_filtered.fna,AOAxC6x2x5CMxS-C3x1x5CMxS_S92_len_filtered.fna,AOAxNEG3-S1x4x5CMxS_S20_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 5051717

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

# subsample AOA database otherwise framebot will crash!
seqtk sample -s100 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOA_amoA/derep_AOA_amoA.fna 2000 > /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOA_amoA/subSampled_AOA.fna

mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/AOA_amoA -N -i 0.5 -l 50 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOA_amoA/subSampled_AOA.fna uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' AOA_amoA_corr_prot.fasta > uc_AOA_amoA_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_AOA_amoA_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_AOA_amoA_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout AOA_AAV.txt -qmask none -dbmask none
