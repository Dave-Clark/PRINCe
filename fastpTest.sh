mkdir {mergedSeqs,discardedSeqs,htmlReports}

cutadapt -g "CCTACGGGNGGCWGCAG;min_overlap=14" -G "GACTACHVGGGTATCTAATCC;min_overlap=18" -j 0 --discard-untrimmed -o out.1.fastq -p out.2.fastq C1x1x5CMxS_S138_L001_R1_001.fastq C1x1x5CMxS_S138_L001_R2_001.fastq

fastp -i out.1.fastq -I out.2.fastq -q 20 -u 20 -l 200 -y -A -5 -c --overlap_len_require 10 --merge --merged_out mergedSeqs/testmerged2.fastq -o discardSeqs/testdiscarded2.fastq -O discardSeqs/testR2_discarded2.fastq -h htmlReports/test_html_report2.html

seqtk seq -a mergedSeqs/testmerged2.fastq | grep "^[^>]" | awk "{print length}" > seqLengths2.txt

# all reads above 275nt
seqtk seq -a mergedSeqs/testmerged.fastq | awk '!/^>/ { next } { getline seq } length(seq) >= 350 { print $0 "\n" seq }'| awk '!/^>/ { next } { getline seq } length(seq) <= 450 { print $0 "\n" seq }' > lenFilteredSeqs.fna

for f in *R1_001.fastq; do
fastp -i $f -I ${f%R1_001.fastq}R2_001.fastq -q 20 -u 20 -l 250 -y -A -5 -c --overlap_len_require 10 --merge --merged_out mergedSeqs/${f%R1_001.fastq}merged.fastq -o discardSeqs/${f%001.fastq}discarded.fastq -O discardSeqs/${f%R1_001.fastq}R2_discarded.fastq -h htmlReports/${f%R1_001.fastq}_html_report.html
done

grep "^[^>]" allSeqs.fa | awk "{print length}" > seqLengths.txt

seqtk seq -a mergedSeqs/testmerged.fastq |
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  |\
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
    sort -k1,1n | cut -f 2- | tr "\t" "\n" > lenSortedSeqs.fa
