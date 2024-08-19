##
## Notes
##
## Prepare reference indices and flatfiles according to ../CONFIG.md
##
## variables:
## ${sample_prefix}:            label prefix for fastq
## ${ref_gtf}:                  reference genome gtf by Ensembl; e.g. mus_musculus.101.mainChr.gtf
## ${star_index}:               star reference index


##
## Repeat the following steps per fastq
##

## Trimming: Trimmomatic
java -jar trimmomatic.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping: STAR
STAR --genomeDir ${star_index} --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterMatchNmin 20 --alignIntronMax 200000 --alignMatesGapMax 2000 --alignEndsProtrude 10 ConcordantPair --outMultimapperOrder Random --limitOutSAMoneReadBytes 10000000 --sjdbOverhang 100 --outFilterMultimapNmax 1

## Generate BigWig normalized to mapped reads per sample: python package deeptools bamCoverage
bamCoverage -b ./star/${sample_prefix}.bam -o ./star/${sample_prefix}.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig

## Count reads per genes: subread featureCounts
featureCounts -T 16 -o ./star/${sample_prefix}_fcounts.txt -t exon -g gene_id -s 1 -a ${ref_gtf} ./star/${sample_prefix}.bam


##
## Combine per-sample results
##

## Join all sample counts to raw count matrix: =counts.matrix
merge_counts.pl ./star/*_fcounts.txt >counts.matrix

## RNA differential analysis: R with DESeq2 package
## de_analysis.R script has to be adapted to include samples/conditions of the respective dataset
R --vanilla -q < de_analysis.R

## Normalize matrix for downstream plots: R with DESeq2 package -> counts.matrix.norm
deseq_norm.R -m counts.matrix
