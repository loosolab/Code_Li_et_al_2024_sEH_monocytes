# Reference Flatfiles and Samples

## RNA-Seq

### Prepare Reference GTF/FASTA

 - organism/assembly:	Ensembl mus musculus release 101, mm10 (https://ftp.ensembl.org/pub/)
 - ref_gtf:				mus_musculus.101.mainChr.gtf
 - star_index: 			Prepare STAR genome index as described in STAR manual

### Prepare Sample FASTQ
 - Copy sample FASTQ files to subfolder ./raw/

### Samples: RNA-Seq 1

```
prefix	fastq
M2c-WT_1	Wittig_Lib_XM_M2c-WT_1_R1.fastq.gz
M2c-WT_2	Wittig_Lib_XM_M2c-WT_2_R1.fastq.gz
M2c-WT_3	Wittig_Lib_XM_M2c-WT_3_R1.fastq.gz
M2c-KO_1	Wittig_Lib_XM_M2c-KO_1_R1.fastq.gz
M2c-KO_2	Wittig_Lib_XM_M2c-KO_2_R1.fastq.gz
M2c-KO_3	Wittig_Lib_XM_M2c-KO_3_R1.fastq.gz
```
