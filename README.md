GeneQC Python Tutorial

Dependencies and preliminary data
The GeneQC package requires Python 3 to execute, includes the blast+ and SAMtools libraries, GeneQC takes Reference Genome, Annotation file, Read alignment (mapped by HISAT2, compressed to BAM file by SAMtools) as input, you may take the following steps to obtain the preliminary files. The GeneQC will generate feature extraction and modeling results(D-score).
The GeneQC also requires Scikit-learn package. Please following this link to install scikit-learn. (http://scikit-learn.org/stable/install.html). Instead installing any package for Python, user can use anaconda to execute the GeneQC (Recommend). Download link: (https://www.anaconda.com/download/)

Instructions

For running GeneQC for plant:
Three inputs: reference genome, annotation file, sam or bam file should in the same folder of GeneQC2.py. The outputs will be generated in this folder as well.

GeneQC2.py [1] [reference genome] [standard gff annotation file] [sam or bam file]

Example: A.thaliana
GeneQC2.py 1 Athaliana_167_TAIR9.fa Athaliana_167_TAIR10.gene.gff3 ERR1297323.bam

For running GeneQC for aniaml:

Step1: Create new defined transcripts and new defined transcripts annotation:
extract_transcript_seq_gff.py [reference genome] [standard gff annotation file] [new defined transcripts sequence file] [new defined transcripts gff annotation file]

Example:
extract_transcript_seq_gff.py GCF_000001405.37_GRCh38.p11_genomic.fna GCF_000001405.37_GRCh38.p11_genomic.gff human_transcripts_seq.fa human_transcripts_seq.gff

Step2: Do RNA-seq mapping work with the new mapping results, use following commands:

hisat2-build -f [new defined transcripts sequence file] ./hisatindex/Humo

For Bulk RNA-seq data:
hisat2 -x ./hisatindex/Humo -k 10 -1 [fastq file 1] -2 [fastq file 2] -S [sam file]

For Single cell RNA-seq data:
hisat2 -x ./hisatindex/Humo -k 10 [fastq file] -S [sam file]

Example:
hisat2-build -f human_transcripts_seq.fa ./hisatindex/Humo
hisat2 -x ./hisatindex/Humo -k 10 SRR491087.fastq -S SRR491087.sam

Step3: Run GeneQC animal:
GeneQC2.py 2 [new defined transcripts sequence file] [new defined transcripts gff annotation file] [sam file]

Example:
GeneQC2.py 2 human_transcripts_seq.fa human_transcripts_seq.gff SRR491087.sam

Outputs:
The GeneQC will generate two results in the same folder of code:
1.	Feature extraction results (bamfilename_out.txt)
2.	Modeling results (bamfilename_out.csv)
