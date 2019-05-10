# GeneQC Tutorial

The GeneQC package requires Python 3 to execute, includes the blast+ and SAMtools libraries, GeneQC takes Reference Genome, Annotation file, Read alignment (mapped by HISAT2, compressed to BAM file by SAMtools) as input, you may take the following steps to obtain the preliminary files. The GeneQC will generate feature extraction and modeling results(D-score).

The package contains:

1. Minimum package of Blast for makeblastdb and Blastn

2. Samtools package version 1.2.1

3. GeneQC python code

4. A special script "extract_transcript_seq_gff.py" which used to transform the genome sequences to transcript sequences and the corresponding annotation file.

## Hardware Requirement:

Hardware requirements are based on the size of Reference Genome and Annotation file. High performance cluster preferred.

## Software Requirement:

GeneQC package requires Python 3 to execute. We recommended user use anaconda 5 to execute the GeneQC. Most of high performance clusters have already installed the anaconda, so user can load anaconda directly using following code.

1. Check available modules in the cluster:
```{r,engine='bash',eval=FALSE}
module avail
```
2. find anaconda (4 or above) in available modules and load it. (for example, if anaconda5/5.0.0-3.6 is in the available modules list):
```{r,engine='bash',eval=FALSE, module}
module load anaconda5/5.0.0-3.6
```

## Installation
The source code of GeneQC is freely available at: https://github.com/maqin2001/GeneQC

To install GeneQC, first download the zip file manually from github, or use the code below in Unix:
```{r,engine='bash',eval=FALSE, download}
cd your_folder_path
wget https://github.com/maqin2001/geneqc/archive/master.zip
```
Unzip the file:
```{r,engine='bash',eval=FALSE, unzip}
unzip master.zip
```

## Input data preperation

Seven sample species of data (Reference Genome, Annotation file, Read alignment (mapped by HISAT2, compressed to BAM file by SAMtools) can be downloaded from our website (http://bmbl.sdstate.edu/GeneQC/result.html).

## Instructions (Plant Genome)

For running GeneQC for plant data: Three inputs (data) are required: reference genome, annotation file, sam or bam file should be uploaded to folder of "GeneQC_Python" under the GeneQC-master folder in the cluster.

Move to the path of folder of "GeneQC_Python"
```{r,engine='bash',eval=FALSE}
cd your_folder_path_of_"GeneQC_Python"
module load anaconda(chose_your_version)
```

Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py [1] [reference genome] [standard gff annotation file] [sam or bam file]
```

Example: A.thaliana
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 1 Athaliana_167_TAIR9.fa Athaliana_167_TAIR10.gene.gff3 ERR1297323.bam
```

The outputs will be generated in this folder as well. ERR1297323_out.txt will be feature extraction results. ERR1297323_out.csv will be D-scoure results.

## Instructions (Amimal Transcript)

For running GeneQC for animal data: Three inputs (data) are required: reference genome, annotation file, fastq file should be uploaded to folder of "GeneQC_Python" under the GeneQC-master folder in the clust

Move to the path of folder of "GeneQC_Python"
```{r,engine='bash',eval=FALSE}
cd your_folder_path_of_"GeneQC_Python"
module load anaconda(chose_your_version)
```

Step1: Create new defined transcripts and new defined transcripts annotation:
```{r,engine='bash',eval=FALSE}
python extract_transcript_seq_gff.py [reference genome] [standard gff annotation file] [new defined transcripts sequence file] [new defined transcripts gff annotation file]
```

Example: Humo sapiens
```{r,engine='bash',eval=FALSE}
python extract_transcript_seq_gff.py GCF_000001405.37_GRCh38.p11_genomic.fna GCF_000001405.37_GRCh38.p11_genomic.gff human_transcripts_seq.fa human_transcripts_seq.gff
```

### For Bulk RNA-seq data:
Step2: Do RNA-seq mapping work with the new mapping results, use following commands (the example used aligner HISAT2):
```{r,engine='bash',eval=FALSE}
module load hisat2
hisat2-build -f human_transcripts_seq.fa ./hisatindex/Humo
hisat2 -x ./hisatindex/Humo -k 10 -p 40 -1 SRR6029567_1.fastq -2 SRR6029567_2.fastq -S SRR6029567.sam
```

Step3: Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 2 human_transcripts_seq.fa human_transcripts_seq.gff SRR6029567.sam
```

### For Single cell RNA-seq data:
Step2: Do RNA-seq mapping work with the new mapping results, use following commands (the example used aligner HISAT2):
```{r,engine='bash',eval=FALSE}
module load hisat2
hisat2-build -f human_transcripts_seq.fa ./hisatindex/Humo
hisat2 -x ./hisatindex/Humo -k 10 SRR491087.fastq -S SRR491087.sam
```

Step3: Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 2 human_transcripts_seq.fa human_transcripts_seq.gff SRR491087.sam
```

The outputs will be generated in this folder as well. SRR6029567_out.txt will be feature extraction results. SRR6029567_out.csv will be D-scoure results.
