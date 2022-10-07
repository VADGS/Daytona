# Daytona_combine
A combine of Daytona and Daytona_cl

## What to do
The pipeline is a combine of Daytona and Daytona_cl. The SARS-CoV-2 sequences from clearlabs or other sequencing machine, such as illumina, can all be analyzed by the pipeline.

## Prerequisites
Nextflow should be installed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed.

Singularity is also needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.

In addition, the below docker container images are needed in the pipeline. These images should be downloaded to the directory /apps/staphb-toolkit/containers/ in your local computer. You can find them from ncbi/sra-human-scrubber (https://hub.docker.com/r/ncbi/sra-human-scrubber) and StaPH-B/docker-builds (https://github.com/StaPH-B/docker-builds).

1. fastqc_0.11.9.sif
2. trimmomatic_0.39.sif
3. bbtools_38.76.sif
4. multiqc_1.8.sif
5. bwa_0.7.17.sif
6. samtools_1.12.sif
7. vadr_1.3.sif
8. pangolin_4.1.2-pdata-1.13.sif
9. nextclade_2021-03-15.sif
10. sra-human-scrubber_1.1.2021-05-05.sif

## How to run
### If the sequence datasets are NOT from clearlabs: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz". Test data can be found in the directory /fastqs/testdata. If you want to use the test data, copy them to the directory /fastqs.
2. open file "parames.yaml", set the parameters. 
3. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_sc2_combine.sh"

### If the sequence datasets are from clearlabs: 
1. put your data files into directories /fastqs, /bams, and /assemblies. Test data can be found in these directories. If you want to use the test data, copy them to the directories /fastqs, /bams, and /assemblies. 
2. open file "parames_clearlabs.yaml", set the parameters. 
3. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_sc2_combine.sh clearlabs"

## Results
All results can be found in the directory /output.

