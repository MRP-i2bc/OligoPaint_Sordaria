#!/bin/bash
#SBATCH --time 1000:00:00
#SBATCH --cpus-per-task 32
#SBATCH --mem 100G

eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/genome_analysis

cd /data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport


python3 CentroReport_.py  /data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta \
    -o /data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport/Results/Centromere_density_report_WT_J2_10000_10_trf_RM_test_10000.html \
    -k 10 \
    -w 10000 \
    -t 32 \
    --run-trf

python3 CentroReport_.py /data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta \
    -o /data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport/Results/Centromere_density_report_WT_J2_10000_10_trf_RM_test_1000.html \
    -k 10 \
    -w 1000 \
    -t 32 \
    --run-trf

python3 CentroReport_.py  /data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta \
    -o /data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport/Results/Centromere_density_report_WT_J2_10000_40_trf_RM_test_10000.html \
    -k 40 \
    -w 10000 \
    -t 32 \
    --run-trf

python3 CentroReport_.py  /data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta \
    -o /data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport/Results/Centromere_density_report_WT_J2_10000_40_trf_RM_test_1000.html \
    -k 40 \
    -w 1000 \
    -t 32 \
    --run-trf