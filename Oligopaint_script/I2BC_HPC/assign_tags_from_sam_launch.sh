#!/bin/bash
#SBATCH --time 1000:00:00
#SBATCH --cpus-per-task 32
#SBATCH --mem 100G

eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/genome_analysis

cd /data/work/I2BC/antoine.bochet

ID_list="/data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/list_of_chrIDs.txt"

while read -r chrID
do  
    #python3 filter_fastq_by_list.py probe_50/${chrID}/${chrID}_50.fastq /data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/purified_probes_list_WT_J2.txt probe_50/${chrID}/${chrID}_50_purif.fastq
    echo $chrID
    #samtools view -bS probe_50/${chrID}/${chrID}_50_kmer_filtered_aligned_WT_J2.sam.unique.sam >  probe_50/${chrID}/${chrID}_aligned_WT_J2.unique.bam
    #bedtools bamtofastq -i  probe_50/${chrID}/${chrID}_aligned_WT_J2.unique.bam -fq  probe_50/${chrID}/${chrID}_aligned_WT_J2.unique.fastq
    python3 assign_tags_from_sam_coords.py \
        probe_50/${chrID}/${chrID}_50_kmer_filtered_aligned_WT_J2.sam.unique.sam \
        /data/work/I2BC/antoine.bochet/Primers_tag.tsv \
        probe_50/${chrID}/${chrID}_aligned_WT_J2.unique.fastq \
        ./Output_Tagged_From_SAM_denisty_FINAL_Chr
done < "$ID_list"

