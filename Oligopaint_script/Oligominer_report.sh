#!/bin/bash
#PBS -l walltime=400:00:00
#PBS -l select=1:ncpus=16:mpiprocs=36:mem=100GB

eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/env_oligominer
# Extract unique chromosome IDs from the FASTQ file
# Initialize report file
ID_list="list_of_chrIDs.txt"
# Initialize report file
cd /data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test
for L in {38..50}
do
REPORT_FILE="pipeline_report_${L}.tsv"
echo -e "Chromosome\tChr Length\tProbe Length\tMelting Temp min\tMelting temp max\tAverage Melting Temp\tblockParse Candidates\tOutputclean\tKmerfilter (SAM)\tRatio probe/kb" > "$REPORT_FILE"
while read -r chrID; do
    Length=$(bioawk -c fastx '{ print length($seq) }' < $chrID/${chrID}.fna)
    BlockParse="probe_${L}/$chrID/${chrID}_${L}.fastq"
    Outputclean="probe_${L}/$chrID/${chrID}_${L}_cleaned.bed"
    Kmerfilter="probe_${L}/$chrID/${chrID}_${L}_kmer_filtered.bed"
# Loop through each chromosome and count occurrences in each file

    # Count occurrences in FASTQ
    BlockParse_COUNT=$(awk -v chr="$chrID" '/^@/ && $0 ~ chr {count++} END {print count+0}' "$BlockParse")
    echo $BlockParse_COUNT

    # Count occurrences in BED
    Outputclean_COUNT=$(grep -c "^$chrID" "${Outputclean}")
    echo $Outputclean_COUNT

    # Count occurrences in SAM
    Kmerfilter_COUNT=$(grep -c "^$chrID" "${Kmerfilter}")
    echo $Kmerfilter_COUNT
    #extract Tm
    avg_tm_60=$(awk -v chr="$chrID" '$1 == chr { tm=$NF; sum += tm; count += 1; } END {if (count > 0) {avg_tm = sum / count; printf "%.2f\n", avg_tm; } else { print "NA"; }}' "$Kmerfilter")
    max_tm_60=$(awk -v chr="$chrID" '$1 == chr { tm=$NF; if (tm > max_tm) { max_tm = tm; } } END { print int(max_tm + 0.5); }' "$Kmerfilter")
    min_tm_60=$(awk -v chr="$chrID" 'BEGIN {min_tm = 1000} $1 == chr { tm=$NF; if (tm < min_tm) { min_tm = tm; } } END { print int(min_tm + 0.5); }' "$Kmerfilter")

    #Determine Ratio 
    Ratio=$(awk "BEGIN {print $Kmerfilter_COUNT/ $Length*1000}" )
    echo $Ratio

    # Append results to the report
    echo -e "${chrID}\t${Length}\t60\t${min_tm_60}\t${max_tm_60}\t${avg_tm_60}\t${BlockParse_COUNT}\t${Outputclean_COUNT}\t${Kmerfilter_COUNT}\t${Ratio}" >> "$REPORT_FILE"
    
done < "$ID_list"
echo "Report generated: $REPORT_FILE"
done

