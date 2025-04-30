#!/bin/bash
#PBS -l walltime=400:00:00
#PBS -l select=1:ncpus=16:mpiprocs=36:mem=100GB
eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/env_oligominer

cd /data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test
#Negative genome generation 
#input files : 
ID_list="/data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/list_of_chrIDs.txt"
SEQ_FILE="/data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/K_Hell_Genome.fna"
JELLYFISH_DB="Sordaria_km.jf"

#Modify the jf nomenclature
jellyfish count -m 17 -s 40M -o Sordaria_km.jf -t 8 -L 2 $SEQ_FILE
cat $SEQ_FILE | grep '>'  | cut -f 1 -d ' ' |  sed 's/>//g'   > list_of_chrIDs.txt
faidx --split-files ${SEQ_FILE}

while read -r chrID
do
    echo $chrID
    mkdir -p $chrID
    mv "${chrID}.fna" "$chrID"
    ID_FILE=$(mktemp)
    echo "$chrID" > "$ID_FILE"
    #Output file name for the negative sequences
    OUTPUT_FILE="${chrID}_negative.fna"

    #Run the awk command to filter sequences not matching the current ID
    awk -v  target_id=">$chrID" '{
        if ($0 ~/^>/) {
            f = ($0 !~ target_id); 
        }
        if (f) {
            print $0;
        }
    }' "$SEQ_FILE" > "$chrID/$OUTPUT_FILE"
    echo "Generated negative file: $OUTPUT_FILE"         
    
    bowtie2-build $chrID/$OUTPUT_FILE $chrID/${chrID}_negative
    echo "Generated negative genome: ${chrID}"
    
    for L in {38..50}
    do
        (
            PROBE_DIR="probe_${L}/${chrID}"
            mkdir -p $PROBE_DIR

            python blockParse.py -l $L -L $L -f "$chrID/${chrID}.fna" -o "$PROBE_DIR/${chrID}_${L}"
            echo "Genetared oligo file: ${chrID}_${L} for temperature $T"

            bowtie2 -x $chrID/${chrID}_negative -U $PROBE_DIR/${chrID}_${L}.fastq -p 8 --no-hd -t -k 100 --very-sensitive-local -S $PROBE_DIR/${chrID}_${L}.sam
            echo "Generated alignment for : ${chrID}_${L} with basic setting (42:47)"

            python outputClean.py -0 -f $PROBE_DIR/${chrID}_${L}.sam -o $PROBE_DIR/${chrID}_${L}_cleaned
            echo "Cleaned oligo for : ${chrID}_${L} with basic setting (42:47)"

            python kmerFilter.py -f $PROBE_DIR/${chrID}_${L}_cleaned.bed -m 17 -j $JELLYFISH_DB -k 4 -o $PROBE_DIR/${chrID}_${L}_kmer_filtered
            echo "kmer filtered oligo for : ${chrID}_${L} with basic setting (42:47)"

            echo "All done for $chrID at temperature $T"
        ) &
    done
    wait
done < "$ID_list"
