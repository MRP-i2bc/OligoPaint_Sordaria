#!/bin/bash
#SBATCH --time 1000:00:00
#SBATCH --cpus-per-task 32
#SBATCH --mem 100G
eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/env_oligominer


#Negative genome generation 
#input files : 
ID_list="/data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/list_of_chrIDs.txt"
SEQ_FILE="/data/work/I2BC/antoine.bochet/OligoPaint/Size_probe_test/K_Hell_Genome.fna"
JELLYFISH_DB="Sordaria_km.jf"
CONTIGS_FILE="/data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta"

# Output purified probes list file
PURIFIED_PROBES_LIST="purified_probes_list_WT_J2.txt"
#!/bin/bash

# --- Configuration ---

# Base directory containing chromosome-specific folders
BASE_PROBE_DIR="probe_50"

# Bowtie2 alignment options
BOWTIE2_OPTS="-N 0 --no-unal -f -a --sam-no-qname-trunc"
BOWTIE2_THREADS=32

# --- Sanity Checks ---
# ... (Vérifications pour bowtie2, CONTIGS_FILE, BASE_PROBE_DIR) ...
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools command not found. Please install samtools."
    exit 1
fi
if ! command -v grep &> /dev/null; then
    echo "Error: grep command not found."
    exit 1
fi
# *** AJOUT: Vérification pour bedtools ***
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools command not found. Please install bedtools."
    exit 1
fi

echo "Starting probe alignment and BED file generation..."

shopt -s nullglob

for chr_dir in "$BASE_PROBE_DIR"/*; do
    if [[ -d "$chr_dir" ]]; then
        chrID=$(basename "$chr_dir")
        echo "Processing chromosome directory: $chrID"

        processed_file=false

        for fasta_file in "$chr_dir"/*_filtered.fasta; do
             if [[ ! -f "$fasta_file" ]]; then
                 continue
             fi

             echo "  Input sequence file: $fasta_file"
             processed_file=true
             base_name=$(basename "$fasta_file" .fasta)

             # --- Modification des noms de fichiers intermédiaires ---
             ALIGNMENT_OUTPUT_SAM="${chr_dir}/${base_name}_aligned_WT_J2.sam"
             # Produire directement un BAM unique
             UNIQUE_BAM_OUTPUT="${chr_dir}/${base_name}_aligned_WT_J2.unique.bam"
             BED_OUTPUT="${chr_dir}/${base_name}_on_WT_J2_TEST.bed"

             # --- Perform alignment ---
             echo "    Running Bowtie2 alignment..."
             bowtie2 -x "$CONTIGS_FILE" $BOWTIE2_OPTS -p $BOWTIE2_THREADS -U "$fasta_file" -S "$ALIGNMENT_OUTPUT_SAM"

             if [[ $? -ne 0 ]]; then
                 echo "    Error: Bowtie2 alignment failed for $fasta_file. Skipping."
                 # rm -f "$ALIGNMENT_OUTPUT_SAM"
                 continue
             fi
             echo "    Alignment finished: $ALIGNMENT_OUTPUT_SAM"

             # --- Filter unique mappings (NH:i:1) AND convert to BAM directly ---
             echo "    Filtering unique mappings (NH:i:1) and converting to BAM..."

             # Combine header, filter alignment lines for NH:i:1, pipe to samtools view to create BAM
             # L'option -bS de samtools view prend SAM en entrée et sort BAM
             { samtools view -H "$ALIGNMENT_OUTPUT_SAM"; \
               samtools view "$ALIGNMENT_OUTPUT_SAM" | grep -w 'NH:i:1'; } | samtools view -bS - > "$UNIQUE_BAM_OUTPUT"

             # Vérifier si le BAM unique a été créé et contient des données
             # `samtools view -c` compte les alignements dans un BAM/SAM/CRAM
             unique_count=$(samtools view -c "$UNIQUE_BAM_OUTPUT" || echo 0) # Compte les alignements
             if [[ $? -ne 0 || $unique_count -eq 0 ]]; then
                 echo "    No uniquely mapping reads (NH:i:1) found or BAM creation failed. Skipping BED generation."
                 rm -f "$UNIQUE_BAM_OUTPUT" # Supprimer le fichier vide ou échoué
                 # rm -f "$ALIGNMENT_OUTPUT_SAM"
                 continue
             fi
             echo "    Unique BAM created: $UNIQUE_BAM_OUTPUT ($unique_count reads)"

             # --- *** REMPLACEMENT DE AWK PAR BEDTOOLS *** ---
             echo "    Converting unique BAM to BED format using bedtools..."
             # -i : Fichier d'entrée BAM
             # -bed6 : Sortir au format BED6 (chr, start, end, name, score, strand)
             # bedtools gère automatiquement les coordonnées 0-based/1-based et le CIGAR
             bedtools bamtobed -i "$UNIQUE_BAM_OUTPUT" -bed6 > "$BED_OUTPUT"
             # --- *** FIN DU REMPLACEMENT *** ---

             # Vérifier la sortie de bedtools
             if [[ $? -ne 0 || ! -s "$BED_OUTPUT" ]]; then
                 echo "    Error: bedtools bamtobed failed or produced empty output for $UNIQUE_BAM_OUTPUT."
                 rm -f "$BED_OUTPUT" # Supprimer le fichier BED échoué/vide
             else
                 echo "    BED file generated: $BED_OUTPUT"
                 # Optional: Clean up intermediate files (garder le SAM original peut être utile)
                 rm -f "$UNIQUE_BAM_OUTPUT" # On n'a plus besoin du BAM unique
                 # rm -f "$ALIGNMENT_OUTPUT_SAM"
             fi

        done # Fin boucle sur les fichiers fasta

        if ! $processed_file ; then
             echo "  No input FASTA files found matching '*_filtered.fasta' in $chr_dir"
        fi
    fi
done

shopt -u nullglob

echo "Script finished."