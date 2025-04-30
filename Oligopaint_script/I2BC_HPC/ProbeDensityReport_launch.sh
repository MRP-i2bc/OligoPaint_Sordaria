#!/bin/bash
#SBATCH --time 1000:00:00
#SBATCH --cpus-per-task 32
#SBATCH --mem 100G

eval "$(~/miniconda3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate ~/miniforge3/envs/genome_analysis

cd /data/work/I2BC/antoine.bochet

# --- Définir les chemins et paramètres (plus lisible) ---
REF_GENOME="/data/work/I2BC/antoine.bochet/Sordaria_macrospora/TEST_PIPELINE/assemblies/WT_J2.fasta"
PROBE_DIR_BASE="/data/work/I2BC/antoine.bochet/Output_Tagged_From_SAM_denisty_FINAL_Chr" # Répertoire parent contenant les FASTAs de sondes
OUTPUT_HTML="Probe_density_WT_J2_combined_FINAL.html" # Nouveau nom pour éviter confusion
OUTPUT_DIR="/data/work/I2BC/antoine.bochet/ProbeDensityOutput_final" # Répertoire de sortie dédié pour cette analyse
GFF_ANNOTATIONS="/data/work/I2BC/antoine.bochet/Sordaria_macrospora/CentroReport/WT_J2.gff3"
GENES_LIST="SMAC4_02785,SMAC4_03194,SMAC4_05555,SMAC4_07252,SMAC4_08584,SMAC4_02486,SMAC4_00050,SMAC4_01862,SMAC4_05405,SMAC4_05404,SMAC4_05403,SMAC4_05402,SMAC4_05401"
WINDOW_SIZE=10000
THREADS=32
BLAST_OPTIONS="-perc_identity 100 -qcov_hsp_perc 100"

# --- Rassembler TOUS les fichiers FASTA de sondes ---
# Utilise find pour chercher tous les .fasta dans le répertoire de base et ses sous-dossiers
# Stocke les résultats dans un tableau bash 'probe_files_array'
# L'option -print0 et read -d '' gère les noms de fichiers avec espaces/caractères spéciaux
probe_files_array=()
while IFS= read -r -d '' file; do
    probe_files_array+=("$file")
done < <(find "$PROBE_DIR_BASE" -name '*.fasta' -print0)

# Vérifier si des fichiers ont été trouvés
if [ ${#probe_files_array[@]} -eq 0 ]; then
    echo "Erreur : Aucun fichier FASTA de sonde trouvé dans '$PROBE_DIR_BASE' ou ses sous-répertoires."
    exit 1
fi

echo "Sondes trouvées :"
printf " - %s\n" "${probe_files_array[@]}"
echo "Nombre total de fichiers de sondes : ${#probe_files_array[@]}"

# --- Exécuter ProbeDensityReport.py UNE SEULE FOIS avec tous les fichiers ---
echo "Lancement de ProbeDensityReport.py..."

python3 ProbeDensityReport.py \
    -o "$OUTPUT_HTML" \
    --output_dir "$OUTPUT_DIR" \
    -w "$WINDOW_SIZE" \
    -t "$THREADS" \
    --blast_options "$BLAST_OPTIONS" \
    --gff-annotations "$GFF_ANNOTATIONS" \
    --genes-list "$GENES_LIST" \
    "$REF_GENOME" \
    "${probe_files_array[@]}" # Passe tous les fichiers de sondes trouvés

# Vérifier le code de sortie du script Python
exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo "ProbeDensityReport.py terminé avec succès."
    echo "Rapport généré : $OUTPUT_HTML"
else
    echo "Erreur : ProbeDensityReport.py a échoué avec le code de sortie $exit_code."
fi

exit $exit_code