import argparse
import os
import sys
import pandas as pd
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import re
import gzip
import random # <-- Ajout de l'import

try:
    import pysam
except ImportError:
    print("ERREUR: Le module 'pysam' est requis mais non trouvé.")
    print("Veuillez l'installer avec : pip install pysam")
    sys.exit(1)

# --- Fonctions Utilitaires ---

def sanitize_filename(name):
    # ... (inchangé) ...
    name = name.replace(" ", "_")
    name = re.sub(r'[^\w\._-]', '', name)
    if not name or name == ".":
        return "unnamed_tag"
    return name

# --- MODIFIÉ : Renvoie (chrom, start, end) ---
def parse_qname_for_g1_coords(qname):
    """
    Parse un identifiant au format 'chromosome:start-end' et retourne (chrom, start, end) tuple.
    Retourne None si le format ne correspond pas. Start/End sont retournés comme entiers (supposés 1-based).
    """
    match = re.match(r'^([^:]+):(\d+)-(\d+)', qname)
    if match:
        chrom, start_str, end_str = match.groups()
        try:
            # Convertir en entier
            start = int(start_str)
            end = int(end_str)
            return chrom, start, end # Retourne les 3 informations
        except ValueError:
            # print(f"AVERTISSEMENT: Erreur de conversion numérique dans QNAME '{qname}'")
            return None # Échec de la conversion
    else:
        # print(f"AVERTISSEMENT: Format QNAME non reconnu pour le parsing des coordonnées G1: '{qname}'")
        return None

def load_tag_intervals(tags_filepath):
    # ... (inchangé par rapport à la version précédente qui fonctionne) ...
    print(f"Chargement des définitions de tags depuis : {tags_filepath}...")
    tag_intervals = defaultdict(IntervalTree)
    try:
        df = pd.read_csv(tags_filepath, sep='\t')
        cols_lower = {c.lower(): c for c in df.columns}
        req_cols_lower = {'tag_name', 'chromosome_1', 'start_1', 'end_1'}
        if not req_cols_lower.issubset(cols_lower.keys()):
             raise ValueError(f"Colonnes manquantes ou mal interprétées dans {tags_filepath}. Attendu (tab-séparé): Tag_name, Chromosome_1, start_1, end_1")

        tag_col = cols_lower['tag_name']
        chr_col = cols_lower['chromosome_1']
        start_col = cols_lower['start_1']
        end_col = cols_lower['end_1']

        count = 0
        for index, row in df.iterrows():
            try:
                chrom_wtj2 = str(row[chr_col])
                start_str = str(row[start_col]).replace(" ", "").replace("\u00A0", "")
                end_str = str(row[end_col]).replace(" ", "").replace("\u00A0", "")
                if not start_str.isdigit() or not end_str.isdigit(): raise ValueError("Coordonnées non numériques")
                start_tag = int(start_str)
                end_tag = int(end_str)
                tag_name = str(row[tag_col])
                if pd.isna(start_tag) or pd.isna(end_tag) or pd.isna(chrom_wtj2) or pd.isna(tag_name): raise ValueError("Valeur NaN")
                tag_intervals[chrom_wtj2].addi(start_tag - 1, end_tag, tag_name) # 0-based [start-1, end)
                count += 1
            except (ValueError, TypeError, KeyError) as e:
                print(f"AVERTISSEMENT: Ignoré ligne invalide {index + 2} dans {tags_filepath}: {row.to_dict()} - Erreur: {e}")
                continue

        print(f"Chargé {count} régions de tags valides pour {len(tag_intervals)} chromosomes WT_J2.")
        if count == 0: print(f"AVERTISSEMENT: Aucune région de tag valide chargée.")
        return dict(tag_intervals)
    except FileNotFoundError:
        print(f"ERREUR: Fichier de définition des tags non trouvé: {tags_filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERREUR lors de la lecture du fichier de tags {tags_filepath}: {e}")
        sys.exit(1)

# --- NOUVELLE FONCTION POUR LIRE LE FASTQ ---
def load_primer_sequences_fastq_by_id(fastq_filepath):
    """
    Charge les séquences depuis un fichier FASTQ dans un dictionnaire.
    La clé du dictionnaire est l'identifiant complet trouvé après '@' (avant espace).
    Cet identifiant DOIT correspondre au QNAME du fichier SAM (ex: 'chromosome:start-end').
    """
    print(f"Chargement des séquences depuis FASTQ: {fastq_filepath} (clé=ID header)...")
    sequences = {}
    record_counter = 0
    sequences_loaded = 0

    open_func = gzip.open if fastq_filepath.endswith(".gz") else open
    try:
        with open_func(fastq_filepath, "rt") as infile:
            while True:
                line1 = infile.readline() # Header line (@...)
                if not line1: break # Fin du fichier
                line2 = infile.readline() # Sequence line
                line3 = infile.readline() # + line
                line4 = infile.readline() # Quality line

                record_counter += 1

                if not (line1 and line2 and line3 and line4):
                    print(f"AVERTISSEMENT: Enregistrement FASTQ incomplet (lu {record_counter}). Arrêt de la lecture.")
                    break

                if line1.startswith('@'):
                    # Clé = identifiant complet après '@' et avant le premier espace
                    header_id = line1[1:].strip().split(maxsplit=1)[0]
                    sequence = line2.strip().upper()

                    if not header_id:
                         print(f"AVERTISSEMENT: En-tête FASTQ vide détecté à l'enregistrement #{record_counter}. Ignoré.")
                         continue # Skip

                    if not sequence:
                         print(f"AVERTISSEMENT: Séquence vide détectée pour l'ID '{header_id}' à l'enregistrement #{record_counter}. Ignoré.")
                         continue # Skip

                    # Vérifier si l'ID existe déjà (peut arriver dans certains FASTQ)
                    # if header_id in sequences:
                    #     print(f"AVERTISSEMENT: ID FASTQ dupliqué '{header_id}' rencontré à l'enregistrement #{record_counter}. La séquence précédente sera écrasée.")
                    sequences[header_id] = sequence
                    sequences_loaded += 1 # Compter seulement si ajouté/écrasé
                else:
                     print(f"AVERTISSEMENT: Ligne ne commençant pas par '@' détectée à la place d'un en-tête FASTQ (enregistrement #{record_counter}). Ignoré.")

        print(f"Terminé la lecture de {record_counter} enregistrements.")
        # Utiliser len(sequences) car sequences_loaded ne compte pas les écrasements de duplicats
        print(f"Chargé {len(sequences)} séquences uniques (basées sur l'ID).")
        if not sequences: print(f"AVERTISSEMENT: Aucune séquence chargée depuis {fastq_filepath}.")
        return sequences
    except FileNotFoundError:
        print(f"ERREUR: Fichier FASTQ des primers non trouvé: {fastq_filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERREUR lors de la lecture du fichier FASTQ {fastq_filepath}: {e}")
        sys.exit(1)
# --- FIN NOUVELLE FONCTION ---

# --- Fonction Principale ---
def main():
    parser = argparse.ArgumentParser(
        description="Assigne des tags aux primers basés sur SAM (WT_J2 coords), filtre pour densité (~2/kb sur G1 coords), et génère FASTA groupés."
    )
    parser.add_argument("input_sam", help="Fichier SAM pré-aligné sur WT_J2. QNAME='chromosome:start-end' (coords G1).")
    parser.add_argument("primer_tags_file", help="Fichier TSV définissant les tags sur WT_J2 (Primer_tag.tsv).")
    parser.add_argument("primer_sequences_fasta", help="Fichier FASTA contenant les séquences. ID='chromosome:start-end'.") # Ou _fastq si besoin
    parser.add_argument("output_base_dir", help="Répertoire de base pour la sortie.")
    parser.add_argument("--target_density", type=float, default=3.0, help="Densité cible en sondes/kb (défaut: 2.0). La sélection se fait à 2 par fenêtre de 1000bp.")
    parser.add_argument("--tag_unassigned", default="Unassigned", help="Tag par défaut (défaut: Unassigned).")

    args = parser.parse_args()
    TARGET_PROBES_PER_WINDOW = 2 # Cible de 2 sondes par fenêtre de 1000bp

    # --- Chargement des données ---
    tag_intervals = load_tag_intervals(args.primer_tags_file)
    primer_sequences = load_primer_sequences_fastq_by_id(args.primer_sequences_fasta) # Adapter si FASTQ

    if not tag_intervals or not primer_sequences:
        print("ERREUR: Chargement des tags ou séquences échoué/vide. Arrêt.")
        sys.exit(1)

    # --- Etape 1: Traitement SAM et Stockage Intermédiaire ---
    print("\n--- Lecture SAM, Assignation Tags et Préparation Filtrage ---")
    # Structure: {target_g1_chr: [(probe_id, g1_start, assigned_tag, sequence)]}
    probes_to_filter = defaultdict(list)
    processed_sam_records = 0
    skipped_qname_parse_error = 0
    skipped_no_seq = 0
    skipped_unmapped_or_secondary = 0
    temp_assigned_count = 0
    temp_unassigned_count = 0

    try:
        with pysam.AlignmentFile(args.input_sam, "r") as samfile:
            for alignment in samfile:
                processed_sam_records += 1
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    skipped_unmapped_or_secondary += 1
                    continue
                probe_id = alignment.query_name
                if not probe_id: continue

                # 1. Parser QNAME pour G1 coords
                g1_coords = parse_qname_for_g1_coords(probe_id)
                if not g1_coords:
                    skipped_qname_parse_error += 1
                    continue
                target_chr, g1_start, g1_end = g1_coords # g1_start est 1-based

                # 2. Récupérer séquence
                sequence = primer_sequences.get(probe_id)
                if not sequence:
                    skipped_no_seq += 1
                    continue

                # 3. Assigner tag basé sur coords WTJ2 (G2)
                assigned_tag = args.tag_unassigned
                chr_wtj2 = alignment.reference_name
                start_wtj2_0based = alignment.reference_start
                end_wtj2_0based_exclusive = alignment.reference_end

                if chr_wtj2 is not None and start_wtj2_0based is not None and end_wtj2_0based_exclusive is not None :
                    midpoint_wtj2_0based = (start_wtj2_0based + end_wtj2_0based_exclusive) // 2
                    if chr_wtj2 in tag_intervals:
                        overlapping_tags = tag_intervals[chr_wtj2].at(midpoint_wtj2_0based)
                        if overlapping_tags:
                            first_interval = min(overlapping_tags, key=lambda x: x.begin)
                            assigned_tag = first_interval.data
                            temp_assigned_count += 1 # Compte avant filtrage
                        else: temp_unassigned_count += 1
                    else: temp_unassigned_count += 1
                else: temp_unassigned_count += 1

                # 4. Stocker pour filtrage
                probes_to_filter[target_chr].append(
                    (probe_id, g1_start, assigned_tag, sequence)
                )
    # ... (Gestion Exceptions SAM) ...
    except FileNotFoundError: print(f"ERREUR: Fichier SAM non trouvé: {args.input_sam}"); sys.exit(1)
    except ValueError as e: print(f"ERREUR: Analyse SAM {args.input_sam}. Valide? Détails: {e}"); sys.exit(1)
    except Exception as e: print(f"ERREUR inattendue SAM: {e}"); import traceback; traceback.print_exc(); sys.exit(1)

    print("Lecture SAM et assignation initiale terminées.")
    print(f"Stats avant filtrage: Traités={processed_sam_records}, Assignés={temp_assigned_count}, NonAssignés={temp_unassigned_count}, SkipQNAME={skipped_qname_parse_error}, SkipNoSeq={skipped_no_seq}, SkipAlign={skipped_unmapped_or_secondary}")

    # --- Etape 2: Filtrage par Densité ---
    print(f"\n--- Filtrage pour densité cible (~{TARGET_PROBES_PER_WINDOW} sondes / kb sur Genome 1) ---")
    # Structure finale: {target_chr: {assigned_tag: [(probe_id, sequence)]}}
    grouped_probes = defaultdict(lambda: defaultdict(list))
    total_probes_before_filtering = 0
    total_probes_after_filtering = 0

    # Itérer sur chaque chromosome cible (G1)
    for target_chr, probe_list in probes_to_filter.items():
        num_probes_chr_before = len(probe_list)
        total_probes_before_filtering += num_probes_chr_before
        print(f"  Filtrage chromosome {target_chr} ({num_probes_chr_before} sondes)...")

        # Grouper les sondes par fenêtre de 1000bp basée sur g1_start (1-based)
        windows = defaultdict(list)
        for probe_data in probe_list:
            # probe_data = (probe_id, g1_start, assigned_tag, sequence)
            g1_start = probe_data[1]
            # Calculer l'index de la fenêtre (0-based)
            window_key = (g1_start - 1) // 1000
            windows[window_key].append(probe_data) # Stocker le tuple complet

        # Sélectionner les sondes dans chaque fenêtre
        selected_probes_for_chr = []
        for window_key, probes_in_window in windows.items():
            if len(probes_in_window) <= TARGET_PROBES_PER_WINDOW:
                # Garder toutes les sondes si <= cible
                selected_probes_for_chr.extend(probes_in_window)
            else:
                # Sélectionner aléatoirement le nombre cible de sondes
                selected = random.sample(probes_in_window, TARGET_PROBES_PER_WINDOW)
                selected_probes_for_chr.extend(selected)

        num_probes_chr_after = len(selected_probes_for_chr)
        print(f"    -> Conservé {num_probes_chr_after} sondes après filtrage.")
        total_probes_after_filtering += num_probes_chr_after

        # Remplir la structure de groupement final avec les sondes sélectionnées
        for selected_data in selected_probes_for_chr:
            probe_id, _, assigned_tag, sequence = selected_data # g1_start n'est plus nécessaire ici
            grouped_probes[target_chr][assigned_tag].append((probe_id, sequence))

    print(f"Filtrage terminé.")
    print(f"Nombre total de sondes avant filtrage: {total_probes_before_filtering}")
    print(f"Nombre total de sondes après filtrage: {total_probes_after_filtering}")


    # --- Etape 3: Écriture des fichiers FASTA (utilise grouped_probes filtré) ---
    print("\n--- Écriture des fichiers FASTA groupés ---")
    # ... (La boucle d'écriture existante est correcte et utilise grouped_probes) ...
    files_written = 0
    total_probes_written = 0
    os.makedirs(args.output_base_dir, exist_ok=True)
    if not grouped_probes: print("AVERTISSEMENT: Aucun groupe de sondes à écrire.")
    for target_chr, tags_dict in grouped_probes.items():
        safe_target_chr = sanitize_filename(target_chr)
        chr_dir = os.path.join(args.output_base_dir, safe_target_chr)
        os.makedirs(chr_dir, exist_ok=True)
        for assigned_tag, probes_list in tags_dict.items():
            if not probes_list: continue
            safe_tag_name = sanitize_filename(assigned_tag)
            output_fasta_path = os.path.join(chr_dir, f"{safe_tag_name}.fasta")
            try:
                with open(output_fasta_path, "w") as outfile:
                    written_in_file = 0
                    for probe_id, sequence in probes_list:
                        outfile.write(f">{probe_id}\n")
                        outfile.write(f"{sequence}\n")
                        total_probes_written += 1
                        written_in_file += 1
                files_written += 1
            except Exception as e:
                 print(f"ERREUR lors de l'écriture du fichier {output_fasta_path}: {e}")

    print(f"Écriture terminée.")
    print(f"Nombre total de fichiers FASTA écrits: {files_written}")
    print(f"Nombre total de sondes écrites (après filtrage): {total_probes_written}")
    print(f"Fichiers générés dans : {args.output_base_dir}")

if __name__ == "__main__":
    import re
    import random # Assurer l'import
    main()