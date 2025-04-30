import argparse
import os
import shlex # Pour construire/afficher les commandes de manière plus sûre
import subprocess
import math
from collections import Counter
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import concurrent.futures
from functools import partial
import time # To measure execution time

# --- Fonctions d'Analyse ---

def parse_fasta(fasta_file):
    """Lit un fichier FASTA et retourne un dictionnaire {id: sequence}."""
    try:
        sequences = {record.id: str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")}
        print(f"Fichier FASTA lu : {len(sequences)} séquences trouvées.")
        if not sequences:
            raise ValueError("Aucune séquence trouvée dans le fichier FASTA.")
        return sequences
    except FileNotFoundError:
        print(f"Erreur : Le fichier FASTA '{fasta_file}' n'a pas été trouvé.")
        exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier FASTA : {e}")
        exit(1)

# --- Fonctions pour Outils Externes (Placeholders) ---

def run_trf(fasta_file, output_dir="trf_output"):
    """
    Exécute TRF (Tandem Repeats Finder) sur le fichier FASTA entier.
    TRF doit être installé et accessible dans le PATH de l'environnement conda.
    Retourne le chemin complet vers le fichier .dat généré ou None si échec.
    """
    print("\n--- Exécution de TRF ---")
    trf_executable = "trf" # Assumes trf est dans le PATH de l'env conda

    # --- PARAMÈTRES TRF (Ces paramètres sont des exemples, ajustez si besoin) ---
    # Format: match mismatch delta PM PI minscore maxperiod
    trf_params = ["2", "7", "7", "80", "10", "50", "2000"]
    trf_flags = [
        "-d", # Nécessaire pour générer le fichier .dat
        "-h", # Empêche la génération du fichier HTML (souvent inutile)
        "-m"  # Optionnel : génère un fichier masque
    ]
    # --- FIN PARAMÈTRES ---

    # 1. Construire le nom de fichier .dat attendu par TRF
    base_fasta_name = os.path.basename(fasta_file)
    # Le nom est basé sur les paramètres numériques uniquement
    expected_dat_suffix = ".".join(trf_params) + ".dat"
    expected_dat_filename = f"{base_fasta_name}.{expected_dat_suffix}"

    # 2. Préparer les chemins absolus (plus sûr pour subprocess)
    abs_fasta_path = os.path.abspath(fasta_file)
    abs_output_dir = os.path.abspath(output_dir)
    expected_dat_filepath = os.path.join(abs_output_dir, expected_dat_filename)

    # 3. Créer le répertoire de sortie
    try:
        os.makedirs(abs_output_dir, exist_ok=True)
    except OSError as e:
        print(f"Erreur lors de la création du répertoire de sortie '{abs_output_dir}': {e}")
        return None

    # 4. Construire la commande complète
    command = [trf_executable, abs_fasta_path] + trf_params + trf_flags
    print(f"Répertoire de travail pour TRF : {abs_output_dir}")
    # Utiliser shlex.quote pour afficher la commande de manière interprétable par le shell
    print(f"Commande TRF : {' '.join(shlex.quote(c) for c in command)}")

    # 5. Exécuter TRF via subprocess
    try:
        # Exécute la commande DEPUIS le répertoire de sortie spécifié
        # TRF créera ses fichiers directement là
        result = subprocess.run(command,
                                cwd=abs_output_dir, # Important: définit le répertoire d'exécution
                                capture_output=True, # Capturer stdout/stderr
                                text=True,           # Décoder la sortie en texte
                                encoding='utf-8',    # Spécifier l'encodage
                                errors='replace')    # Gérer les erreurs d'encodage potentielles

        # Afficher la sortie standard (peut contenir des infos utiles)
        print("--- Sortie Standard TRF ---")
        print(result.stdout if result.stdout else "(Aucune)")
        print("--------------------------")
        # Afficher la sortie d'erreur même si succès (peut contenir des warnings)
        if result.stderr:
            print("--- Sortie Erreur TRF (Warnings possibles) ---")
            print(result.stderr)
            print("--------------------------------------------")

        # 6. Vérifier si le fichier .dat attendu existe
        if os.path.exists(expected_dat_filepath):
            print(f"Exécution de TRF réussie.")
            print(f"Fichier de sortie TRF trouvé : {expected_dat_filepath}")
            return expected_dat_filepath # Retourner le chemin complet du fichier .dat
        else:
            # Cas étrange où TRF réussit mais le fichier n'est pas là
            print(f"Erreur : TRF s'est terminé avec succès (code 0), mais le fichier .dat attendu '{expected_dat_filepath}' est introuvable.")
            print("Vérifiez les paramètres TRF, les permissions d'écriture et les conventions de nommage de TRF.")
            return None

    # 7. Gérer les erreurs d'exécution
    except FileNotFoundError:
        # L'exécutable TRF n'a pas été trouvé dans le PATH
        print(f"Erreur CRITIQUE : L'exécutable '{trf_executable}' est introuvable.")
        print("Vérifiez que TRF est bien installé dans l'environnement Conda et que cet environnement est activé.")
        return None
    except subprocess.CalledProcessError as e:
        # TRF a retourné un code d'erreur (non-zéro)
        print(f"Erreur CRITIQUE : L'exécution de TRF a échoué avec le code de sortie {e.returncode}.")
        print(f"Commande exécutée : {' '.join(shlex.quote(c) for c in e.cmd)}")
        print("--- Sortie Standard (TRF Échec) ---")
        print(e.stdout if e.stdout else "(Aucune)")
        print("--- Sortie Erreur (TRF Échec) ---")
        print(e.stderr if e.stderr else "(Aucune)")
        print("----------------------------------")
        return None
    except Exception as e:
        # Autres erreurs potentielles (permissions, etc.)
        print(f"Erreur inattendue lors de l'exécution de TRF : {e}")
        return None

def parse_trf_output(dat_file_path, sequence_ids):
    """
    Parse le fichier .dat unique généré par TRF pour extraire les positions
    des répétitions pour chaque séquence originale.

    Args:
        dat_file_path (str): Chemin vers le fichier .dat généré par run_trf.
                               Peut être None si run_trf a échoué.
        sequence_ids (list): Liste ou set des IDs de séquence originaux du fichier FASTA.
                               Nécessaire pour savoir quelles séquences chercher.

    Returns:
        dict: Dictionnaire {seq_id: list_of_tuples[(start, end)]}.
              Retourne un dict vide si le fichier est introuvable, vide ou en cas d'erreur.
    """
    print(f"\n--- Parsing du fichier de résultats TRF : {dat_file_path} ---")
    if not dat_file_path or not os.path.exists(dat_file_path):
        print("Avertissement : Fichier .dat de TRF non fourni ou introuvable. Aucune donnée TRF ne sera chargée.")
        return {seq_id: [] for seq_id in sequence_ids} # Retourne vide pour toutes les séquences

    # Initialiser le dictionnaire des résultats avec une liste vide pour chaque ID attendu
    all_positions = {seq_id: [] for seq_id in sequence_ids}
    current_seq_id = None
    sequences_found_in_dat = set()

    try:
        with open(dat_file_path, 'r', encoding='utf-8', errors='replace') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue # Ignorer lignes vides

                # Détecter la ligne indiquant une nouvelle séquence
                if line.startswith("Sequence:"):
                    parts = line.split(maxsplit=2) # Sépare "Sequence:", "seqid", "description..."
                    if len(parts) >= 2:
                        trf_seq_id = parts[1]
                        # Vérifier si cet ID est l'un de ceux que nous attendons
                        if trf_seq_id in all_positions:
                            current_seq_id = trf_seq_id
                            sequences_found_in_dat.add(current_seq_id)
                            # print(f"  Début des données pour la séquence : {current_seq_id}") # Debug
                        else:
                            # Si l'ID du fichier TRF n'est pas dans notre liste initiale (ne devrait pas arriver si FASTA est cohérent)
                            print(f"Avertissement Ligne {line_num}: ID de séquence '{trf_seq_id}' trouvé dans '{dat_file_path}' mais pas dans les IDs attendus. Section ignorée.")
                            current_seq_id = None # Ignorer les données suivantes jusqu'à la prochaine ligne Sequence: valide
                    else:
                         # Ligne commençant par Sequence: mais mal formée
                         print(f"Avertissement Ligne {line_num}: Ligne 'Sequence:' mal formée ignorée: '{line}'")
                         current_seq_id = None

                # Si on est dans une section de séquence valide et la ligne commence par un chiffre (données)
                elif current_seq_id and line[0].isdigit():
                    parts = line.split()
                    # Format TRF typique: start end period_size copy_num ...
                    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                        try:
                            # Les coordonnées TRF sont 1-based, inclusives
                            start_1based = int(parts[0])
                            end_1based = int(parts[1])

                            # Convertir en 0-based, semi-ouvert [start, end) pour Python
                            start_0based = start_1based - 1
                            end_0based_exclusive = end_1based # L'indice de fin est déjà exclusif

                            if start_0based < end_0based_exclusive: # Vérification de cohérence
                                all_positions[current_seq_id].append((start_0based, end_0based_exclusive))
                            else:
                                print(f"Avertissement Ligne {line_num}: Position TRF invalide ignorée (start >= end après conversion 0-based): '{line}' pour {current_seq_id}")
                        except ValueError:
                            print(f"Avertissement Ligne {line_num}: Ligne de données TRF contenant des non-nombres ignorée : '{line}' pour {current_seq_id}")
                    # else: Ignorer les autres lignes dans une section (commentaires, etc.)

        # Vérifier si des séquences attendues n'ont pas été trouvées dans le fichier .dat
        missing_seqs = set(sequence_ids) - sequences_found_in_dat
        if missing_seqs:
            print(f"Avertissement : Les IDs de séquence suivants du FASTA n'ont pas été trouvés dans le fichier de sortie TRF '{dat_file_path}': {', '.join(missing_seqs)}")

        # Compter et afficher le résumé
        total_found = sum(len(v) for v in all_positions.values())
        print(f"Parsing TRF terminé. {total_found} répétitions trouvées au total pour {len(sequences_found_in_dat)} séquences.")
        # for seq_id in sorted(sequences_found_in_dat): # Optionnel: afficher le compte par séquence
        #     print(f"  - {seq_id}: {len(all_positions[seq_id])} répétitions")

        return all_positions

    except Exception as e:
        print(f"Erreur majeure et inattendue lors du parsing du fichier {dat_file_path}: {e}")
        # Retourner un dictionnaire vide pour éviter de planter la suite
        return {seq_id: [] for seq_id in sequence_ids}

# --- Fonctions K-mer (Parallélisées) ---

def count_kmers_for_sequence(sequence_tuple, k):
    """Compte les k-mers pour une seule séquence (chromosome/contig). Worker Function."""
    seq_id, seq = sequence_tuple
    local_counts = Counter()
    if len(seq) >= k:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            # Gérer les Ns: ignorer les k-mers contenant 'N'
            if 'N' not in kmer:
                local_counts[kmer] += 1
    # print(f"  Comptage terminé pour {seq_id}") # Peut être très verbeux
    return local_counts

def count_genome_kmers_parallel(sequences, k, num_workers):
    """Compte les k-mers sur tout le génome en parallèle."""
    print(f"\n--- Comptage parallèle des {k}-mers ({num_workers} workers) ---")
    start_time = time.time()
    genome_counts = Counter()
    # 'partial' fixe l'argument 'k' pour la fonction worker
    worker_func = partial(count_kmers_for_sequence, k=k)

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Utiliser sequences.items() pour obtenir des tuples (id, seq)
        future_to_seqid = {executor.submit(worker_func, item): item[0] for item in sequences.items()}
        
        for future in concurrent.futures.as_completed(future_to_seqid):
            seq_id = future_to_seqid[future]
            try:
                local_counts = future.result()
                genome_counts.update(local_counts)
                # print(f"  Résultats fusionnés pour {seq_id}") # Verbeux
            except Exception as exc:
                print(f'Erreur lors du comptage pour {seq_id}: {exc}')

    end_time = time.time()
    print(f"Comptage des k-mers terminé en {end_time - start_time:.2f} secondes.")
    print(f"Total de {len(genome_counts):,} {k}-mers distincts trouvés.")
    return genome_counts

# --- Fonctions de Calcul de Densité (Worker pour Parallélisation) ---

def calculate_density(seq_len, window_size, feature_positions=None, kmer_starts=None):
    """
    Calcule la densité (couverture ou comptage normalisé) par fenêtre.
    Retourne: (window_centers, densities)
    """
    if seq_len == 0: return [], []
    window_starts = list(range(0, seq_len, window_size))
    num_windows = len(window_starts)
    densities = [0.0] * num_windows
    window_centers = [s + min(window_size / 2, (seq_len - s)/2) for s in window_starts] # Centre de la fenêtre

    if feature_positions is not None: # Calcul basé sur la couverture (ex: TRF)
        # Trier les features par début peut aider mais n'est pas crucial ici
        for i, win_start in enumerate(window_starts):
            win_end = min(win_start + window_size, seq_len)
            win_len = win_end - win_start
            if win_len <= 0: continue

            covered_length = 0
            for feat_start, feat_end in feature_positions:
                # Calculer l'overlap entre la feature [feat_start, feat_end) et la fenêtre [win_start, win_end)
                overlap_start = max(win_start, feat_start)
                overlap_end = min(win_end, feat_end)
                covered_length += max(0, overlap_end - overlap_start)
            densities[i] = covered_length / win_len

    elif kmer_starts is not None: # Calcul basé sur le comptage (ex: K-mers)
        # Compter combien de k-mers débutent dans chaque fenêtre
        kmer_indices = [pos // window_size for pos in kmer_starts if pos // window_size < num_windows]
        counts_per_window = Counter(kmer_indices)
        for i, win_start in enumerate(window_starts):
            win_end = min(win_start + window_size, seq_len)
            win_len = win_end - win_start
            if win_len <= 0: continue
            # Normaliser par la taille de la fenêtre (k-mers par bp)
            # Alternative: utiliser le compte brut counts_per_window.get(i, 0)
            densities[i] = counts_per_window.get(i, 0) / win_len

    return window_centers, densities


def calculate_all_densities_for_chromosome(args_tuple):
    """
    Fonction worker pour calculer les densités k-mers, TRF (optionnel)
    et RepeatMasker (TEs, RetroTEs, Satellites) pour un chromosome.
    """
    # Dépaqueter les arguments (SSR enlevé)
    (seq_id, seq, k, unique_kmers_set, repeated_kmers_set,
     repeat_threshold, window_size, trf_positions_for_chr,
     rm_te_positions_for_chr, rm_retro_positions_for_chr,
     rm_sat_positions_for_chr) = args_tuple # <--- Signature mise à jour

    # print(f"  Calcul des densités pour {seq_id}...") # Verbeux
    seq_len = len(seq)
    densities_for_chr = {}
    centers = [] # Garder une référence aux centres de fenêtres

    # --- Densité K-mers Uniques ---
    unique_kmer_starts = [i for i in range(seq_len - k + 1) if seq[i:i+k] in unique_kmers_set]
    centers, unique_density = calculate_density(seq_len, window_size, kmer_starts=unique_kmer_starts)
    densities_for_chr[f'Kmers Uniques (k={k})'] = (centers, unique_density)

    # --- Densité K-mers Répétés ---
    repeated_kmer_starts = [i for i in range(seq_len - k + 1) if seq[i:i+k] in repeated_kmers_set]
    if centers: # S'assurer qu'on a des centres de fenêtres
        _, repeated_density = calculate_density(seq_len, window_size, kmer_starts=repeated_kmer_starts)
    else:
        repeated_density = []
    densities_for_chr[f'Kmers Répétés (k={k}, seuil≥{repeat_threshold})'] = (centers, repeated_density)

    # --- Densité TRF (Optionnel) ---
    # Note: Le nom de trace sera créé même si trf_positions_for_chr est vide,
    # result sera (centers, [0.0]*len(centers))
    if centers:
        _, trf_density = calculate_density(seq_len, window_size, feature_positions=trf_positions_for_chr)
    else:
        trf_density = []
    densities_for_chr['Répétitions Tandem (TRF)'] = (centers, trf_density)


    # --- Densité RepeatMasker : Transposons (Total) ---
    if centers:
        _, rm_te_density = calculate_density(seq_len, window_size, feature_positions=rm_te_positions_for_chr)
    else:
        rm_te_density = []
    densities_for_chr['Transposons (RepeatMasker)'] = (centers, rm_te_density)

    # --- Densité RepeatMasker : Rétrotransposons ---
    if centers:
        _, rm_retro_density = calculate_density(seq_len, window_size, feature_positions=rm_retro_positions_for_chr)
    else:
        rm_retro_density = []
    densities_for_chr['Rétrotransposons (RM)'] = (centers, rm_retro_density)

    # --- Densité RepeatMasker : Satellites ---
    if centers:
        _, rm_sat_density = calculate_density(seq_len, window_size, feature_positions=rm_sat_positions_for_chr)
    else:
        rm_sat_density = []
    densities_for_chr['Satellites (RM)'] = (centers, rm_sat_density)

    # PAS DE CALCUL SSR ICI

    # Retourner les données nécessaires pour le graphique
    # print(f"  Densités calculées pour {seq_id}.") # Verbeux
    return seq_id, seq_len, centers, densities_for_chr


# Potentiellement utile pour parser le GFF : from BCBio import GFF # pip install bcbio-gff
# Mais on peut faire un parsing simple pour ce besoin spécifique.

def run_repeatmasker(fasta_file, output_dir="repeatmasker_output", threads=1, lib_path=None): # lib_path est maintenant utilisé
    """
    Exécute RepeatMasker sur le fichier FASTA.
    RepeatMasker et ses dépendances doivent être installés (conda recommandé).
    Retourne le chemin vers le fichier GFF de sortie ou None si échec.
    """
    print("\n--- Exécution de RepeatMasker ---")
    rm_executable = "RepeatMasker"

    abs_fasta_path = os.path.abspath(fasta_file)
    abs_output_dir = os.path.abspath(output_dir)
    # Nom du fichier GFF de sortie peut varier, on vérifiera .out.gff
    gff_expected_name = os.path.basename(fasta_file) + ".out.gff"
    gff_expected_path = os.path.join(abs_output_dir, gff_expected_name)

    command = [
        rm_executable,
        "-pa", str(threads), # Nombre de CPUs
        "-gff",             # Sortie GFF
        "-dir", abs_output_dir # Répertoire de sortie
    ]

    # --- Utilisation de la bibliothèque fournie ---
    if lib_path:
        abs_lib_path = os.path.abspath(lib_path)
        if not os.path.exists(abs_lib_path):
             print(f"Erreur CRITIQUE : Le fichier bibliothèque spécifié '{abs_lib_path}' est introuvable.")
             return None
        command.extend(["-lib", abs_lib_path]) # Utilise l'argument de la fonction
        print(f"Utilisation de la bibliothèque de répétitions : {abs_lib_path}")
    else:
        # Si aucune lib n'est fournie, RepeatMasker utilisera sa config par défaut (souvent vertébrés)
        # Il faudrait peut-être une erreur ici si on s'attendait à une lib custom.
        print("Avertissement : Aucune bibliothèque custom (-lib) fournie. RepeatMasker utilisera sa configuration par défaut.")
        # Alternativement, lever une erreur si une lib était attendue :
        # print("Erreur CRITIQUE : L'exécution de RepeatMasker a été demandée sans fournir de bibliothèque custom via --repeatmasker-lib.")
        # return None

    command.append(abs_fasta_path) # Fichier d'entrée à la fin

    print(f"Répertoire de travail pour RepeatMasker : {abs_output_dir}")
    print(f"Commande RepeatMasker : {' '.join(shlex.quote(c) for c in command)}")

    try:
        os.makedirs(abs_output_dir, exist_ok=True)
        result = subprocess.run(command,
                                capture_output=True,
                                text=True,
                                check=True, # Laisser check=True, RM doit finir proprement
                                encoding='utf-8',
                                errors='replace')

        print("--- Sortie Standard RepeatMasker (Résumé Fin) ---")
        print(result.stdout[-1500:] if result.stdout else "(Vide)") # Afficher la fin
        print("-------------------------------------------")
        if result.stderr:
             print("--- Sortie Erreur RepeatMasker ---")
             print(result.stderr) # Afficher si non vide
             print("--------------------------------")

        # Vérifier l'existence du fichier GFF (nommé .out.gff par RM)
        if os.path.exists(gff_expected_path):
            print(f"Exécution de RepeatMasker réussie.")
            print(f"Fichier GFF trouvé : {gff_expected_path}")
            return gff_expected_path
        else:
            print(f"Erreur : RepeatMasker s'est terminé (code 0), mais le fichier GFF attendu '{gff_expected_path}' est introuvable.")
            print(f"Vérifiez le contenu du répertoire de sortie : {abs_output_dir}")
            # Essayer de trouver un autre .gff au cas où ?
            gff_files_in_dir = [f for f in os.listdir(abs_output_dir) if f.endswith('.gff')]
            if gff_files_in_dir:
                print(f"Fichiers GFF trouvés dans le répertoire : {gff_files_in_dir}")
                # On pourrait retourner le premier trouvé, mais c'est risqué
            return None

    except FileNotFoundError:
        print(f"Erreur CRITIQUE : L'exécutable '{rm_executable}' est introuvable.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Erreur CRITIQUE : L'exécution de RepeatMasker a échoué (code {e.returncode}).")
        print(f"Commande : {' '.join(shlex.quote(c) for c in e.cmd)}")
        print("--- Sortie Standard (Échec RM) ---")
        print(e.stdout)
        print("--- Sortie Erreur (Échec RM) ---")
        print(e.stderr)
        print("--------------------------------")
        return None
    except Exception as e:
        print(f"Erreur inattendue lors de l'exécution de RepeatMasker : {e}")
        return None


def parse_repeatmasker_gff(gff_file, sequence_ids):
    """
    Parse le fichier GFF de RepeatMasker pour classifier les répétitions.
    Retourne un dictionnaire contenant les positions pour différentes classes.
    """
    print(f"\n--- Parsing du fichier GFF RepeatMasker : {gff_file} ---")
    results = {
        'transposons': {seq_id: [] for seq_id in sequence_ids},
        'retrotransposons': {seq_id: [] for seq_id in sequence_ids},
        'satellites': {seq_id: [] for seq_id in sequence_ids}
        # On pourrait ajouter 'dna_transposons' si besoin
    }
    counts = {'transposons': 0, 'retrotransposons': 0, 'satellites': 0}

    if not gff_file or not os.path.exists(gff_file):
        print("Avertissement : Fichier GFF de RepeatMasker non fourni ou introuvable.")
        return results

    try:
        with open(gff_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue # Ignorer commentaires et lignes vides

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    # print(f"Avertissement Ligne {line_num}: Ligne GFF incomplète ignorée.")
                    continue

                seq_id = fields[0]
                source = fields[1] # Souvent 'RepeatMasker'
                feature_type = fields[2] # Ex: 'similarity', 'match', 'repeat'
                start_str = fields[3]
                end_str = fields[4]
                attributes_str = fields[8]

                # Ignorer si ce n'est pas une feature de type 'repeat' ou similaire trouvée par RM
                # (Certains GFF RM contiennent aussi 'gene', 'exon' etc. si on masque un génome annoté)
                # On se concentre sur les hits de similarité avec la bibliothèque de répétitions.
                # Le type exact peut varier, on se base surtout sur les attributs.

                if seq_id not in results['transposons']:
                    # print(f"Avertissement Ligne {line_num}: ID de séquence '{seq_id}' trouvé dans le GFF mais pas dans le FASTA original. Ignoré.")
                    continue

                try:
                    # GFF est 1-based, inclusif. Convertir en 0-based, semi-ouvert [start, end)
                    start_0based = int(start_str) - 1
                    end_0based_exclusive = int(end_str)
                    if start_0based >= end_0based_exclusive:
                         # print(f"Avertissement Ligne {line_num}: Position GFF invalide (start >= end): {line.strip()}")
                         continue
                except ValueError:
                    # print(f"Avertissement Ligne {line_num}: Coordonnées GFF invalides : {line.strip()}")
                    continue

                # --- Classification basée sur les attributs ---
                # Format des attributs: "Target=Motif:Classe/Famille Start End Score;..."
                # Exemple: Target "Motif:LTR/Gypsy" 1 150 550; ID "RM_1";
                is_retro = False
                is_satellite = False
                repeat_class_family = ""

                # Extraction simple de la classification
                if "Target" in attributes_str:
                    target_part = attributes_str.split("Target=")[1].split(";")[0].strip()
                    # Enlever les guillemets s'ils existent
                    if target_part.startswith('"') and target_part.endswith('"'):
                        target_part = target_part[1:-1]
                    # Format typique: "Motif:Classe/Famille" ou "Simple_repeat" ou "Satellite"
                    if ":" in target_part:
                        repeat_class_family = target_part.split(":")[-1] # Ex: "LTR/Gypsy", "LINE/L1", "Satellite", "Simple_repeat"
                    else:
                        repeat_class_family = target_part # Ex: "Satellite", "(AT)n"

                    # Classer comme Rétrotransposon
                    if any(rt in repeat_class_family for rt in ["LTR", "LINE", "SINE", "Retroposon"]):
                        is_retro = True
                    # Classer comme Satellite
                    elif "Satellite" in repeat_class_family:
                        is_satellite = True

                # Ajouter aux listes appropriées
                # Tous les hits de RM sont considérés comme des 'transposons' au sens large (éléments répétitifs)
                coords = (start_0based, end_0based_exclusive)
                results['transposons'][seq_id].append(coords)
                counts['transposons'] += 1

                if is_retro:
                    results['retrotransposons'][seq_id].append(coords)
                    counts['retrotransposons'] += 1

                if is_satellite:
                    results['satellites'][seq_id].append(coords)
                    counts['satellites'] += 1

        print(f"Parsing GFF RepeatMasker terminé.")
        print(f"  - Transposons (total répétitions RM): {counts['transposons']:,}")
        print(f"  - Rétrotransposons (LTR/LINE/SINE...): {counts['retrotransposons']:,}")
        print(f"  - Satellites (classés par RM): {counts['satellites']:,}")
        return results

    except Exception as e:
        print(f"Erreur majeure lors du parsing du fichier GFF {gff_file}: {e}")
        # Retourner un dictionnaire vide pour ne pas planter
        return {cat: {seq_id: [] for seq_id in sequence_ids} for cat in results}




# --- Fonctions de Plotting et Rapport ---

def create_plots(seq_id, seq_len, window_centers, densities_dict, window_size):
    """Crée un graphique Plotly pour un chromosome avec plusieurs traces de densité."""
    fig = go.Figure()
    if not window_centers: # Si le chromosome est trop court pour une fenêtre
        print(f"Avertissement : Chromosome {seq_id} trop court (longueur {seq_len}) pour la taille de fenêtre {window_size}, graphique vide généré.")
        fig.update_layout(title=f"Densité des Séquences - Chromosome {seq_id} (Trop court)")
        return fig

    plot_added = False
    for name, data in densities_dict.items():
        if data: # Vérifier si l'analyse a produit des données
            centers, values = data
            if centers and values: # S'assurer qu'il y a des données à tracer
                fig.add_trace(go.Scatter(x=centers, y=values, mode='lines', name=name))
                plot_added = True

    if not plot_added:
         print(f"Avertissement : Aucune donnée de densité à tracer pour {seq_id}.")
         # Ajouter une annotation pour indiquer qu'il n'y a pas de données
         fig.add_annotation(x=seq_len/2, y=0, text="Aucune donnée de densité calculée", showarrow=False)


    fig.update_layout(
        title=f"Densité des Séquences - Chromosome {seq_id} (Longueur: {seq_len:,} bp)",
        xaxis_title=f"Position sur le chromosome (Fenêtres de {window_size:,} bp)",
        yaxis_title="Densité (varie selon la métrique)",
        xaxis=dict(range=[0, seq_len]),
        yaxis=dict(rangemode='tozero'), # Assure que l'axe Y commence à 0
        hovermode='x unified' # Améliore l'affichage au survol
    )
    return fig

def generate_html_report(figures_dict, output_file):
    """Génère le rapport HTML avec les graphiques Plotly."""
    print(f"\n--- Génération du rapport HTML : {output_file} ---")
    start_time = time.time()
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("<!DOCTYPE html>\n<html>\n<head>\n")
            f.write("<meta charset='utf-8' />\n")
            f.write("<title>Rapport Densité Centromérique</title>\n")
            # Inclure Plotly.js depuis CDN pour l'interactivité
            f.write("<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n")
            f.write("<style> body { font-family: sans-serif; margin: 20px; } h1, h2 { color: #333; } hr { border: 0; height: 1px; background: #ccc; } .plot-container { margin-bottom: 30px; border: 1px solid #eee; padding: 10px; box-shadow: 2px 2px 5px #eee; } </style>\n")
            f.write("</head>\n<body>\n")
            f.write("<h1>Analyse de Densité pour Localisation Centromérique</h1>\n")

            # Trier les chromosomes/contigs (essayer un tri "naturel" si possible)
            try:
                # Tente de trier numériquement si les ID ressemblent à des nombres (ex: 'chr1', 'chr10')
                sorted_ids = sorted(figures_dict.keys(), key=lambda x: int(''.join(filter(str.isdigit, x))) if any(char.isdigit() for char in x) else float('inf'))
            except:
                # Tri alphabétique simple en cas d'échec du tri numérique
                sorted_ids = sorted(figures_dict.keys())

            for seq_id in sorted_ids:
                fig = figures_dict[seq_id]
                f.write(f"<div class='plot-container'>\n")
                f.write(f"<h2>Chromosome / Contig: {seq_id}</h2>\n")
                # Exporter le graphique en HTML (div uniquement) et l'intégrer
                plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
                f.write(plot_html)
                f.write("</div>\n")
                # f.write("<hr>\n") # Séparateur peut-être redondant avec le container

            f.write("</body>\n</html>\n")
        end_time = time.time()
        print(f"Rapport généré en {end_time - start_time:.2f} secondes.")
    except Exception as e:
        print(f"Erreur lors de la génération du rapport HTML : {e}")

# --- Fonction Principale ---

import argparse
import os
import subprocess
import math
from collections import Counter
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import concurrent.futures
from functools import partial
import time
import shlex

# ASSUMER QUE LES FONCTIONS SUIVANTES SONT DÉFINIES PRÉCÉDEMMENT:
# parse_fasta, run_trf, parse_trf_output,
# count_kmers_for_sequence, count_genome_kmers_parallel,
# run_repeatmasker, parse_repeatmasker_gff,
# calculate_density, calculate_all_densities_for_chromosome (version ci-dessus),
# create_plots, generate_html_report

def main():
    parser = argparse.ArgumentParser(
        description="Analyse la densité des k-mers, répétitions TRF (optionnel) et répétitions RepeatMasker pour aider à identifier les régions centromériques.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # --- Arguments Généraux ---
    parser.add_argument("fasta_file", help="Fichier FASTA du génome.")
    parser.add_argument("-o", "--output_html", default="centromere_density_report.html", help="Nom du fichier HTML de sortie.")
    parser.add_argument("-w", "--window_size", type=int, default=50000, help="Taille des fenêtres d'analyse (en bp).")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Nombre de processus workers (0 = tous les cœurs).")

    # --- Arguments K-mer ---
    parser.add_argument("-k", type=int, default=9, help="Taille des k-mers pour l'analyse.")
    parser.add_argument("-r", "--repeat_threshold", type=int, default=5, help="Seuil minimum d'occurrence pour qu'un k-mer soit 'répété'.")

    # --- Arguments TRF (Optionnel) ---
    parser.add_argument("--run-trf", action='store_true', help="Exécuter TRF (Tandem Repeats Finder).")
    parser.add_argument("--provide-trf-dat", metavar='FILE', help="Fournir un fichier .dat TRF pré-calculé au lieu d'exécuter TRF.")
    parser.add_argument("--output-dir-trf", default="trf_analysis_output", help="Répertoire de sortie pour TRF.")

    # --- Arguments RepeatMasker (Optionnel) ---
    parser.add_argument("--run-repeatmasker", action='store_true', help="Exécuter RepeatMasker (peut être long!).")
    parser.add_argument("--repeatmasker-lib", metavar='FILE', help="Chemin vers une bibliothèque de répétitions custom FASTA pour RepeatMasker (requis si --run-repeatmasker est utilisé sans --provide-rm-gff).") # Précision ajoutée
    parser.add_argument("--provide-rm-gff", metavar='FILE', help="Fournir un fichier GFF RepeatMasker pré-calculé.")
    parser.add_argument("--output-dir-rm", default="repeatmasker_output", help="Répertoire de sortie pour RepeatMasker.")

    args = parser.parse_args()

    # --- Vérification des arguments pour RepeatMasker ---
    if args.run_repeatmasker and not args.provide_rm_gff and not args.repeatmasker_lib:
        parser.error("L'option --run-repeatmasker requiert soit --repeatmasker-lib pour spécifier une bibliothèque, soit --provide-rm-gff pour utiliser des résultats existants.")
    if args.repeatmasker_lib and not os.path.exists(args.repeatmasker_lib):
         parser.error(f"Le fichier bibliothèque RepeatMasker spécifié (--repeatmasker-lib) est introuvable : {args.repeatmasker_lib}")
    # --- Fin Vérification ---


    # --- Début de l'exécution ---
    script_start_time = time.time()
    print("--- Initialisation de l'Analyse ---")
    # (Affichage des paramètres...)
    num_workers = max(1, os.cpu_count() if args.threads <= 0 else min(args.threads, os.cpu_count()))
    print(f"Utilisation de {num_workers} workers pour les tâches parallélisables.")

    # 1. Lire le génome
    sequences = parse_fasta(args.fasta_file)
    sequence_ids_list = list(sequences.keys())

    # 2. Analyse K-mers
    # (Code inchangé...)
    genome_kmer_counts = count_genome_kmers_parallel(sequences, args.k, num_workers)
    unique_kmers_set = {kmer for kmer, count in genome_kmer_counts.items() if count == 1}
    repeated_kmers_set = {kmer for kmer, count in genome_kmer_counts.items() if count >= args.repeat_threshold}
    print(f"{len(unique_kmers_set):,} k-mers uniques identifiés.")
    print(f"{len(repeated_kmers_set):,} k-mers répétés (≥{args.repeat_threshold}) identifiés.")
    del genome_kmer_counts # Libérer mémoire


    # 3. Analyses des caractéristiques répétées
    trf_positions = {seq_id: [] for seq_id in sequence_ids_list}
    rm_results = {'transposons': {}, 'retrotransposons': {}, 'satellites': {}} # Initialiser vide

    # --- TRF ---
    # (Code inchangé...)
    trf_dat_file_path = args.provide_trf_dat
    if not trf_dat_file_path and args.run_trf:
        trf_dat_file_path = run_trf(args.fasta_file, output_dir=args.output_dir_trf)
    if trf_dat_file_path:
        trf_positions = parse_trf_output(trf_dat_file_path, sequence_ids_list)
    # (Gestion des messages inchangée...)


    # --- RepeatMasker ---
    rm_gff_file_path = args.provide_rm_gff
    if not rm_gff_file_path and args.run_repeatmasker:
         print("\nNOTE: L'exécution de RepeatMasker peut être très longue.")
         # Passe la valeur de args.repeatmasker_lib à la fonction
         rm_gff_file_path = run_repeatmasker(args.fasta_file,
                                             output_dir=args.output_dir_rm,
                                             threads=num_workers,
                                             lib_path=args.repeatmasker_lib) # <--- Utilisation correcte
    if rm_gff_file_path:
        rm_results = parse_repeatmasker_gff(rm_gff_file_path, sequence_ids_list)
    # (Gestion des messages et structure vide inchangée...)


    # 4. Préparation des données pour le calcul parallèle des densités
    # (Code inchangé...)
    print(f"\n--- Préparation du calcul parallèle des densités ({num_workers} workers) ---")
    tasks_args = []
    for seq_id, seq in sequences.items():
        trf_data = trf_positions.get(seq_id, [])
        rm_te_data = rm_results.get('transposons', {}).get(seq_id, [])
        rm_retro_data = rm_results.get('retrotransposons', {}).get(seq_id, [])
        rm_sat_data = rm_results.get('satellites', {}).get(seq_id, [])
        tasks_args.append((
            seq_id, seq, args.k, unique_kmers_set, repeated_kmers_set,
            args.repeat_threshold, args.window_size,
            trf_data, rm_te_data, rm_retro_data, rm_sat_data
        ))


    # 5. Exécution du calcul parallèle des densités
    # (Code inchangé...)
    print(f"\n--- Lancement du calcul parallèle des densités ---")
    start_time_density = time.time()
    all_plot_data = {}
    processed_count = 0
    total_seqs = len(sequences)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_seqid = {executor.submit(calculate_all_densities_for_chromosome, args_tuple): args_tuple[0] for args_tuple in tasks_args}
        for future in concurrent.futures.as_completed(future_to_seqid):
            seq_id_res = future_to_seqid[future]
            try:
                seq_id_res, seq_len, centers, densities_dict = future.result()
                all_plot_data[seq_id_res] = (seq_len, centers, densities_dict)
                processed_count += 1
                print(f"  Densités calculées pour {seq_id_res} ({processed_count}/{total_seqs})", end='\r')
            except Exception as exc:
                print(f"\nErreur lors du calcul de densité pour {seq_id_res}: {exc}")
    print(f"\nCalcul des densités terminé en {time.time() - start_time_density:.2f} secondes.")


    # 6. Génération des Graphiques (Séquentiel)
    # (Code inchangé...)
    print("\n--- Génération des graphiques ---")
    all_figures = {}
    start_time_plots = time.time()
    # S'assurer qu'il y a des données à tracer avant de générer le rapport
    if not all_plot_data:
         print("Erreur CRITIQUE: Aucune donnée de densité n'a été calculée pour aucun chromosome. Impossible de générer les graphiques.")
    else:
        for seq_id, plot_data_tuple in all_plot_data.items():
            seq_len, centers, densities_dict = plot_data_tuple
            if not centers: # Vérifier si les centres de fenêtre existent
                print(f"Avertissement: Aucune donnée de densité (pas de fenêtres/centres) pour {seq_id}. Graphique vide.")
            # Créer le graphique même si les densités sont nulles, mais les centres doivent exister
            all_figures[seq_id] = create_plots(seq_id, seq_len, centers, densities_dict, args.window_size)
        end_time_plots = time.time()
        print(f"Graphiques générés en {end_time_plots - start_time_plots:.2f} secondes.")


    # 7. Génération du Rapport HTML (Séquentiel)
    # (Code inchangé...)
    if all_figures: # Ne générer le rapport que si des figures ont été créées
        generate_html_report(all_figures, args.output_html)
    else:
        print("Erreur: Aucune figure n'a été générée, le rapport HTML ne sera pas créé.")

    # --- Fin de l'exécution ---
    script_end_time = time.time()
    print("\n--- Analyse Terminée ---")
    print(f"Temps d'exécution total : {script_end_time - script_start_time:.2f} secondes.")
    if all_figures: # Afficher le message seulement si le rapport a été généré
        print(f"Rapport disponible : {args.output_html}")


if __name__ == "__main__":
    main()
