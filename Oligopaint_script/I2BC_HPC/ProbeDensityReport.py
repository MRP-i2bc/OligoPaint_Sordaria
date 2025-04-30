#!/usr/bin/env python3
import argparse
import os
import subprocess
import math
from collections import Counter, defaultdict
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import concurrent.futures
from functools import partial
import time
import shlex

# --- Fonction de Parsing FASTA (Similaire à CentroReport.py) ---
def parse_fasta(fasta_file):
    """Lit un fichier FASTA et retourne un dictionnaire {id: record}."""
    try:
        sequences = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
        print(f"Fichier FASTA '{os.path.basename(fasta_file)}' lu : {len(sequences)} séquences trouvées.")
        if not sequences:
            raise ValueError(f"Aucune séquence trouvée dans le fichier FASTA: {fasta_file}")
        return sequences
    except FileNotFoundError:
        print(f"Erreur CRITIQUE : Le fichier FASTA '{fasta_file}' n'a pas été trouvé.")
        exit(1)
    except Exception as e:
        print(f"Erreur CRITIQUE lors de la lecture du fichier FASTA '{fasta_file}': {e}")
        exit(1)

# --- NOUVELLE Fonction de Parsing GFF (Issue de CentroReport.py) ---
def parse_gff_for_features(gff_file, sequence_ids, target_gene_names=None):
    """
    Parse un fichier GFF/GFF3 pour extraire les positions de gènes spécifiques
    et identifier le locus ARNr.

    Args:
        gff_file (str): Chemin vers le fichier GFF/GFF3.
        sequence_ids (list): Liste des IDs de séquence attendus (du FASTA).
        target_gene_names (list, optional): Liste de noms/IDs de gènes à afficher spécifiquement.
                                            Si None, aucun gène spécifique n'est extrait.

    Returns:
        tuple: (gene_locations, rrna_locations)
            gene_locations (dict): {seq_id: [(start, end, name)], ...} pour les gènes cibles.
            rrna_locations (dict): {seq_id: (min_start, max_end), ...} pour le span ARNr.
    """
    print(f"\n--- Parsing du fichier GFF d'annotations : {gff_file} ---")
    gene_locations = {seq_id: [] for seq_id in sequence_ids}
    rrna_raw_locations = {seq_id: [] for seq_id in sequence_ids} # Pour calculer le span global

    if not gff_file or not os.path.exists(gff_file):
        print("Avertissement : Fichier GFF non fourni ou introuvable.")
        return gene_locations, {seq_id: None for seq_id in sequence_ids} # Retourne structures vides

    target_gene_set = set(target_gene_names) if target_gene_names else set()
    found_target_genes = set()

    try:
        with open(gff_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                seq_id = fields[0]
                feature_type = fields[2].lower() # Convertir en minuscule pour la comparaison
                attributes_str = fields[8]

                if seq_id not in sequence_ids: # Ignorer si le contig/chr n'est pas dans notre FASTA
                    continue

                try:
                    start_0based = int(fields[3]) - 1
                    end_0based_exclusive = int(fields[4])
                    if start_0based >= end_0based_exclusive: continue
                except ValueError:
                    continue # Ignorer si coordonnées invalides

                # --- Extraction des attributs (simplifiée) ---
                attributes = {}
                for part in attributes_str.split(';'):
                    if '=' in part:
                        key_value = part.strip().split('=', 1)
                        if len(key_value) == 2:
                           key, value = key_value
                           attributes[key.lower()] = value.strip().strip('"') # Clés en minuscules

                gene_name = attributes.get('name', attributes.get('id', attributes.get('gene_id', None)))

                # --- Trouver les Gènes Cibles ---
                if target_gene_set and feature_type in ['gene', 'pseudogene', 'mrna', 'cds']: # Chercher divers types de features de gènes
                    if gene_name and gene_name in target_gene_set:
                        # Éviter d'ajouter plusieurs fois le même gène (ex: si trouvé comme 'gene' puis 'mRNA')
                        # On pourrait stocker dans un set ou vérifier avant d'ajouter
                        is_already_added = any(g[2] == gene_name for g in gene_locations[seq_id])
                        if not is_already_added:
                             gene_locations[seq_id].append((start_0based, end_0based_exclusive, gene_name))
                             found_target_genes.add(gene_name)
                             # print(f"  DEBUG: Trouvé gène cible {gene_name} à {seq_id}:{start_0based}-{end_0based_exclusive}")

                # --- Trouver les composants ARNr ---
                is_rrna_component = False
                if feature_type in ['rrna', 'rrna_gene', 'rrna_operon']:
                    is_rrna_component = True
                elif gene_name and any(rrna_tag in gene_name.upper() for rrna_tag in ['5S', '18S', '28S', '5.8S', 'LSU', 'SSU', 'RIBOSOMAL_RNA']):
                    is_rrna_component = True
                elif attributes.get('product') and 'ribosomal rna' in attributes.get('product', '').lower():
                     is_rrna_component = True

                if is_rrna_component:
                    rrna_raw_locations[seq_id].append((start_0based, end_0based_exclusive))
                    # print(f"  DEBUG: Trouvé composant ARNr {gene_name or feature_type} à {seq_id}:{start_0based}-{end_0based_exclusive}")


        # --- Calculer le Span Global du Locus ARNr ---
        rrna_locations = {}
        for seq_id, positions in rrna_raw_locations.items():
            if positions:
                min_start = min(p[0] for p in positions)
                max_end = max(p[1] for p in positions)
                rrna_locations[seq_id] = (min_start, max_end)
                print(f"  Locus ARNr identifié sur {seq_id}: {min_start:,} - {max_end:,} bp")
            else:
                rrna_locations[seq_id] = None

        # --- Vérifier les gènes cibles non trouvés ---
        not_found_genes = target_gene_set - found_target_genes
        if target_gene_set and not_found_genes: # Afficher seulement si on cherchait des gènes spécifiques
            print(f"Avertissement : Les gènes cibles suivants n'ont pas été trouvés dans le fichier GFF : {', '.join(not_found_genes)}")

        print(f"Parsing GFF terminé.")
        if target_gene_set: print(f"  {len(found_target_genes)}/{len(target_gene_set)} gènes cibles trouvés.")
        return gene_locations, rrna_locations

    except Exception as e:
        print(f"Erreur majeure lors du parsing du fichier GFF {gff_file}: {e}")
        return ({seq_id: [] for seq_id in sequence_ids}, {seq_id: None for seq_id in sequence_ids})


# --- Fonctions pour l'Alignement (BLAST - inchangées) ---
def run_makeblastdb(reference_fasta, db_path_prefix, output_dir):
    # ... (code inchangé) ...
    """Crée une base de données BLAST à partir du génome de référence."""
    print("\n--- Création de la base de données BLAST ---")
    makeblastdb_executable = "makeblastdb"
    abs_ref_path = os.path.abspath(reference_fasta)
    abs_db_prefix = os.path.abspath(db_path_prefix) # Doit être juste le préfixe, pas le dir

    command = [
        makeblastdb_executable,
        "-in", abs_ref_path,
        "-dbtype", "nucl",
        "-out", abs_db_prefix, # Préfixe pour les fichiers de la DB
        "-parse_seqids" # Recommandé pour garder les IDs originaux
    ]
    print(f"Commande makeblastdb : {' '.join(shlex.quote(c) for c in command)}")

    try:
        # Pas besoin de spécifier cwd si abs_db_prefix inclut le chemin
        result = subprocess.run(command, capture_output=True, text=True, check=True, encoding='utf-8', errors='replace')
        print(result.stdout)
        if result.stderr:
            print("--- Sortie Erreur makeblastdb ---")
            print(result.stderr)
            print("-------------------------------")
        print("Base de données BLAST créée avec succès.")
        return True
    except FileNotFoundError:
        print(f"Erreur CRITIQUE : '{makeblastdb_executable}' introuvable. Vérifiez l'installation de BLAST+.")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Erreur CRITIQUE : makeblastdb a échoué (code {e.returncode}).")
        print(f"Commande : {' '.join(shlex.quote(c) for c in e.cmd)}")
        print(f"Sortie Erreur : {e.stderr}")
        return False
    except Exception as e:
        print(f"Erreur inattendue lors de l'exécution de makeblastdb : {e}")
        return False

def run_blastn(query_fasta, db_path_prefix, output_file, threads=1, blast_options=""):
    # ... (code inchangé) ...
    """Exécute blastn pour aligner les sondes sur la base de données du génome."""
    print(f"\n--- Exécution de blastn pour {os.path.basename(query_fasta)} ---")
    blastn_executable = "blastn"
    abs_query_path = os.path.abspath(query_fasta)
    abs_db_prefix = os.path.abspath(db_path_prefix)
    abs_output_file = os.path.abspath(output_file)

    # Format de sortie tabulaire standard
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    command = [
        blastn_executable,
        "-query", abs_query_path,
        "-db", abs_db_prefix,
        "-out", abs_output_file,
        "-outfmt", outfmt,
        "-num_threads", str(threads)
    ]
    # Ajouter les options blast supplémentaires si fournies
    if blast_options:
        command.extend(shlex.split(blast_options)) # Permet des options comme "-evalue 1e-5 -perc_identity 95"

    print(f"Commande blastn : {' '.join(shlex.quote(c) for c in command)}")

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True, encoding='utf-8', errors='replace')
        # blastn écrit généralement peu sur stdout/stderr si succès, mais on vérifie
        if result.stdout: print(f"Sortie standard blastn:\n{result.stdout}")
        if result.stderr: print(f"Sortie erreur blastn:\n{result.stderr}")

        if os.path.exists(abs_output_file):
             # Vérifier si le fichier n'est pas vide (succès mais aucun hit?)
             if os.path.getsize(abs_output_file) > 0:
                print(f"blastn terminé. Résultats dans : {abs_output_file}")
             else:
                print(f"Avertissement : blastn terminé, mais le fichier de sortie '{abs_output_file}' est vide (aucun alignement trouvé ?).")
             return True
        else:
             print(f"Erreur : blastn terminé (code 0) mais fichier de sortie '{abs_output_file}' non trouvé.")
             return False

    except FileNotFoundError:
        print(f"Erreur CRITIQUE : '{blastn_executable}' introuvable. Vérifiez l'installation de BLAST+.")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Erreur CRITIQUE : blastn a échoué (code {e.returncode}).")
        print(f"Commande : {' '.join(shlex.quote(c) for c in e.cmd)}")
        print(f"Sortie Erreur : {e.stderr}")
        return False
    except Exception as e:
        print(f"Erreur inattendue lors de l'exécution de blastn : {e}")
        return False

def parse_blast_output(blast_output_file, sequence_ids):
    # ... (code inchangé) ...
    """
    Parse le fichier de sortie BLAST tabulaire (outfmt 6).
    Retourne un dictionnaire {seq_id: [(start, end)]}.
    """
    print(f"--- Parsing du fichier BLAST : {os.path.basename(blast_output_file)} ---")
    # Initialise avec des listes vides pour tous les chromosomes/contigs attendus
    probe_hits = {seq_id: [] for seq_id in sequence_ids}
    hits_parsed = 0
    lines_read = 0

    if not os.path.exists(blast_output_file):
        print(f"Avertissement : Fichier BLAST '{blast_output_file}' introuvable.")
        return probe_hits

    try:
        with open(blast_output_file, 'r', encoding='utf-8') as f:
            for line in f:
                lines_read += 1
                if line.startswith('#'): continue # Ignorer commentaires éventuels
                fields = line.strip().split('\t')
                # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                if len(fields) < 12: continue

                sseqid = fields[1] # ID du chromosome/contig de référence

                if sseqid not in probe_hits:
                    # print(f"Avertissement: Hit sur '{sseqid}' non trouvé dans les séquences de référence. Ignoré.")
                    continue

                try:
                    sstart = int(fields[8])
                    send = int(fields[9])
                    # Convertir en coordonnées 0-based [start, end)
                    # Gérer le cas où sstart > send (alignement sur le brin complémentaire)
                    start_0based = min(sstart, send) - 1
                    end_0based_exclusive = max(sstart, send)

                    if start_0based < end_0based_exclusive:
                        probe_hits[sseqid].append((start_0based, end_0based_exclusive))
                        hits_parsed += 1
                    else:
                         print(f"Avertissement : Coordonnées BLAST invalides ignorées (start >= end): {line.strip()}")

                except ValueError:
                    print(f"Avertissement : Ligne BLAST mal formée (coordonnées non numériques?) ignorée : {line.strip()}")
                    continue

        print(f"Parsing BLAST terminé. {hits_parsed} alignements trouvés sur {lines_read} lignes lues.")
        return probe_hits

    except Exception as e:
        print(f"Erreur majeure lors du parsing du fichier BLAST {blast_output_file}: {e}")
        return {seq_id: [] for seq_id in sequence_ids} # Retourner structure vide

# --- Fonctions de Calcul de Densité et Worker (inchangées) ---
def calculate_density_probes(seq_len, window_size, probe_hits_positions):
    # ... (code inchangé) ...
    """
    Calcule la densité de 'hits' de sondes par fenêtre.
    Densité = Nombre de hits dont le début est dans la fenêtre / longueur fenêtre.
    Retourne: (window_centers, densities)
    """
    if seq_len == 0: return [], []
    window_starts = list(range(0, seq_len, window_size))
    num_windows = len(window_starts)
    densities = [0.0] * num_windows
    window_centers = [s + min(window_size / 2, (seq_len - s) / 2) for s in window_starts]

    if not probe_hits_positions: # Si aucune sonde n'a mappé sur ce chromosome
        return window_centers, densities

    # Extraire seulement les positions de début pour le comptage par fenêtre
    hit_starts = [start for start, end in probe_hits_positions]

    # Compter combien de hits débutent dans chaque fenêtre
    hit_indices = [pos // window_size for pos in hit_starts if pos // window_size < num_windows]
    counts_per_window = Counter(hit_indices)

    for i, win_start in enumerate(window_starts):
        win_end = min(win_start + window_size, seq_len)
        win_len = win_end - win_start
        if win_len <= 0: continue
        # Densité = Nombre de hits / longueur fenêtre (hits par bp)
        # Alternative: Densité = Nombre de hits bruts (counts_per_window.get(i, 0))
        densities[i] = counts_per_window.get(i, 0) *1000.0/ win_len

    return window_centers, densities


def calculate_densities_for_chromosome_worker(args_tuple):
    # ... (code inchangé) ...
    """
    Fonction worker pour calculer les densités de chaque set de sondes pour un chromosome.
    """
    seq_id, seq_len, window_size, all_probe_hits_for_seq = args_tuple
    # all_probe_hits_for_seq est un dict: {probe_set_name: [(start, end), ...]}

    densities_dict = {}
    centers = [] # Initialiser au cas où il n'y a pas de hits

    # Calculer la densité pour chaque set de sondes fourni
    for probe_set_name, probe_hits_list in all_probe_hits_for_seq.items():
        # Utiliser la fonction de densité spécifique aux probes (comptage de hits)
        centers, density_values = calculate_density_probes(seq_len, window_size, probe_hits_list)
        densities_dict[probe_set_name] = (centers, density_values) # Stocker avec le nom du set

    # S'assurer qu'on a des 'centers' même si un set de sondes n'a pas de hits
    if not centers and seq_len > 0:
         window_starts = list(range(0, seq_len, window_size))
         centers = [s + min(window_size / 2, (seq_len - s) / 2) for s in window_starts]

    # Retourner les données nécessaires pour le graphique
    return seq_id, seq_len, centers, densities_dict

# --- Fonctions Plotting (MODIFIÉE) et Rapport HTML ---

# MODIFICATION de create_plots pour accepter les annotations
# Dans ProbeDensityReport.py

def create_plots(seq_id, seq_len, window_centers, densities_dict, window_size,
                 gene_positions=None, rrna_locus=None):
    """Crée un graphique Plotly pour un chromosome avec les densités de sondes
       (uniquement celles présentes sur ce chr) et ajoute des annotations."""

    fig = go.Figure()
    if not window_centers:
        print(f"Avertissement : Chromosome {seq_id} trop court (longueur {seq_len}) pour la taille de fenêtre {window_size}, graphique vide généré.")
        fig.update_layout(title=f"Densité des Sondes - Chromosome {seq_id} (Trop court)")
        return fig

    plot_added = False
    probe_sets_plotted = [] # Garder une trace des sondes ajoutées

    # --- MODIFICATION : Filtrer les densités avant de plotter ---
    print(f"  Génération graphique pour {seq_id}. Analyse des densités pour les sets : {', '.join(densities_dict.keys())}") # Debug
    for name, data in densities_dict.items(): # 'name' est le probe_set_name
        if data:
            centers, values = data
            # Vérifier si la liste des valeurs de densité contient des valeurs non nulles
            # Utiliser sum() > epsilon est plus robuste pour les flottants que any(v > 0)
            if centers and values and sum(values) > 1e-9: # Seuil très petit pour considérer comme non nul
                # Si la densité est non nulle, ajouter la trace
                fig.add_trace(go.Scatter(x=centers, y=values, mode='lines', name=name))
                plot_added = True
                probe_sets_plotted.append(name) # Ajouter à la liste des tracés
            # else:
                # Optionnel : log si une sonde est ignorée sur ce chr
                # print(f"    -> Ignoré {name} sur {seq_id} (densité nulle).")

    if probe_sets_plotted:
         print(f"    -> Affichage des densités pour : {', '.join(probe_sets_plotted)}")
    # --- FIN MODIFICATION ---

    if not plot_added:
         print(f"Avertissement : Aucune donnée de densité de sonde non nulle à tracer pour {seq_id}.")
         fig.add_annotation(x=seq_len/2, y=0.5, text="Aucune densité de sonde significative trouvée", showarrow=False, yref="paper") # Message mis à jour

    # --- Ajout du Locus ARNr (inchangé) ---
    if rrna_locus:
        # ... (code identique) ...
        rrna_start, rrna_end = rrna_locus
        fig.add_vrect(
            x0=rrna_start, x1=rrna_end,
            fillcolor="LightSkyBlue", opacity=0.2, layer="below", line_width=0, # Couleur différente ?
            annotation_text="ARNr Locus", annotation_position="top left",
            annotation=dict(font=dict(size=10, color="DarkBlue"))
        )

    # --- Ajout des Positions des Gènes (inchangé) ---
    if gene_positions:
        # ... (code identique) ...
        y_annotation = 1.05
        for gene_start, gene_end, gene_name in gene_positions:
            gene_midpoint = (gene_start + gene_end) / 2
            fig.add_vline(
                x=gene_midpoint,
                line_width=1, line_dash="dash", line_color="ForestGreen",
            )
            fig.add_annotation(
                x=gene_midpoint, y=y_annotation, yref="paper", text=gene_name,
                showarrow=False, font=dict(size=9, color="ForestGreen"), textangle=-45
            )

    # --- Mise en page (inchangée) ---
    fig.update_layout(
        title=f"Densité des Sondes - Chromosome {seq_id} (Longueur: {seq_len:,} bp)",
        xaxis_title=f"Position sur le chromosome (Fenêtres de {window_size:,} bp)",
        yaxis_title="Densité de Hits de Sonde (hits/kb)",
        xaxis=dict(range=[0, seq_len]),
        yaxis=dict(range=[0, None], autorange=True),
        hovermode='x unified',
        legend_title_text='Sondes Actives', # On peut ajuster le titre de la légende
        margin=dict(t=80)
    )
    return fig


def generate_html_report(figures_dict, output_file):
    """Génère le rapport HTML avec les graphiques Plotly."""
    # (Identique à la version précédente de CentroReport.py)
    print(f"\n--- Génération du rapport HTML : {output_file} ---")
    start_time = time.time()
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("<!DOCTYPE html>\n<html>\n<head>\n")
            f.write("<meta charset='utf-8' />\n")
            # Mettre à jour le titre de la page HTML
            f.write("<title>Rapport Densité des Sondes</title>\n")
            f.write("<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n")
            f.write("<style> body { font-family: sans-serif; margin: 20px; } h1, h2 { color: #333; } hr { border: 0; height: 1px; background: #ccc; } .plot-container { margin-bottom: 30px; border: 1px solid #eee; padding: 10px; box-shadow: 2px 2px 5px #eee; } </style>\n")
            f.write("</head>\n<body>\n")
            # Mettre à jour le titre principal
            f.write("<h1>Analyse de Densité des Sondes sur le Génome</h1>\n")

            try:
                sorted_ids = sorted(figures_dict.keys(), key=lambda x: int(''.join(filter(str.isdigit, x))) if any(char.isdigit() for char in x) else float('inf'))
            except:
                sorted_ids = sorted(figures_dict.keys())

            for seq_id in sorted_ids:
                fig = figures_dict[seq_id]
                f.write(f"<div class='plot-container'>\n")
                f.write(f"<h2>Chromosome / Contig: {seq_id}</h2>\n")
                plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
                f.write(plot_html)
                f.write("</div>\n")

            f.write("</body>\n</html>\n")
        end_time = time.time()
        print(f"Rapport généré en {end_time - start_time:.2f} secondes.")
    except Exception as e:
        print(f"Erreur lors de la génération du rapport HTML : {e}")


# --- Fonction Principale (MODIFIÉE) ---

def main():
    parser = argparse.ArgumentParser(
        description="Calcule et visualise la densité de sondes sur un génome, avec annotations optionnelles.", # Description mise à jour
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # --- Arguments Alignement/Densité ---
    parser.add_argument("reference_fasta", help="Fichier FASTA du génome de référence.")
    parser.add_argument("probe_files", nargs='+', help="Un ou plusieurs fichiers FASTA contenant les séquences de sondes.")
    parser.add_argument("-o", "--output_html", default="probe_density_report.html", help="Nom du fichier HTML de sortie.")
    parser.add_argument("--output_dir", default="probe_analysis_output", help="Répertoire pour les fichiers intermédiaires (DB BLAST, résultats BLAST).")
    parser.add_argument("-w", "--window_size", type=int, default=50000, help="Taille des fenêtres d'analyse (en bp).")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Nombre de threads/processus (pour BLAST et calculs parallèles). 0 = tous les cœurs.")
    parser.add_argument("--blast_options", default="", help="Options supplémentaires pour blastn (ex: '-evalue 1e-10 -perc_identity 95').")

    # --- NOUVEAUX Arguments pour Annotations GFF ---
    parser.add_argument("--gff-annotations", metavar='FILE', help="Fichier GFF/GFF3 contenant les annotations (gènes, ARNr) à afficher sur les graphiques.")
    parser.add_argument("--genes-list", metavar='GENE,GENE,...', help="Liste de noms/IDs de gènes (séparés par virgule) à afficher, tels qu'ils apparaissent dans les attributs du GFF.")


    args = parser.parse_args()

    # --- Début de l'exécution ---
    script_start_time = time.time()
    print("--- Initialisation de l'Analyse de Densité de Sondes ---")
    print(f"Génome Référence: {args.reference_fasta}")
    print(f"Fichiers Sondes: {', '.join(args.probe_files)}")
    if args.gff_annotations:
        print(f"Fichier Annotations: {args.gff_annotations}")
        if args.genes_list:
            print(f"Gènes à afficher: {args.genes_list}")
    print(f"Sortie HTML: {args.output_html}")
    print(f"Répertoire de travail: {args.output_dir}")
    print(f"Taille Fenêtre: {args.window_size:,} bp")

    num_workers = max(1, os.cpu_count() if args.threads <= 0 else min(args.threads, os.cpu_count()))
    print(f"Utilisation de {num_workers} workers.")

    # 1. Créer le répertoire de sortie
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except OSError as e:
        print(f"Erreur CRITIQUE: Impossible de créer le répertoire de sortie '{args.output_dir}': {e}")
        exit(1)

    # 2. Lire le génome de référence
    ref_sequences = parse_fasta(args.reference_fasta)
    ref_sequence_ids = list(ref_sequences.keys())

    # --- NOUVEAU: Parser les annotations GFF SI fourni ---
    gene_locations = {seq_id: [] for seq_id in ref_sequence_ids}
    rrna_locations = {seq_id: None for seq_id in ref_sequence_ids}
    if args.gff_annotations:
        gene_target_list = args.genes_list.split(',') if args.genes_list else None
        gene_locations, rrna_locations = parse_gff_for_features(args.gff_annotations, ref_sequence_ids, gene_target_list)
    else:
        print("Info: Aucune annotation GFF fournie (--gff-annotations), les graphiques n'afficheront pas les gènes ou l'ARNr.")
    # --- Fin Parsing Annotations ---

    # 3. Créer la base de données BLAST
    db_prefix = os.path.join(args.output_dir, os.path.basename(args.reference_fasta) + "_blastdb")
    if not run_makeblastdb(args.reference_fasta, db_prefix, args.output_dir):
        print("Échec de la création de la base de données BLAST. Arrêt.")
        exit(1)

    # 4. Aligner chaque fichier de sondes et parser les résultats
    all_probe_hits = defaultdict(lambda: {seq_id: [] for seq_id in ref_sequence_ids})
    blast_success_count = 0
    for probe_fasta_file in args.probe_files:
        probe_set_name = os.path.splitext(os.path.basename(probe_fasta_file))[0]
        blast_output_file = os.path.join(args.output_dir, f"{probe_set_name}_vs_{os.path.basename(args.reference_fasta)}.blast.tsv")

        if run_blastn(probe_fasta_file, db_prefix, blast_output_file, threads=num_workers, blast_options=args.blast_options):
            parsed_hits = parse_blast_output(blast_output_file, ref_sequence_ids)
            # Utiliser defaultdict simplifie l'ajout
            for seq_id, hits in parsed_hits.items():
                 all_probe_hits[probe_set_name][seq_id].extend(hits)
            blast_success_count += 1
        else:
            print(f"Erreur lors de l'exécution de blastn pour {probe_fasta_file}. Ce set de sondes sera ignoré.")

    if blast_success_count == 0 :
         print("Erreur CRITIQUE : Aucun alignement de sonde n'a pu être obtenu (tous les BLAST ont échoué). Arrêt.")
         exit(1)
    elif blast_success_count < len(args.probe_files):
         print("Avertissement: Certains alignements de sondes ont échoué.")

    # 5. Préparation des données pour le calcul parallèle des densités
    print(f"\n--- Préparation du calcul parallèle des densités ({num_workers} workers) ---")
    tasks_args = []
    for seq_id, record in ref_sequences.items():
        seq_len = len(record.seq)
        probe_hits_for_this_seq = {
            probe_set_name: hits_dict.get(seq_id, [])
            for probe_set_name, hits_dict in all_probe_hits.items()
        }
        # Le tuple d'arguments doit correspondre à la signature de calculate_densities_for_chromosome_worker
        tasks_args.append((
            seq_id, seq_len, args.window_size, probe_hits_for_this_seq
        ))

    # 6. Exécution du calcul parallèle des densités
    print(f"\n--- Lancement du calcul parallèle des densités ---")
    start_time_density = time.time()
    all_plot_data = {}
    processed_count = 0
    total_seqs = len(ref_sequences)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_seqid = {executor.submit(calculate_densities_for_chromosome_worker, args_tuple): args_tuple[0] for args_tuple in tasks_args}
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


    # 7. Génération des Graphiques (MODIFIÉ pour passer les annotations GFF)
    print("\n--- Génération des graphiques ---")
    all_figures = {}
    start_time_plots = time.time()
    if not all_plot_data:
         print("Erreur CRITIQUE: Aucune donnée de densité n'a été calculée.")
    else:
        for seq_id, plot_data_tuple in all_plot_data.items():
            seq_len, centers, densities_dict = plot_data_tuple
            # Récupérer les annotations GFF pour ce chromosome
            genes_for_seq = gene_locations.get(seq_id, [])
            rrna_for_seq = rrna_locations.get(seq_id, None)

            # Passer les annotations à create_plots
            all_figures[seq_id] = create_plots(seq_id, seq_len, centers, densities_dict, args.window_size,
                                               gene_positions=genes_for_seq,
                                               rrna_locus=rrna_for_seq)
        end_time_plots = time.time()
        print(f"Graphiques générés en {end_time_plots - start_time_plots:.2f} secondes.")

    # 8. Génération du Rapport HTML
    if all_figures:
        generate_html_report(all_figures, args.output_html)
    else:
        print("Erreur: Aucune figure n'a été générée, le rapport HTML ne sera pas créé.")

    # --- Fin de l'exécution ---
    script_end_time = time.time()
    print("\n--- Analyse Terminée ---")
    print(f"Temps d'exécution total : {script_end_time - script_start_time:.2f} secondes.")
    if all_figures:
        print(f"Rapport disponible : {args.output_html}")

if __name__ == "__main__":
    main()