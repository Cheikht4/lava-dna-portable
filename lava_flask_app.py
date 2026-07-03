#!/usr/bin/env python3
"""
LAVA-DNA Interface Flask
Interface web stable pour les scripts LAVA STEM et LOOP
"""

from flask import Flask, render_template, request, redirect, url_for, flash, jsonify, send_file, session, g, make_response
from datetime import datetime
import os
import subprocess
import tempfile
import time
import uuid
import threading
from jinja2 import pass_context
from werkzeug.utils import secure_filename
import json
import re


# Dictionnaire de traductions
TRANSLATIONS = {
    'fr': {
        'title': 'Interface LAVA-DNA',
        'subtitle': 'Design de primers LAMP avec LAVA - Interface web stable',
        'upload_title': '1. Upload du fichier FASTA',
        'upload_label': 'Fichier de séquences FASTA :',
        'upload_formats': 'Formats acceptés : .fas, .fasta, .fa, .txt (1GB max)',
        'upload_current': 'Fichier actuel :',
        'upload_button': 'Uploader',
        'config_title': '2. Configuration des paramètres LAVA',
        'primer_type': 'Type de primers',
        'stem_primers': 'STEM primers',
        'loop_primers': 'LOOP primers',
        'lamp_mode': 'Mode LAMP',
        'lamp_classic': 'LAMP Classique (6 primers)',
        'lamp_classic_desc': 'F3, F2, F1, B1, B2, B3 seulement',
        'lamp_enriched': 'LAMP Enrichi (8 primers)',
        'lamp_enriched_stem': 'F3, F2, F1, FSTEM, BSTEM, B1, B2, B3',
        'lamp_enriched_loop': 'F3, F2, F1, FLOOP, BLOOP, B1, B2, B3',
        'general_params': 'Paramètres généraux',
        'max_signature_length': 'Longueur max signature',
        'max_primers_generated': 'Max primers générés',
        'match_minimum': '% Match minimum',
        'match_minimum_desc': 'Correspondance stricte minimum',
        'match_after_iupac': '% Match après IUPAC',
        'match_after_iupac_desc': 'Match après bases dégénérées',
        'primer_elimination': '% Élimination primer',
        'primer_elimination_desc': 'Seuil d\'élimination final',
        'min_noise_freq': 'Fréquence min. bruit',
        'min_noise_freq_desc': 'Ignorer bases < X% (défaut 5%)',
        'outer_primers': 'Outer Primers (F3/B3)',
        'middle_primers': 'Middle Primers (F2/B2)',
        'inner_primers': 'Inner Primers (F1/B1)',
        'target_length': 'Longueur cible',
        'min_length': 'Longueur min',
        'max_length': 'Longueur max',
        'target_tm': 'Tm cible (°C)',
        'min_tm': 'Tm min (°C)',
        'max_tm': 'Tm max (°C)',
        'loop_advanced_params': 'Paramètres LOOP avancés',
        'loop_target_length': 'Longueur cible LOOP',
        'loop_min_length': 'Longueur min LOOP',
        'loop_max_length': 'Longueur max LOOP',
        'loop_target_tm': 'Tm cible LOOP (°C)',
        'loop_min_tm': 'Tm min LOOP (°C)',
        'loop_max_tm': 'Tm max LOOP (°C)',
        'loop_min_gap': 'Gap min LOOP',
        'loop_min_gap_desc': 'Minimum 15nt entre F2 et F1',
        'stem_advanced_params': 'Paramètres STEM avancés',
        'stem_target_length': 'Longueur cible STEM',
        'stem_min_length': 'Longueur min STEM',
        'stem_max_length': 'Longueur max STEM',
        'stem_target_tm': 'Tm cible STEM (°C)',
        'stem_min_tm': 'Tm min STEM (°C)',
        'stem_max_tm': 'Tm max STEM (°C)',
        'saved_successfully': 'Paramètres sauvegardés avec succès',
        'advanced_params': 'Paramètres Avancés',
        'thermo_conditions': 'Thermodynamique & Conditions',
        'architecture_geometry': 'Architecture & Géométrie',
        'diversity_mismatches': 'Diversité & Mismatches',
        'execution_config': 'Configuration Exécution',
        'salt_mono': 'Sel Monovalent (mM)',
        'salt_div': 'Sel Divalent (mM)',
        'dntp_conc': 'dNTP (mM)',
        'dna_conc': 'ADN (nM)',
        'max_tm_diff': 'Différence Tm Max (°C)',
        'min_signatures': 'Couverture Cible Minimum (%)',
        'max_overlap': 'Chevauchement Max (%)',
        'entropy_threshold': 'Seuil Entropie',
        'max_poly_bases': 'Max Poly-Bases',
        'max_total_degenerate_bases': 'Max Bases Dég. (Total)',
        'max_consecutive_degenerate_bases': 'Max Bases Dég. (Consécutive)',
        'max_3prime_degenerate_bases': 'Max Bases Dég. (3\')',
        'max_tolerated_mismatches': 'Tolérance Mismatches (Biol)',
        'three_prime_zone': 'Zone 3\' (bases)',
        'max_dist_outer_middle': 'Dist. Max Outer-Middle',
        'max_dist_middle_inner': 'Dist. Max Middle-Inner',
        'penalty_plateau': 'Plateau de Pénalité (0.1-0.5)',
        'penalty_slope': 'Pente Sigmoïde (0.05-0.5)',
        # Clés manquantes / Missing keys
        'resolve_overlap_label': 'Priorité de Dédoublonnage (Chevauchement)',
        'resolve_overlap_penalty': 'Pénalité Biochimique et Géométrique (Défaut/Sûr)',
        'resolve_overlap_coverage': 'Pourcentage de Couverture (Universalité)',
        'resolve_overlap_desc': 'Critère utilisé pour garder le "champion" d\'une région ciblée par plusieurs signatures.',
        'penalty_plateau_desc': 'Ratio zone "confort" (ex: 0.25)',
        'penalty_slope_desc': 'Pente sigmoïde (ex: 0.15)',
        'min_signatures_desc': '% minimum de séquences cibles que la signature doit amplifier. Ex: 1 = tolérant, 70 = strict.',
        'spatial_reduction_title': 'Réduction spatiale des candidats',
        'window_size_label': 'Taille de fenêtre (nt)',
        'window_size_desc': 'Largeur de fenêtre génomique en nt. 0 = désactivé. Ex: 10 = fenêtre de 10nt.',
        'max_per_window_label': 'Max candidats / fenêtre',
        'max_per_window_desc': 'Nombre max d\'amorces gardées par fenêtre. 0 = désactivé. Ex: 3 = 3 meilleurs par fenêtre.',
        'spatial_reduction_info': 'Garde les K meilleurs candidats par fenêtre de W nucléotides. Réduit le temps de calcul sans sacrifier la diversité spatiale. Valeurs recommandées : Fenêtre=10, Max=3.',
        'save_params': 'Sauvegarder les paramètres',
        'execution_title': '3. Exécution LAVA',
        'output_name': 'Nom du fichier de sortie',
        'execute_button': 'Lancer l\'exécution',
        'language': 'Langue',
        'params_updated': 'Paramètres mis à jour',
        'file_uploaded': 'Fichier uploadé avec succès',
        'no_file_uploaded': 'Aucun fichier FASTA uploadé',
        'invalid_file': 'Type de fichier non autorisé',
        'execution_list': 'Liste des exécutions',
        'new_execution': 'Nouvelle exécution',
        'refresh': 'Actualiser',
        'no_execution': 'Aucune exécution',
        'launch_first': 'Lancez votre première exécution LAVA depuis l\'accueil',
        'status': 'Statut',
        'help': 'Aide',
        'stem_vs_loop': 'STEM vs LOOP',
        'stem_desc': 'Primers avec structure tige-boucle',
        'loop_desc': 'Primers en boucle simple',
        'accepted_formats': 'Formats acceptés',
        'file_size_max': 'Taille max : 1GB',
        'parameters': 'Paramètres',
        'adjust_needs': 'Ajustez selon vos besoins',
        'default_values': 'Les valeurs par défaut conviennent généralement',
        'running_executions': 'Exécutions en cours',
        'fasta_file': 'Fichier FASTA',
        'name_no_extension': 'Nom sans extension',
        'files_will_be_named': 'Les fichiers seront nommés : nom.primers, nom.all_signatures, etc.',
        'upload_first_message': 'Veuillez d\'abord uploader un fichier FASTA pour pouvoir lancer l\'exécution.',
        'no_file': 'Aucun fichier',
        'home': 'Accueil',
        'executions': 'Exécutions',
        'monitoring': 'Monitoring de l\'exécution',
        'script': 'Script',
        'file': 'Fichier',
        'output': 'Sortie',
        'status': 'Statut',
        'back_to_executions': 'Retour aux exécutions',
        'execution_completed': 'Exécution terminée',
        'detailed_info': 'Informations détaillées',
        'technical_details': 'Détails techniques',
        'hide_technical_details': 'Masquer détails techniques',
        'show_technical_details': 'Voir détails techniques',
        'real_time_logs': 'Logs en temps réel',
        'last_100_lines': '(dernières 100 lignes)',
        'all_logs_display': 'Tous les logs disponibles',
        'view_all_logs': 'Voir tous les logs',
        'normal_display': 'Affichage normal',
        'download_logs': 'Télécharger logs',
        'results_available': 'Résultats disponibles',
        'download_all': 'Tout télécharger',
        'select_all': 'Sélectionner tout',
        'download_selected': 'Télécharger sélection',
        'error_detected': 'Erreur détectée',
        'need_help': 'Besoin d\'aide ?',
        'refreshing': 'Actualisation...',
        'status_running': 'En cours',
        'status_completed': 'Terminé',
        'status_completed_no_results': 'Terminé (0 signature)',
        'status_completed_unknown': 'Terminé (?)',
        'status_error': 'Erreur',
        'status_stopped': 'Arrêté',
        'status_starting': 'Démarrage',
        'error_input_not_aligned': "Le fichier fourni ne semble pas être un alignement multiple : les séquences ont des longueurs différentes (ou une seule séquence fournie). Veuillez aligner vos séquences (par exemple avec MAFFT ou Clustal) avant de les soumettre.",
        'msg_signatures_found': "✅ {count} signature(s) trouvée(s)",
        'msg_completed_no_results': "Exécution terminée : aucune signature trouvée avec les paramètres actuels. Essayez d'assouplir les seuils (couverture, dégénérescence) ou de vérifier l'alignement d'entrée.",
        'msg_completed_unknown': "ℹ️ Exécution terminée : le statut des signatures n'a pas pu être déterminé dans les logs.",
        'msg_exec_error': "❌ Erreur d'exécution (code {code})",
        'msg_error_detail': "\nDétail : {detail}",
        'sugg_memory': "\n💡 Suggestion: Essayez avec un fichier plus petit ou des paramètres moins stricts",
        'sugg_file': "\n💡 Suggestion: Vérifiez que le fichier FASTA existe et est accessible",
        'sugg_perms': "\n💡 Suggestion: Problème de permissions, contactez l'administrateur",
        'full_id': "ID complet",
        'created': "Créé",
        'started': "Démarré",
        'ended': "Terminé",
        'duration': "Durée",
        'log_lines': "Lignes de logs",
        'command': "Commande",
        'stop_confirm': "Arrêter l'exécution ?",
        'help_tip_1': "Vérifiez que votre fichier FASTA est correctement formaté",
        'help_tip_2': "Essayez avec des paramètres plus permissifs",
        'help_tip_3': "Contactez le support si le problème persiste",
        'progress_init': "Initialisation...",
        'progress_done': "Terminé",
        'total_lines_gen': "lignes au total / {total} générées",
        'lines_count': "lignes",
        'view_btn': "Voir",
        'stop_btn': "Stop",
        'results_btn': "Résultats",
        'footer_text': 'Interface LAVA-DNA Flask - Interface stable pour le design d\'amorces LAMP'
    },
    'en': {
        'title': 'LAVA-DNA Interface',
        'subtitle': 'LAMP primer design with LAVA - Stable web interface',
        'upload_title': '1. FASTA File Upload',
        'upload_label': 'FASTA sequence file:',
        'upload_formats': 'Accepted formats: .fas, .fasta, .fa, .txt (1GB max)',
        'upload_current': 'Current file:',
        'upload_button': 'Upload',
        'config_title': '2. LAVA Parameters Configuration',
        'primer_type': 'Primer type',
        'stem_primers': 'STEM primers',
        'loop_primers': 'LOOP primers',
        'lamp_mode': 'LAMP Mode',
        'lamp_classic': 'Classic LAMP (6 primers)',
        'lamp_classic_desc': 'F3, F2, F1, B1, B2, B3 only',
        'lamp_enriched': 'Enriched LAMP (8 primers)',
        'lamp_enriched_stem': 'F3, F2, F1, FSTEM, BSTEM, B1, B2, B3',
        'lamp_enriched_loop': 'F3, F2, F1, FLOOP, BLOOP, B1, B2, B3',
        'general_params': 'General parameters',
        'max_signature_length': 'Max signature length',
        'max_primers_generated': 'Max primers generated',
        'match_minimum': '% Minimum match',
        'match_minimum_desc': 'Minimum strict match',
        'match_after_iupac': '% Match after IUPAC',
        'match_after_iupac_desc': 'Match after degenerate bases',
        'primer_elimination': '% Primer elimination',
        'primer_elimination_desc': 'Final elimination threshold',
        'min_noise_freq': 'Min noise frequency',
        'min_noise_freq_desc': 'Ignore bases < X% (default 5%)',
        'outer_primers': 'Outer Primers (F3/B3)',
        'middle_primers': 'Middle Primers (F2/B2)',
        'inner_primers': 'Inner Primers (F1/B1)',
        'target_length': 'Target length',
        'min_length': 'Min length',
        'max_length': 'Max length',
        'target_tm': 'Target Tm (°C)',
        'min_tm': 'Min Tm (°C)',
        'max_tm': 'Max Tm (°C)',
        'loop_advanced_params': 'Advanced LOOP parameters',
        'loop_target_length': 'LOOP target length',
        'loop_min_length': 'LOOP min length',
        'loop_max_length': 'LOOP max length',
        'loop_target_tm': 'LOOP target Tm (°C)',
        'loop_min_tm': 'LOOP min Tm (°C)',
        'loop_max_tm': 'LOOP max Tm (°C)',
        'loop_min_gap': 'LOOP min gap',
        'loop_min_gap_desc': 'Minimum 15nt between F2 and F1',
        'stem_advanced_params': 'Advanced STEM parameters',
        'stem_target_length': 'STEM target length',
        'stem_min_length': 'STEM min length',
        'stem_max_length': 'STEM max length',
        'stem_target_tm': 'STEM target Tm (°C)',
        'stem_min_tm': 'STEM min Tm (°C)',
        'stem_max_tm': 'STEM max Tm (°C)',
        'saved_successfully': 'Parameters saved successfully',
        'advanced_params': 'Advanced Parameters',
        'thermo_conditions': 'Thermodynamics & Conditions',
        'architecture_geometry': 'Architecture & Geometry',
        'diversity_mismatches': 'Diversity & Mismatches',
        'execution_config': 'Execution Configuration',
        'salt_mono': 'Monovalent Salt (mM)',
        'salt_div': 'Divalent Salt (mM)',
        'dntp_conc': 'dNTP (mM)',
        'dna_conc': 'DNA (nM)',
        'max_tm_diff': 'Max Tm Diff (°C)',
        'min_signatures': 'Min Target Coverage (%)',
        'max_overlap': 'Max Overlap (%)',
        'entropy_threshold': 'Entropy Threshold',
        'max_poly_bases': 'Max Poly-Bases',
        'max_total_degenerate_bases': 'Max Degen. Bases (Total)',
        'max_consecutive_degenerate_bases': 'Max Degen. Bases (Consecutive)',
        'max_3prime_degenerate_bases': 'Max Degen. Bases (3\')',
        'max_tolerated_mismatches': 'Biology Mismatch Tolerance',
        'three_prime_zone': '3\' Zone (bases)',
        'max_dist_outer_middle': 'Max Dist. Outer-Middle',
        'max_dist_middle_inner': 'Max Dist. Middle-Inner',
        'penalty_plateau': 'Penalty Plateau (0.1-0.5)',
        'penalty_slope': 'Sigmoid Slope (0.05-0.5)',
        # Clés manquantes EN / Missing EN keys
        'resolve_overlap_label': 'Overlap Deduplication Priority',
        'resolve_overlap_penalty': 'Biochemical & Geometric Penalty (Default/Safe)',
        'resolve_overlap_coverage': 'Coverage Percentage (Universality)',
        'resolve_overlap_desc': 'Criterion used to keep the "champion" of a region targeted by multiple signatures.',
        'penalty_plateau_desc': '"Comfort zone" ratio (e.g. 0.25)',
        'penalty_slope_desc': 'Sigmoid slope (e.g. 0.15)',
        'min_signatures_desc': 'Minimum % of target sequences the signature must amplify. E.g. 1 = tolerant, 70 = strict.',
        'spatial_reduction_title': 'Spatial Candidate Reduction',
        'window_size_label': 'Window size (nt)',
        'window_size_desc': 'Genomic window width in nt. 0 = disabled. E.g. 10 = 10nt window.',
        'max_per_window_label': 'Max candidates / window',
        'max_per_window_desc': 'Max primers kept per window. 0 = disabled. E.g. 3 = 3 best per window.',
        'spatial_reduction_info': 'Keeps the K best candidates per W-nucleotide window. Drastically reduces computation time without sacrificing spatial diversity. Recommended: Window=10, Max=3.',
        'save_params': 'Save parameters',
        'execution_title': '3. LAVA Execution',
        'output_name': 'Output file name',
        'execute_button': 'Launch execution',
        'language': 'Language',
        'params_updated': 'Parameters updated',
        'file_uploaded': 'File uploaded successfully',
        'no_file_uploaded': 'No FASTA file uploaded',
        'invalid_file': 'File type not allowed',
        'execution_list': 'Execution list',
        'new_execution': 'New execution',
        'refresh': 'Refresh',
        'no_execution': 'No execution',
        'launch_first': 'Launch your first LAVA execution from home',
        'status': 'Status',
        'help': 'Help',
        'stem_vs_loop': 'STEM vs LOOP',
        'stem_desc': 'Primers with stem-loop structure',
        'loop_desc': 'Simple loop primers',
        'accepted_formats': 'Accepted formats',
        'file_size_max': 'Max size: 1GB',
        'parameters': 'Parameters',
        'adjust_needs': 'Adjust according to your needs',
        'default_values': 'Default values generally work well',
        'running_executions': 'Running executions',
        'fasta_file': 'FASTA file',
        'name_no_extension': 'Name without extension',
        'files_will_be_named': 'Files will be named: name.primers, name.all_signatures, etc.',
        'upload_first_message': 'Please first upload a FASTA file to launch execution.',
        'no_file': 'No file',
        'home': 'Home',
        'executions': 'Executions',
        'monitoring': 'Execution monitoring',
        'script': 'Script',
        'file': 'File',
        'output': 'Output',
        'status': 'Status',
        'back_to_executions': 'Back to executions',
        'execution_completed': 'Execution completed',
        'detailed_info': 'Detailed information',
        'technical_details': 'Technical details',
        'hide_technical_details': 'Hide technical details',
        'show_technical_details': 'Show technical details',
        'real_time_logs': 'Real-time logs',
        'last_100_lines': '(last 100 lines)',
        'all_logs_display': 'All available logs',
        'view_all_logs': 'View all logs',
        'normal_display': 'Normal display',
        'download_logs': 'Download logs',
        'results_available': 'Results available',
        'download_all': 'Download all',
        'select_all': 'Select all',
        'download_selected': 'Download selected',
        'error_detected': 'Error detected',
        'need_help': 'Need help?',
        'refreshing': 'Refreshing...',
        'status_running': 'Running',
        'status_completed': 'Completed',
        'status_completed_no_results': 'Completed (0 signatures)',
        'status_completed_unknown': 'Completed (?)',
        'status_error': 'Error',
        'status_stopped': 'Stopped',
        'status_starting': 'Starting',
        'error_input_not_aligned': "The provided file does not appear to be a multiple sequence alignment: sequences have different lengths (or only one sequence provided). Please align your sequences (e.g., using MAFFT or Clustal) before submitting.",
        'msg_signatures_found': "✅ {count} signature(s) found",
        'msg_completed_no_results': "Execution completed: no signature found with current parameters. Try relaxing thresholds (coverage, degeneracy) or checking input alignment.",
        'msg_completed_unknown': "ℹ️ Execution completed: signature status could not be determined from logs.",
        'msg_exec_error': "❌ Execution error (code {code})",
        'msg_error_detail': "\nDetail: {detail}",
        'sugg_memory': "\n💡 Suggestion: Try with a smaller file or less strict parameters",
        'sugg_file': "\n💡 Suggestion: Check that the FASTA file exists and is accessible",
        'sugg_perms': "\n💡 Suggestion: Permission issue, contact administrator",
        'full_id': "Full ID",
        'created': "Created",
        'started': "Started",
        'ended': "Finished",
        'duration': "Duration",
        'log_lines': "Log lines",
        'command': "Command",
        'stop_confirm': "Stop execution?",
        'help_tip_1': "Check that your FASTA file is correctly formatted",
        'help_tip_2': "Try using more permissive parameters",
        'help_tip_3': "Contact support if the issue persists",
        'progress_init': "Initializing...",
        'progress_done': "Completed",
        'total_lines_gen': "total lines / {total} generated",
        'lines_count': "lines",
        'view_btn': "View",
        'stop_btn': "Stop",
        'results_btn': "Results",
        'footer_text': 'LAVA-DNA Flask Interface - Stable interface for LAMP primer design'
    }
}

app = Flask(__name__)
flask_env = os.environ.get('FLASK_ENV', 'development')
secret_key = os.environ.get('SECRET_KEY')
if not secret_key:
    if flask_env == 'production':
        raise RuntimeError("ERREUR CRITIQUE DE SECURITE : Variable d'environnement SECRET_KEY obligatoire en mode production.")
    secret_key = os.urandom(24)
app.secret_key = secret_key
app.config['MAX_CONTENT_LENGTH'] = 1 * 1024 * 1024 * 1024  # 1GB max file size

# Securisation des cookies de session (Priorite 2)
app.config.update(
    SESSION_COOKIE_HTTPONLY=True,
    SESSION_COOKIE_SAMESITE='Lax'
)
if flask_env == 'production' or os.environ.get('SESSION_COOKIE_SECURE', 'False').lower() in ('true', '1', 't'):
    app.config['SESSION_COOKIE_SECURE'] = True

import uuid
import time as time_module
from collections import defaultdict

# Rate limiter leger par IP (Priorite 5)
ip_request_history = defaultdict(list)

def check_rate_limit(max_requests=15, window_seconds=60):
    client_ip = request.remote_addr or 'unknown'
    now = time_module.time()
    ip_request_history[client_ip] = [t for t in ip_request_history[client_ip] if now - t < window_seconds]
    if len(ip_request_history[client_ip]) >= max_requests:
        return False
    ip_request_history[client_ip].append(now)
    return True

def check_execution_ownership(execution_id):
    """Verifie que l'execution appartient bien au visiteur courant"""
    if execution_id not in running_executions:
        from flask import abort
        abort(404)
    execution = running_executions[execution_id]
    owner_id = execution.get('owner_id')
    if owner_id and owner_id != session.get('user_id'):
        from flask import abort
        abort(403)

@app.errorhandler(500)
def handle_internal_error(error):
    """Masquer les traces techniques en production"""
    if os.environ.get('FLASK_ENV') == 'production':
        return "Erreur interne du serveur. Veuillez reessayer ou contacter l'administrateur.", 500
    return str(error), 500

@app.after_request
def add_header(response):
    """Empecher la mise en cache par le navigateur pour garantir que le changement de langue s'affiche immediatement."""
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '-1'
    return response

# Ajouter le filtre basename pour les templates
@app.template_filter('basename')
def basename_filter(path):
    return os.path.basename(path)

@app.before_request
def before_request():
    """Definir la langue globale pour la requete et attribuer un user_id anonyme"""
    if 'user_id' not in session:
        session['user_id'] = str(uuid.uuid4())
    g.lang = session.get('language', 'fr')
    if g.lang not in TRANSLATIONS:
        g.lang = 'fr'

@app.template_filter('t')
@pass_context
def translate_filter(context, key):
    """Filtre Jinja2 pour les traductions"""
    lang = getattr(g, 'lang', None) or session.get('language', 'fr')
    if lang not in TRANSLATIONS:
        lang = 'fr'
    if key is None:
        return ""
    return TRANSLATIONS.get(lang, {}).get(key, key)

@app.context_processor
def inject_globals():
    """Injecter automatiquement la langue dans tous les templates"""
    lang = getattr(g, 'lang', None) or session.get('language', 'fr')
    if lang not in TRANSLATIONS:
        lang = 'fr'
    return {
        'lang': lang
    }

# Configuration
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'
ALLOWED_EXTENSIONS = {'fas', 'fasta', 'fa', 'txt'}

# Créer les dossiers nécessaires
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs('templates', exist_ok=True)
os.makedirs('static', exist_ok=True)

# Stockage des exécutions en cours
running_executions = {}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# Ensemble des paramètres flottants connus (utilisé par _convert_param_value)
# Known float parameters set (used by _convert_param_value)
FLOAT_PARAMS = {
    'min_base_frequency', 'entropy_threshold', 'penalty_plateau', 
    'penalty_slope', 'dntp_conc', 'dna_conc', 'salt_monovalent', 
    'salt_divalent', 'max_tm_diff', 'max_primer_gen',
    'outer_primer_target_tm', 'outer_primer_min_tm', 'outer_primer_max_tm',
    'middle_primer_target_tm', 'middle_primer_min_tm', 'middle_primer_max_tm',
    'inner_primer_target_tm', 'inner_primer_min_tm', 'inner_primer_max_tm',
    'loop_primer_target_tm', 'loop_primer_min_tm', 'loop_primer_max_tm',
    'stem_primer_target_tm', 'stem_primer_min_tm', 'stem_primer_max_tm',
    'primer_min_match_percent', 'primer_iupac_min_percent', 
    'min_primer_coverage', 'max_overlap_percent'
}

def _convert_param_value(key, value):
    """Convertit une valeur de paramètre au bon type (float, int ou str).
    Converts a parameter value to the correct type (float, int, or str)."""
    is_float = (key in FLOAT_PARAMS or 
                any(hint in key for hint in ('tm', 'percent', 'coverage', 'conc', 
                                              'salt', 'frequency', 'slope', 
                                              'plateau', 'threshold')))
    if is_float:
        try:
            return float(value)
        except ValueError:
            return value
    else:
        try:
            return int(value)
        except ValueError:
            return value

def _apply_lamp_mode(params, lamp_mode, script_type):
    """Applique le mode LAMP (classic/enriched) aux paramètres include_*.
    Applies LAMP mode (classic/enriched) to include_* parameters."""
    params['lamp_mode'] = lamp_mode
    if lamp_mode == 'enriched':
        if script_type == 'STEM':
            params['include_stem_primers'] = True
            params['include_loop_primers'] = False
        else:
            params['include_loop_primers'] = True
            params['include_stem_primers'] = False
    else:  # classic
        params['include_stem_primers'] = False
        params['include_loop_primers'] = False

def get_language():
    """Récupérer la langue actuelle (via g défini dans before_request)"""
    return getattr(g, 'lang', 'fr')

def get_text(key):
    """Récupérer un texte traduit selon la langue actuelle"""
    lang = get_language()
    return TRANSLATIONS.get(lang, {}).get(key, key)

def get_default_params():
    """Paramètres par défaut pour LAVA"""
    return {
        'script_type': 'STEM',
        'lamp_mode': 'classic',  # Par défaut mode classique
        'signature_max_length': 400,
        'max_primer_gen': 5000,
        'primer_min_match_percent': 80,
        'primer_iupac_min_percent': 90,
        'min_primer_coverage': 80,

        'min_base_frequency': 0.05,
        'entropy_threshold': 1.5,
        'max_total_degenerate_bases': 2,
        'max_consecutive_degenerate_bases': 2,
        'max_3prime_degenerate_bases': 2,
        'max_tolerated_mismatches': 0,
        'three_prime_zone_size': 5,
        'max_poly_bases': 2,
        'min_signatures_for_success': 1,
        'max_overlap_percent': 0,
        # Reduction spatiale par fenetre / Spatial window reduction
        'window_size': 0,        # 0 = desactive, ex: 5 = fenetre de 5nt
        'max_per_window': 0,     # 0 = desactive, ex: 3 = 3 candidats max par fenetre
        
        # Thermodynamics & Conditions
        'max_tm_diff': 5.0,
        'dntp_conc': 1.4,
        'dna_conc': 400.0,
        'salt_monovalent': 50.0,
        'salt_divalent': 8.0,

        # Architecture specific
        'max_dist_outer_middle': 30,
        'max_dist_middle_inner': 30,
        'penalty_plateau': 0.25,
        'penalty_slope': 0.15,
        
        # Outer primers
        'outer_primer_target_length': 20,
        'outer_primer_min_length': 18,
        'outer_primer_max_length': 22,
        'outer_primer_target_tm': 58.0,
        'outer_primer_min_tm': 57.0,
        'outer_primer_max_tm': 59.0,
        
        # Middle primers
        'middle_primer_target_length': 20,
        'middle_primer_min_length': 18,
        'middle_primer_max_length': 22,
        'middle_primer_target_tm': 60.0,
        'middle_primer_min_tm': 59.0,
        'middle_primer_max_tm': 65.0,
        
        # Inner primers
        'inner_primer_target_length': 20,
        'inner_primer_min_length': 15,
        'inner_primer_max_length': 22,
        'inner_primer_target_tm': 60.0,
        'inner_primer_min_tm': 59.0,
        'inner_primer_max_tm': 65.0,
        
        # STEM/LOOP specific
        'stem_primer_target_length': 20,
        'stem_primer_min_length': 18,
        'stem_primer_max_length': 22,
        'stem_primer_target_tm': 60.0,
        'stem_primer_min_tm': 59.0,
        'stem_primer_max_tm': 61.0,
        'include_stem_primers': True,
        
        'loop_primer_target_length': 20,
        'loop_primer_min_length': 18,
        'loop_primer_max_length': 22,
        'loop_primer_target_tm': 60.0,
        'loop_primer_min_tm': 59.0,
        'loop_primer_max_tm': 61.0,
        'loop_min_gap': 25,
        'include_loop_primers': True
    }

@app.route('/')
def index():
    """Page d'accueil"""
    if 'params' not in session:
        session['params'] = get_default_params()
    
    return render_template('index.html', 
                         params=session['params'],
                         running_count=len(running_executions))

@app.route('/set_language/<language>')
def set_language(language):
    """Changer la langue de l'interface / Change the interface language"""
    next_url = request.args.get('next') or url_for('index')
    
    if language in TRANSLATIONS:
        session['language'] = language
        session.modified = True
        g.lang = language
        print(f"🌐 Langue changee vers: {language}")
    else:
        print(f"❌ Langue non supportee: {language}")
        
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta http-equiv="refresh" content="0;url={next_url}">
</head>
<body>
    <script>window.location.replace("{next_url}");</script>
</body>
</html>"""
    resp = make_response(html_content, 200)
    resp.headers['Content-Type'] = 'text/html; charset=utf-8'
    resp.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    resp.headers['Pragma'] = 'no-cache'
    resp.headers['Expires'] = '0'
    return resp

@app.route('/upload', methods=['POST'])
def upload_file():
    """Upload du fichier FASTA"""
    if not check_rate_limit(max_requests=15, window_seconds=60):
        flash('Trop de requetes. Veuillez patienter une minute.', 'error')
        return redirect(url_for('index'))
    if 'fasta_file' not in request.files:
        flash('Aucun fichier sélectionné', 'error')
        return redirect(url_for('index'))
    
    file = request.files['fasta_file']
    if file.filename == '':
        flash('Aucun fichier sélectionné', 'error')
        return redirect(url_for('index'))
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        unique_filename = f"{timestamp}_{filename}"
        filepath = os.path.join(UPLOAD_FOLDER, unique_filename)
        file.save(filepath)
        
        session['uploaded_file'] = filepath
        session['uploaded_filename'] = filename
        flash(f'Fichier {filename} uploadé avec succès', 'success')
    else:
        flash('Type de fichier non autorisé. Utilisez .fas, .fasta, .fa ou .txt', 'error')
    
    return redirect(url_for('index'))

@app.route('/update_params', methods=['POST'])
def update_params():
    """Mise à jour des paramètres"""
    if 'params' not in session:
        session['params'] = get_default_params()
    
    # Traiter script_type en PREMIER
    if 'script_type' in request.form:
        session['params']['script_type'] = request.form['script_type']
    
    # Traiter lamp_mode en SECOND (après script_type)
    if 'lamp_mode' in request.form:
        _apply_lamp_mode(session['params'], request.form['lamp_mode'], 
                        session['params']['script_type'])
    
    # Traiter tous les autres paramètres du formulaire
    # Process all other form parameters
    for key, value in request.form.items():
        if key not in ['script_type', 'lamp_mode']:
            session['params'][key] = _convert_param_value(key, value)
    
    session.modified = True
    flash(get_text('params_updated'), 'success')
    return redirect(url_for('index'))

def translate_error_to_user_friendly(error_message, lang='fr'):
    """Traduit les erreurs techniques en messages compréhensibles pour l'utilisateur"""
    error_str = str(error_message).lower()
    is_en = (lang == 'en')
    
    # Erreurs d'encodage
    if "utf-8" in error_str and "decode" in error_str:
        if "0xc3" in error_str or "continuation byte" in error_str:
            if is_en:
                return ("❌ FILE ENCODING PROBLEM\n\n"
                       "The FASTA file contains unsupported special characters.\n\n"
                       "SOLUTIONS:\n"
                       "• Open your file in a text editor (Notepad++, VSCode)\n"
                       "• Save it with 'UTF-8' encoding\n"
                       "• Check that sequence names contain no accented characters")
            return ("❌ PROBLÈME D'ENCODAGE DU FICHIER\n\n"
                   "Le fichier FASTA contient des caractères spéciaux non compatibles.\n\n"
                   "SOLUTIONS :\n"
                   "• Ouvrez votre fichier dans un éditeur de texte (Notepad++, VSCode)\n"
                   "• Sauvegardez-le avec l'encodage 'UTF-8'\n"
                   "• Vérifiez qu'il n'y a pas de caractères accentués dans les noms de séquences\n"
                   "• Remplacez les caractères comme à, é, è, ç par a, e, e, c")
        else:
            if is_en:
                return ("❌ ENCODING ERROR\n\n"
                       "The file cannot be read properly.\n"
                       "Ensure it is saved as plain text (UTF-8).")
            return ("❌ ERREUR D'ENCODAGE\n\n"
                   "Le fichier ne peut pas être lu correctement.\n"
                   "Assurez-vous qu'il est en format texte simple (UTF-8).")
    
    # Erreurs de fichier non trouvé
    if "no such file" in error_str or "file not found" in error_str:
        if is_en:
            return ("❌ FILE NOT FOUND\n\n"
                   "The FASTA file could not be located.\n"
                   "Please upload it again.")
        return ("❌ FICHIER NON TROUVÉ\n\n"
               "Le fichier FASTA n'a pas pu être localisé.\n"
               "Veuillez le télécharger à nouveau.")
    
    # Erreurs de permissions
    if "permission denied" in error_str:
        if is_en:
            return ("❌ PERMISSION ERROR\n\n"
                   "The system does not have permission to read the file.\n"
                   "Contact administrator.")
        return ("❌ ERREUR DE PERMISSIONS\n\n"
               "Le système n'a pas les droits pour lire le fichier.\n"
               "Contactez l'administrateur.")
    
    # Erreurs de format FASTA
    if "fasta" in error_str and ("format" in error_str or "invalid" in error_str):
        if is_en:
            return ("❌ INVALID FASTA FORMAT\n\n"
                   "The file does not respect standard FASTA format.\n\n"
                   "EXPECTED FORMAT:\n"
                   ">sequence_name1\n"
                   "ATCGATCG...\n"
                   ">sequence_name2\n"
                   "GCTAGCTA...")
        return ("❌ FORMAT FASTA INVALIDE\n\n"
               "Le fichier ne respecte pas le format FASTA.\n\n"
               "FORMAT ATTENDU :\n"
               ">nom_sequence1\n"
               "ATCGATCG...\n"
               ">nom_sequence2\n"
               "GCTAGCTA...")
    
    # Erreurs de mémoire
    if "memory" in error_str or "memoryerror" in error_str:
        if is_en:
            return ("❌ INSUFFICIENT MEMORY\n\n"
                   "The file is too large to process.\n"
                   "Try again with a smaller file.")
        return ("❌ MÉMOIRE INSUFFISANTE\n\n"
               "Le fichier est trop volumineux pour être traité.\n"
               "Essayez avec un fichier plus petit.")
    
    # Erreurs Perl spécifiques
    if "can't locate" in error_str and ".pm" in error_str:
        if is_en:
            return ("❌ LAVA CONFIGURATION ERROR\n\n"
                   "A required Perl module is missing.\n"
                   "Contact administrator to reinstall LAVA.")
        return ("❌ ERREUR DE CONFIGURATION LAVA\n\n"
               "Un module Perl requis est manquant.\n"
               "Contactez l'administrateur pour réinstaller LAVA.")
    
    # Erreurs de paramètres LAVA
    if "primer" in error_str and ("length" in error_str or "target" in error_str):
        if is_en:
            return ("❌ INVALID PARAMETERS\n\n"
                   "Primer lengths or target settings are incorrect.\n"
                   "Verify that:\n"
                   "• Min length < Max length\n"
                   "• Values are positive integers")
        return ("❌ PARAMÈTRES INVALIDES\n\n"
               "Les longueurs de primers ou les cibles sont incorrectes.\n"
               "Vérifiez que :\n"
               "• La longueur minimum < longueur maximum\n"
               "• Les valeurs sont des nombres entiers positifs")
    
    # Si l'erreur n'est pas reconnue, retourner un message générique plus utile
    if is_en:
        return (f"❌ TECHNICAL ERROR\n\n"
               f"Original error message:\n{error_message}\n\n"
               f"RECOMMENDED ACTIONS:\n"
               f"• Check your FASTA file format\n"
               f"• Try with a smaller file\n"
               f"• Contact support with this error message")
    return (f"❌ ERREUR TECHNIQUE\n\n"
           f"Message d'erreur original :\n{error_message}\n\n"
           f"ACTIONS RECOMMANDÉES :\n"
           f"• Vérifiez le format de votre fichier FASTA\n"
           f"• Essayez avec un fichier plus petit\n"
           f"• Contactez le support avec ce message d'erreur")

def execute_lava_background(execution_id, script_type, input_file, output_name, params):
    """Exécution LAVA en arrière-plan"""
    try:
        # Préparer la commande
        script_name = f"lava_{script_type.lower()}_primer.pl"
        output_base = os.path.join(RESULTS_FOLDER, output_name)
        output_file = f"{output_base}.primers"
        
        cmd = [
            "perl", script_name,
            "--alignment_fasta", input_file,
            "--output_file", output_file
        ]
        
        # Mapping des noms de paramètres Flask vers les noms Perl
        param_mapping = {
            'iupac_match_percent': 'primer_iupac_min_percent',
            'minimum_primer_coverage': 'min_primer_coverage', 
            'minimum_signature_coverage': 'min_signature_coverage',
            'mismatch_tolerance': 'primer_min_match_percent',
            'signature_max_length': 'signature_max_length',
            'max_primer_gen': 'max_primer_gen',
            'primer_min_match_percent': 'primer_min_match_percent',
            'primer_iupac_min_percent': 'primer_iupac_min_percent',
            'min_primer_coverage': 'min_primer_coverage',
            'min_base_frequency': 'min_base_frequency'
        }
        
        # Paramètres valides pour les scripts Perl (pour filtrer les invalides)
        # Paramètres communs (utilisés par LOOP et STEM)
        common_params = {
            'signature_max_length', 'max_primer_gen', 'primer_min_match_percent',
            'primer_iupac_min_percent', 'min_primer_coverage', 'min_base_frequency',
            'min_signatures_for_success', 'max_overlap_percent', 'resolve_overlap_by',
            'primer3_executable', 'thermodynamic_path', 'alignment_format',
            'dntp_conc', 'dna_conc', 'salt_monovalent', 'salt_divalent',
            'max_poly_bases', 'entropy_threshold', 
            'max_total_degenerate_bases', 'max_consecutive_degenerate_bases',
            'max_3prime_degenerate_bases', 'three_prime_zone_size',
            'max_tolerated_mismatches',
            'penalty_plateau', 'penalty_slope', 'max_tm_diff',
            'outer_pair_target_length', 'middle_pair_target_length', 'inner_pair_target_length',
             # Outer primers
            'outer_primer_target_length', 'outer_primer_min_length', 'outer_primer_max_length', 
            'outer_primer_target_tm', 'outer_primer_min_tm', 'outer_primer_max_tm',
            # Middle primers  
            'middle_primer_target_length', 'middle_primer_min_length', 'middle_primer_max_length',
            'middle_primer_target_tm', 'middle_primer_min_tm', 'middle_primer_max_tm',
            # Inner primers
            'inner_primer_target_length', 'inner_primer_min_length', 'inner_primer_max_length',
            'inner_primer_target_tm', 'inner_primer_min_tm', 'inner_primer_max_tm',
            # Calcul dynamique des longueurs (commun STEM+LOOP depuis Phase 36)
            # Dynamic length calculation (common STEM+LOOP since Phase 36)
            'max_dist_outer_middle', 'max_dist_middle_inner',
            # Reduction spatiale par fenetre / Spatial window reduction
            'window_size', 'max_per_window',
        }

        # Paramètres spécifiques à LOOP
        loop_only_params = {
            'include_loop_primers', 'loop_min_gap', 'min_primer_spacing', 'min_inner_pair_spacing',
            'loop_primer_target_length', 'loop_primer_min_length', 'loop_primer_max_length',
            'loop_primer_target_tm', 'loop_primer_min_tm', 'loop_primer_max_tm',
        }

        # Paramètres spécifiques à STEM
        stem_only_params = {
            'include_stem_primers', 'min_primer_spacing', 'min_inner_pair_spacing',
            'stem_primer_target_length', 'stem_primer_min_length', 'stem_primer_max_length',
            'stem_primer_target_tm', 'stem_primer_min_tm', 'stem_primer_max_tm',
        }
        
        # Sélectionner les paramètres valides selon le script
        valid_perl_params = common_params.copy()
        if script_type.upper() == 'LOOP':
            valid_perl_params.update(loop_only_params)
        elif script_type.upper() == 'STEM':
            valid_perl_params.update(stem_only_params)

        
        # Ajouter les paramètres
        for param_name, param_value in params.items():
            if param_value is not None and param_name not in ['script_type', 'lamp_mode']:
                # Convertir les booléens en entiers pour les scripts Perl
                if isinstance(param_value, bool):
                    param_value = 1 if param_value else 0
                
                # Utiliser le nom mappé ou le nom original
                perl_param_name = param_mapping.get(param_name, param_name)
                
                # Seulement ajouter si c'est un paramètre valide pour Perl
                if perl_param_name in valid_perl_params:
                    cmd.extend([f"--{perl_param_name}", str(param_value)])
                else:
                    print(f"⚠️  Paramètre ignoré (non supporté par Perl): {param_name} -> {perl_param_name}")
        
        # Sauvegarder les paramètres dans un fichier texte pour la traçabilité
        params_file_path = f"{output_base}.params.txt"
        with open(params_file_path, 'w') as pf:
            pf.write(f"=== PARAMÈTRES D'EXÉCUTION LAVA ===\n")
            pf.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            pf.write(f"Script: {script_name}\n")
            pf.write(f"Input Fasta: {input_file}\n")
            pf.write(f"Output Base: {output_base}\n")
            pf.write(f"Mode LAMP: {params.get('lamp_mode', 'classic')}\n")
            pf.write(f"\n--- Paramètres Perl passés ---\n")
            
            # Récupérer les paires de paramètres depuis cmd
            for i in range(2, len(cmd), 2):
                if i + 1 < len(cmd) and cmd[i].startswith('--'):
                    pf.write(f"{cmd[i]}: {cmd[i+1]}\n")

        # Variables d'environnement
        env = os.environ.copy()
        env['PERL5LIB'] = './lib'
        
        # Mettre à jour le statut
        running_executions[execution_id]['status'] = 'running'
        running_executions[execution_id]['command'] = ' '.join(cmd)
        running_executions[execution_id]['start_time'] = datetime.now()
        
        # Lancer le processus
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=os.getcwd(),
            env=env,
            bufsize=1,
            universal_newlines=True
        )
        
        running_executions[execution_id]['process'] = process
        
        # Lire les logs avec capture COMPLÈTE (pas de troncature pendant l'exécution)
        output_lines = []
        buffer_size = 1000  # Taille du buffer pour l'affichage temps réel
        
        while True:
            line = process.stdout.readline()
            if line:
                stripped = line.strip()
                
                # ── Détection des lignes de progression LAVA-PROGRESS ─────────────
                # Format Perl : [LAVA-PROGRESS] label|done|total|extra|rate|eta
                # Detect LAVA-PROGRESS lines emitted by Perl scripts
                if stripped.startswith('[LAVA-PROGRESS]'):
                    try:
                        parts = stripped[len('[LAVA-PROGRESS]'):].strip().split('|')
                        if len(parts) >= 3:
                            running_executions[execution_id]['progress'] = {
                                'label': parts[0].strip(),
                                'done':  int(parts[1]),
                                'total': int(parts[2]),
                                'extra': parts[3].strip() if len(parts) > 3 else '',
                                'rate':  parts[4].strip() if len(parts) > 4 else '',
                                'eta':   int(parts[5]) if len(parts) > 5 and parts[5].strip().isdigit() else 0,
                                'pct':   round(int(parts[1]) / max(int(parts[2]), 1) * 100, 1),
                            }
                    except Exception:
                        pass  # Ligne malformée ignorée silencieusement
                    # NE PAS ajouter ces lignes aux logs affichés
                    continue
                # ─────────────────────────────────────────────────────────────────
                
                output_lines.append(stripped)
                # Pour l'affichage temps réel, garder un mix : début + fin récente
                if len(output_lines) > buffer_size:
                    # Garder les 100 premières + les 400 dernières lignes pour l'affichage
                    display_logs = output_lines[:100] + ["... [logs intermédiaires masqués] ..."] + output_lines[-(buffer_size-101):]
                else:
                    display_logs = output_lines
                
                running_executions[execution_id]['logs'] = display_logs
                running_executions[execution_id]['total_lines'] = len(output_lines)
            else:
                # Pas de nouvelle ligne, vérifier si le processus est terminé
                if process.poll() is not None:
                    break
                # Petite pause pour éviter une boucle CPU intensive
                import time
                time.sleep(0.01)
        
        # Processus terminé - capturer tout ce qui reste
        import time
        time.sleep(0.2)  # Pause plus longue pour vider les buffers
        
        remaining_lines = process.stdout.readlines()
        for remaining_line in remaining_lines:
            if remaining_line.strip():
                output_lines.append(remaining_line.strip())
        
        # Maintenant stocker TOUS les logs pour l'affichage final
        running_executions[execution_id]['logs'] = output_lines  # TOUS les logs
        running_executions[execution_id]['total_lines'] = len(output_lines)
        
        # Récupérer le code de retour
        return_code = process.wait()
        running_executions[execution_id]['return_code'] = return_code
        running_executions[execution_id]['end_time'] = datetime.now()
        
        # Analyser les logs pour déterminer le vrai statut
        all_logs_raw = '\n'.join(output_lines)
        all_logs_text = all_logs_raw.lower()
        
        if return_code == 0:
            # Chercher TOUS les fichiers de résultats générés
            import glob
            result_files = []
            
            # Chercher tous les fichiers qui commencent par le nom de base
            base_pattern = f"{output_base}*"
            potential_files = glob.glob(base_pattern)
            
            # Filtrer et organiser les fichiers
            for file_path in potential_files:
                if os.path.isfile(file_path):
                    # Exclure les fichiers temporaires ou non désirés
                    filename = os.path.basename(file_path)
                    if not filename.endswith('.tmp') and not filename.startswith('.'):
                        result_files.append(file_path)
            
            # Trier par nom pour un affichage ordonné
            result_files.sort()
            
            # Chercher le nombre de signatures trouvées dans les logs
            # Format attendu : "After reduction: N final signatures"
            sig_matches = re.findall(r'after reduction:\s*(\d+)\s*final signatures', all_logs_text, re.IGNORECASE)
            
            user_lang = running_executions[execution_id].get('lang', 'fr')
            t_dict = TRANSLATIONS.get(user_lang, TRANSLATIONS['fr'])
            
            if sig_matches:
                sig_count = int(sig_matches[-1])
                running_executions[execution_id]['signature_count'] = sig_count
                if sig_count > 0:
                    running_executions[execution_id]['status'] = 'completed'
                    running_executions[execution_id]['completion_message'] = t_dict.get('msg_signatures_found').format(count=sig_count)
                else:
                    running_executions[execution_id]['status'] = 'completed_no_results'
                    running_executions[execution_id]['completion_message'] = t_dict.get('msg_completed_no_results')
            else:
                running_executions[execution_id]['status'] = 'completed'
                running_executions[execution_id]['completion_message'] = t_dict.get('msg_completed_unknown')
            
            running_executions[execution_id]['result_files'] = result_files
        else:
            # Code de retour non-zéro = erreur
            user_lang = running_executions[execution_id].get('lang', 'fr')
            t_dict = TRANSLATIONS.get(user_lang, TRANSLATIONS['fr'])
            if 'input_not_aligned' in all_logs_text:
                error_msg = t_dict.get('error_input_not_aligned')
            else:
                error_msg = t_dict.get('msg_exec_error').format(code=return_code)
                
                # Remonter la dernière ligne significative d'erreur depuis les logs
                significant_error_line = None
                error_markers = ['error', 'erreur', 'died', "can't", 'at ', ' line ', 'undefined', 'confess', 'fatal', 'exception']
                for line in reversed(output_lines):
                    line_clean = line.strip()
                    if not line_clean:
                        continue
                    line_lower = line_clean.lower()
                    if any(marker in line_lower for marker in error_markers):
                        significant_error_line = line_clean
                        break
                
                if significant_error_line:
                    error_msg += t_dict.get('msg_error_detail').format(detail=significant_error_line)
                
                # Analyser les logs pour donner plus d'informations et suggestions
                if 'out of memory' in all_logs_text or 'memory' in all_logs_text:
                    error_msg += t_dict.get('sugg_memory')
                elif 'no such file' in all_logs_text:
                    error_msg += t_dict.get('sugg_file')
                elif 'permission denied' in all_logs_text:
                    error_msg += t_dict.get('sugg_perms')
            
            running_executions[execution_id]['status'] = 'error'
            running_executions[execution_id]['error'] = error_msg
    
    except Exception as e:
        user_lang = running_executions[execution_id].get('lang', 'fr')
        running_executions[execution_id]['status'] = 'error'
        user_friendly_error = translate_error_to_user_friendly(e, user_lang)
        running_executions[execution_id]['error'] = user_friendly_error
        # Garder aussi l'erreur technique pour le debug
        running_executions[execution_id]['technical_error'] = str(e)

@app.route('/execute', methods=['POST'])
def execute_lava():
    """Lancer l'exécution LAVA"""
    if not check_rate_limit(max_requests=10, window_seconds=60):
        flash("Trop de requêtes. Veuillez patienter une minute avant de lancer un nouveau calcul.", "error")
        return redirect(url_for('index'))

    if 'uploaded_file' not in session:
        flash('Aucun fichier FASTA uploadé', 'error')
        return redirect(url_for('index'))
        
    # Vérification des quotas de concurrence
    max_global = int(os.environ.get('MAX_CONCURRENT_RUNS', 5))
    max_user = int(os.environ.get('MAX_USER_CONCURRENT_RUNS', 2))
    current_user_id = session.get('user_id')

    running_global = sum(1 for e in running_executions.values() if e.get('status') in ['starting', 'running'])
    running_user = sum(1 for e in running_executions.values() if (e.get('owner_id') == current_user_id or e.get('user_id') == current_user_id) and e.get('status') in ['starting', 'running'])

    if running_user >= max_user:
        flash(f"Vous avez déjà {running_user} calculs en cours. Veuillez attendre leur fin avant d'en lancer un nouveau.", "warning")
        return redirect(url_for('list_executions'))

    if running_global >= max_global:
        flash("Le serveur est actuellement à pleine capacité. Veuillez réessayer dans quelques instants.", "warning")
        return redirect(url_for('list_executions'))
    
    # Récupérer les paramètres du formulaire actuel
    script_type = request.form.get('script_type', 'STEM').upper()
    
    # Liste Blanche
    if script_type not in ['STEM', 'LOOP']:
        flash('Type de script invalide (seuls STEM ou LOOP sont autorisés).', 'error')
        return redirect(url_for('index'))
        
    output_name = request.form.get('output_name', 'lava_result')
    if not output_name:
        output_name = 'lava_result'
        
    # S'assurer que session['params'] existe
    if 'params' not in session:
        session['params'] = get_default_params()
    else:
        # Si le type de script a changé, on réinitialise ou on garde les paramètres communs
        if session['params'].get('script_type') != script_type:
            session['params'] = get_default_params()
        else:
            # S'assurer que tous les paramètres existent
            defaults = get_default_params()
            for k, v in defaults.items():
                if k not in session['params']:
                    session['params'][k] = v
                    
    session['params']['script_type'] = script_type
    if 'lamp_mode' in request.form:
        _apply_lamp_mode(session['params'], request.form['lamp_mode'], script_type)
    elif 'lamp_mode' not in session['params']:
        session['params']['lamp_mode'] = 'classic'
    
    for key, value in request.form.items():
        if key not in ['script_type', 'lamp_mode', 'output_name']:
            if key in ['include_stem_primers', 'include_loop_primers']:
                session['params'][key] = value == 'on'
            else:
                session['params'][key] = _convert_param_value(key, value)
    

    execution_id = str(uuid.uuid4())
    
    running_executions[execution_id] = {
        'id': execution_id,
        'user_id': session.get('user_id'),
        'owner_id': session.get('user_id'),
        'lang': getattr(g, 'lang', None) or session.get('language', 'fr'),
        'status': 'starting',
        'input_file': session['uploaded_file'],
        'output_name': output_name,
        'script_type': script_type,
        'logs': [],
        'total_lines': 0,
        'created_time': datetime.now()
    }
    
    # Lancer l'exécution en arrière-plan
    thread = threading.Thread(
        target=execute_lava_background,
        args=(execution_id, script_type, 
              session['uploaded_file'], output_name, session['params'])
    )
    thread.daemon = True
    thread.start()
    
    flash(f'Exécution LAVA lancée (ID: {execution_id[:8]})', 'success')
    return redirect(url_for('monitor', execution_id=execution_id))

@app.route('/monitor/<execution_id>')
def monitor(execution_id):
    """Page de monitoring d'une execution"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    return render_template('monitor.html', execution=execution)

@app.route('/api/status/<execution_id>')
def api_status(execution_id):
    """API pour recuperer le statut d'une execution"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    
    # Paramètre pour demander tous les logs
    show_all = request.args.get('all_logs', 'false').lower() == 'true'
    
    # Formatter les timestamps
    # Choisir le nombre de logs à envoyer
    logs = execution.get('logs', [])
    
    if show_all:
        # Mode "tous les logs" - limiter à 2000 lignes max pour éviter les problèmes de performance
        if len(logs) > 2000:
            # Envoyer début + fin si trop long
            logs_to_send = logs[:500] + ["... [logs intermédiaires masqués pour performance] ..."] + logs[-1500:]
        else:
            logs_to_send = logs  # Tous les logs si raisonnable
    else:
        # Mode normal - les 100 dernières lignes
        logs_to_send = logs[-100:]
    
    data = {
        'id': execution['id'],
        'status': execution['status'],
        'total_lines': execution.get('total_lines', 0),
        'logs': logs_to_send,
        # Données de progression LAVA-PROGRESS (None si aucune en cours)
        # LAVA-PROGRESS data (None if none currently active)
        'progress': execution.get('progress', None),
    }
    
    if 'start_time' in execution:
        data['start_time'] = execution['start_time'].strftime('%H:%M:%S')
    
    if 'end_time' in execution:
        data['end_time'] = execution['end_time'].strftime('%H:%M:%S')
        duration = execution['end_time'] - execution['start_time']
        data['duration'] = str(duration).split('.')[0]  # Supprimer les microsecondes
    
    if 'error' in execution:
        data['error'] = execution['error']
    
    if 'technical_error' in execution:
        data['technical_error'] = execution['technical_error']
    
    if 'completion_message' in execution:
        data['completion_message'] = execution['completion_message']
    
    if 'result_files' in execution:
        # Créer des informations détaillées sur chaque fichier
        result_files_info = []
        for file_path in execution['result_files']:
            if os.path.exists(file_path):
                file_stat = os.stat(file_path)
                file_size_mb = file_stat.st_size / (1024 * 1024)
                
                # Déterminer le type de fichier
                filename = os.path.basename(file_path)
                if filename.endswith('.primers'):
                    file_type = 'Primers principaux'
                    file_icon = 'fas fa-dna'
                elif 'all_signatures' in filename:
                    file_type = 'Toutes les signatures'
                    file_icon = 'fas fa-list'
                elif filename.endswith('.dash'):
                    file_type = 'Format tableau'
                    file_icon = 'fas fa-table'
                elif 'amplifiees' in filename and filename.endswith('.fasta'):
                    file_type = 'Séquences amplifiées'
                    file_icon = 'fas fa-plus-circle'
                elif 'exclues' in filename and filename.endswith('.fasta'):
                    file_type = 'Séquences exclues'
                    file_icon = 'fas fa-minus-circle'
                elif 'noms' in filename and filename.endswith('.txt'):
                    file_type = 'Noms des séquences'
                    file_icon = 'fas fa-file-alt'
                else:
                    file_type = 'Résultat LAVA'
                    file_icon = 'fas fa-file'
                
                result_files_info.append({
                    'name': filename,
                    'path': file_path,
                    'size_mb': round(file_size_mb, 2),
                    'size_bytes': file_stat.st_size,
                    'type': file_type,
                    'icon': file_icon
                })
        
        data['result_files'] = result_files_info
    
    # Calculer la taille estimée des logs
    if logs:
        estimated_size_mb = len('\n'.join(logs).encode('utf-8')) / (1024 * 1024)
        data['logs_size_mb'] = round(estimated_size_mb, 1)
    
    return jsonify(data)

@app.route('/download-logs/<execution_id>')
def download_logs(execution_id):
    """Telecharger tous les logs d'une execution"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    logs = execution.get('logs', [])
    
    if not logs:
        flash('Aucun log disponible', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    # Calculer la taille estimée
    log_content = '\n'.join(logs)
    estimated_size_mb = len(log_content.encode('utf-8')) / (1024 * 1024)
    
    script_type = execution.get('script_type', 'UNKNOWN')
    timestamp = execution.get('created_time', datetime.now()).strftime('%Y%m%d_%H%M%S')
    
    # Gérer les gros fichiers avec compression
    compress = request.args.get('compress', 'auto')
    
    if compress == 'auto':
        # Auto-compression si > 5MB
        should_compress = estimated_size_mb > 5.0
    else:
        should_compress = compress.lower() == 'true'
    
    if should_compress:
        # Compression gzip
        import gzip
        import io
        
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode='wb') as gz_file:
            gz_file.write(log_content.encode('utf-8'))
        
        compressed_data = buffer.getvalue()
        compressed_size_mb = len(compressed_data) / (1024 * 1024)
        
        filename = f"lava_{script_type.lower()}_logs_{timestamp}.txt.gz"
        
        from flask import Response
        response = Response(
            compressed_data,
            mimetype='application/gzip',
            headers={
                'Content-Disposition': f'attachment; filename={filename}',
                'X-Original-Size': f'{estimated_size_mb:.1f}MB',
                'X-Compressed-Size': f'{compressed_size_mb:.1f}MB'
            }
        )
    else:
        # Fichier texte normal
        filename = f"lava_{script_type.lower()}_logs_{timestamp}.txt"
        
        from flask import Response
        response = Response(
            log_content,
            mimetype='text/plain',
            headers={
                'Content-Disposition': f'attachment; filename={filename}',
                'X-File-Size': f'{estimated_size_mb:.1f}MB'
            }
        )
    
    return response

@app.route('/download-selected/<execution_id>', methods=['POST'])
def download_selected_files(execution_id):
    """Telecharger les fichiers selectionnes en ZIP"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    if 'result_files' not in execution:
        flash('Aucun résultat disponible', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    # Récupérer les fichiers sélectionnés
    selected_files = request.form.getlist('selected_files')
    if not selected_files:
        flash('Aucun fichier sélectionné', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    # Créer un fichier ZIP
    import zipfile
    import io
    
    script_type = execution.get('script_type', 'UNKNOWN')
    timestamp = execution.get('created_time', datetime.now()).strftime('%Y%m%d_%H%M%S')
    zip_filename = f"lava_{script_type.lower()}_results_{timestamp}.zip"
    
    # Créer le ZIP en mémoire
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        files_added = 0
        for file_path in execution['result_files']:
            filename = os.path.basename(file_path)
            if filename in selected_files and os.path.exists(file_path):
                zip_file.write(file_path, filename)
                files_added += 1
    
    if files_added == 0:
        flash('Aucun fichier valide trouvé', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    zip_buffer.seek(0)
    zip_size_mb = len(zip_buffer.getvalue()) / (1024 * 1024)
    
    from flask import Response
    response = Response(
        zip_buffer.getvalue(),
        mimetype='application/zip',
        headers={
            'Content-Disposition': f'attachment; filename={zip_filename}',
            'X-Files-Count': str(files_added),
            'X-Zip-Size': f'{zip_size_mb:.1f}MB'
        }
    )
    
    return response

@app.route('/download-all/<execution_id>')
def download_all_files(execution_id):
    """Telecharger tous les fichiers en ZIP"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    if 'result_files' not in execution or not execution['result_files']:
        flash('Aucun résultat disponible', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    # Créer un fichier ZIP avec tous les fichiers
    import zipfile
    import io
    
    script_type = execution.get('script_type', 'UNKNOWN')
    timestamp = execution.get('created_time', datetime.now()).strftime('%Y%m%d_%H%M%S')
    zip_filename = f"lava_{script_type.lower()}_all_results_{timestamp}.zip"
    
    # Créer le ZIP en mémoire
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        files_added = 0
        for file_path in execution['result_files']:
            if os.path.exists(file_path):
                filename = os.path.basename(file_path)
                zip_file.write(file_path, filename)
                files_added += 1
    
    zip_buffer.seek(0)
    zip_size_mb = len(zip_buffer.getvalue()) / (1024 * 1024)
    
    from flask import Response
    response = Response(
        zip_buffer.getvalue(),
        mimetype='application/zip',
        headers={
            'Content-Disposition': f'attachment; filename={zip_filename}',
            'X-Files-Count': str(files_added),
            'X-Zip-Size': f'{zip_size_mb:.1f}MB'
        }
    )
    
    return response

@app.route('/download/<execution_id>/<filename>')
def download_result(execution_id, filename):
    """Telecharger un fichier de resultat"""
    check_execution_ownership(execution_id)
    execution = running_executions[execution_id]
    if 'result_files' not in execution:
        flash('Aucun résultat disponible', 'error')
        return redirect(url_for('monitor', execution_id=execution_id))
    
    # Chercher le fichier demandé
    for result_file in execution['result_files']:
        if os.path.basename(result_file) == filename:
            return send_file(result_file, as_attachment=True)
    
    flash('Fichier non trouvé', 'error')
    return redirect(url_for('monitor', execution_id=execution_id))

@app.route('/executions')
def list_executions():
    """Liste toutes les executions"""
    current_user = session.get('user_id')
    filtered_execs = {k: v for k, v in running_executions.items() if v.get('owner_id') == current_user or v.get('user_id') == current_user}
    return render_template('executions.html', executions=filtered_execs)

@app.route('/stop/<execution_id>', methods=['POST'])
def stop_execution(execution_id):
    """Arreter une execution"""
    check_execution_ownership(execution_id)
    if execution_id in running_executions:
        execution = running_executions[execution_id]
        if 'process' in execution and execution['process'].poll() is None:
            execution['process'].terminate()
            execution['status'] = 'stopped'
        flash('Execution arretee', 'warning')
    
    return redirect(url_for('monitor', execution_id=execution_id))



# Routes conservatives pour détection fermeture
active_sessions = {}
shutdown_scheduled = False

@app.route('/heartbeat-simple', methods=['POST'])
def heartbeat_simple():
    """Heartbeat simple sans effets de bord"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    session_id = session.get('session_id', str(time.time()))
    session['session_id'] = session_id
    active_sessions[session_id] = time.time()
    
    global shutdown_scheduled
    if shutdown_scheduled:
        print("✅ Activité détectée - Annulation arrêt programmé")
        shutdown_scheduled = False
    
    return {'status': 'alive'}

@app.route('/really-closing', methods=['POST'])
def really_closing():
    """Signal de fermeture très conservatif"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    print("🔍 Signal de fermeture potentiel reçu...")
    
    def conservative_shutdown():
        global shutdown_scheduled
        shutdown_scheduled = True
        
        # Attendre 2 minutes avant de vraiment fermer
        time.sleep(120)
        
        # Vérifier s'il y a eu de l'activité
        if not shutdown_scheduled:
            print("✅ Arrêt annulé par activité")
            return
            
        # Vérifier les sessions actives
        current_time = time.time()
        for session_id, last_seen in active_sessions.items():
            if current_time - last_seen < 180:  # Activité dans les 3 dernières minutes
                print("✅ Session active trouvée - Pas d'arrêt")
                shutdown_scheduled = False
                return
        
        print("👋 Aucune activité confirmée - Arrêt de LAVA")
        os._exit(0)
    
    threading.Thread(target=conservative_shutdown, daemon=True).start()
    return {'status': 'checking'}

@app.route('/long-inactivity', methods=['POST'])
def long_inactivity():
    """Inactivité très prolongée détectée"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    print("💤 Inactivité prolongée (10+ minutes)")
    
    def inactivity_shutdown():
        time.sleep(300)  # Attendre encore 5 minutes
        
        # Vérifier une dernière fois
        current_time = time.time()
        for session_id, last_seen in active_sessions.items():
            if current_time - last_seen < 300:
                print("✅ Activité récente trouvée")
                return
        
        print("😴 Inactivité confirmée - Arrêt de LAVA")
        os._exit(0)
    
    threading.Thread(target=inactivity_shutdown, daemon=True).start()
    return {'status': 'monitoring'}



# Variables globales pour le tracking
browser_sessions = {}
shutdown_timer = None

@app.route('/stay-alive', methods=['POST'])
def stay_alive():
    """Signal de vie du navigateur"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    session_id = session.get('session_id', str(time.time()))
    session['session_id'] = session_id
    session['last_alive'] = datetime.now()
    
    browser_sessions[session_id] = datetime.now()
    
    # Annuler l'arrêt programmé s'il y en a un
    global shutdown_timer
    if shutdown_timer:
        shutdown_timer.cancel()
        shutdown_timer = None
    
    return {'status': 'alive', 'session': session_id}

@app.route('/maybe-closing', methods=['POST'])
def maybe_closing():
    """Signal possible de fermeture"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    print("🤔 Signal possible de fermeture reçu...")
    
    def check_if_really_closed():
        time.sleep(15)  # Attendre 15 secondes
        
        # Vérifier s'il y a eu de l'activité récente
        session_id = session.get('session_id')
        if session_id and session_id in browser_sessions:
            last_alive = browser_sessions[session_id]
            time_since = (datetime.now() - last_alive).total_seconds()
            
            if time_since > 30:  # Pas d'activité depuis 30 secondes
                print("🔒 Fermeture confirmée - Arrêt de l'interface")
                os._exit(0)
            else:
                print("✅ Activité détectée - Pas de fermeture")
    
    global shutdown_timer
    if shutdown_timer:
        shutdown_timer.cancel()
    
    shutdown_timer = threading.Timer(0, check_if_really_closed)
    shutdown_timer.start()
    
    return {'status': 'checking'}

@app.route('/page-hidden', methods=['POST'])
def page_hidden():
    """Page cachée depuis longtemps"""
    if not os.environ.get('ALLOW_SELF_SHUTDOWN', 'False').lower() in ('true', '1', 't'):
        return {'status': 'disabled'}, 404
    print("👁️  Page cachée détectée...")
    
    def delayed_check():
        time.sleep(60)  # Attendre 1 minute
        
        session_id = session.get('session_id')
        if session_id and session_id in browser_sessions:
            last_alive = browser_sessions[session_id]
            time_since = (datetime.now() - last_alive).total_seconds()
            
            if time_since > 120:  # Pas d'activité depuis 2 minutes
                print("💤 Inactivité prolongée - Arrêt de l'interface")
                os._exit(0)
    
    threading.Thread(target=delayed_check, daemon=True).start()
    return {'status': 'monitoring'}

def background_data_cleanup():
    """Purge synchronisée des fichiers et des exécutions terminées depuis longtemps"""
    import shutil
    while True:
        time_module.sleep(3600)  # Vérification toutes les heures
        try:
            retention_hours = float(os.environ.get('DATA_RETENTION_HOURS', 48))
            retention_seconds = retention_hours * 3600
            now = datetime.now()

            to_remove = []
            for exec_id, exec_data in list(running_executions.items()):
                if exec_data.get('status') in ['starting', 'running']:
                    continue
                end_time_str = exec_data.get('end_time')
                if end_time_str:
                    try:
                        end_time = datetime.strptime(end_time_str, '%Y-%m-%d %H:%M:%S')
                        if (now - end_time).total_seconds() > retention_seconds:
                            to_remove.append(exec_id)
                    except Exception:
                        pass

            for exec_id in to_remove:
                exec_data = running_executions.pop(exec_id, None)
                if exec_data:
                    upload_path = exec_data.get('input_file') or exec_data.get('target_file')
                    if upload_path and os.path.exists(upload_path):
                        try: os.remove(upload_path)
                        except Exception: pass
                    output_dir = exec_data.get('output_dir')
                    if output_dir and os.path.exists(output_dir):
                        try: shutil.rmtree(output_dir, ignore_errors=True)
                        except Exception: pass
        except Exception as e:
            print(f"Erreur purge automatique: {e}")

threading.Thread(target=background_data_cleanup, daemon=True).start()

if __name__ == '__main__':
    print("🧬 LAVA-DNA Flask Interface")
    print("=" * 40)
    print("✅ Interface web stable pour LAVA")
    print("🔗 Accès: http://localhost:5001")
    print("=" * 40)
    
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() in ('true', '1', 't')
    listen_host = os.environ.get('FLASK_HOST', '127.0.0.1')
    app.run(debug=debug_mode, host=listen_host, port=int(os.environ.get('FLASK_PORT', 5001)))
