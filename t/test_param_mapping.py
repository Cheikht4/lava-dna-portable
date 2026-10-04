#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Audit bidirectionnel des parametres du moteur LAVA.
Bidirectional audit of LAVA engine parameters.

Detecte les deux moities du probleme :
1. Parametre declare mais jamais lu (parametre mort).
2. Parametre lu mais jamais declare dans GetOptions (parametre orphelin).
3. Verification de l'absence totale des anciens parametres supprimes.
"""

import os
import re
import sys
import glob

def get_code_files():
    """
    Retourne la liste des fichiers de code executables du moteur et des librairies.
    Returns the list of executable code files of the engine and libraries.
    """
    files = ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    files.extend(glob.glob("lib/LLNL/LAVA/**/*.pm", recursive=True))
    files.extend(glob.glob("lib/Bio/**/*.pm", recursive=True))
    return sorted(files)

def extract_option_map(filepath):
    """
    Extrait les cles et variables cibles du hash %optionMap passe a GetOptions.
    Extracts keys and target variables from %optionMap passed to GetOptions.
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    match = re.search(r"my\s+%optionMap\s*=\s*\((.*?)\n\s*\);", content, re.DOTALL)
    if not match:
        raise ValueError(f"Impossible de localiser %optionMap dans {filepath}")

    block = match.group(1)
    declared_options = {}

    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Format Getopt::Long: "cle=type" => \$options{"cible"} ou "cle|alias=type" => ...
        m_target = re.search(r'\"([^\"]+)\"\s*=>\s*.*?(?:\$options|\$opt)\{[\'\"]?([a-zA-Z0-9_]+)[\'\"]?\}', line)
        if m_target:
            spec = m_target.group(1)
            target = m_target.group(2)
            base = spec.split("=")[0].split("!")[0].split(":")[0]
            aliases = [a.strip() for a in base.split("|") if a.strip()]
            for alias in aliases:
                declared_options[alias] = {
                    "target": target,
                    "spec": spec,
                    "file": filepath
                }
        else:
            m_simple = re.search(r'\"([^\"]+)\"', line)
            if m_simple and "=>" in line:
                spec = m_simple.group(1)
                base = spec.split("=")[0].split("!")[0].split(":")[0]
                aliases = [a.strip() for a in base.split("|") if a.strip()]
                for alias in aliases:
                    declared_options[alias] = {
                        "target": alias,
                        "spec": spec,
                        "file": filepath
                    }

    return declared_options

def run_audit():
    print("=" * 70)
    print("AUDIT BIDIRECTIONNEL DES PARAMETRES DU MOTEUR LAVA")
    print("=" * 70)

    loop_options = extract_option_map("lava_loop_primer.pl")
    stem_options = extract_option_map("lava_stem_primer.pl")
    
    all_declared = {}
    all_declared.update(loop_options)
    all_declared.update(stem_options)
    
    print(f"[1] Options declarees dans GetOptions :")
    print(f"    - lava_loop_primer.pl : {len(loop_options)} options")
    print(f"    - lava_stem_primer.pl : {len(stem_options)} options")
    print(f"    - Total unique        : {len(all_declared)} options")

    code_files = get_code_files()
    
    # Contenus sources hors declaration %optionMap
    file_contents = {}
    for cf in code_files:
        with open(cf, "r", encoding="utf-8", errors="replace") as f:
            c = f.read()
        if cf in ["lava_loop_primer.pl", "lava_stem_primer.pl"]:
            c = re.sub(r"my\s+%optionMap\s*=\s*\(.*?\n\s*\);", "", c, flags=re.DOTALL)
        file_contents[cf] = c

    # -------------------------------------------------------------------------
    # PARTIE 1 : Detection des parametres morts (declares mais jamais lus)
    # PART 1: Detection of dead parameters (declared but never read)
    # -------------------------------------------------------------------------
    dead_parameters = []
    
    # Parametres de transition connus et documentes :
    # - output_dir : reserve pour la gestion des repertoires dans lava_flask_app
    # - assembly_batch_size : conserve pour compatibilite apres refonte map-reduce
    KNOWN_TRANSITIONAL_PARAMS = {"output_dir", "assembly_batch_size"}

    for opt_name, info in all_declared.items():
        if opt_name in KNOWN_TRANSITIONAL_PARAMS:
            continue
            
        target = info["target"]
        found = False

        for cf, content in file_contents.items():
            patterns = [
                r'options(?:_r|_ref)?->\{[\'\"]?' + re.escape(opt_name) + r'[\'\"]?\}',
                r'options(?:_r|_ref)?->\{[\'\"]?' + re.escape(target) + r'[\'\"]?\}',
                r'options\{[\'\"]?' + re.escape(opt_name) + r'[\'\"]?\}',
                r'options\{[\'\"]?' + re.escape(target) + r'[\'\"]?\}',
                r'optionWithDefault\([^,]+,\s*[\'\"]' + re.escape(opt_name) + r'[\'\"]',
                r'optionWithDefault\([^,]+,\s*[\'\"]' + re.escape(target) + r'[\'\"]',
                r'optionRequired\([^,]+,\s*[\'\"]' + re.escape(opt_name) + r'[\'\"]',
                r'optionRequired\([^,]+,\s*[\'\"]' + re.escape(target) + r'[\'\"]',
                r'optionDefaults\{[\'\"]?' + re.escape(opt_name) + r'[\'\"]?\}',
                r'optionDefaults\{[\'\"]?' + re.escape(target) + r'[\'\"]?\}',
                r'getOption\([\'\"]' + re.escape(opt_name) + r'[\'\"]\)',
                r'getOption\([\'\"]' + re.escape(target) + r'[\'\"]\)',
                # Prise en charge specifique du chargement d options XML (Options.pm)
                r'optionTitle\s*=\s*[\'\"]' + re.escape(target) + r'[\'\"]'
            ]
            for pat in patterns:
                if re.search(pat, content):
                    found = True
                    break
            if found:
                break

        if not found:
            dead_parameters.append((opt_name, info["file"]))

    # -------------------------------------------------------------------------
    # PARTIE 2 : Detection des parametres orphelins (lus mais jamais declares)
    # PART 2: Detection of orphan parameters (read but never declared)
    # -------------------------------------------------------------------------
    engine_read_patterns = [
        r'options_r->\{[\'\"]([a-zA-Z0-9_]+)[\'\"]\}',
        r'options_ref->\{[\'\"]([a-zA-Z0-9_]+)[\'\"]\}',
        r'optionWithDefault\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*[\'\"]([a-zA-Z0-9_]+)[\'\"]',
        r'optionRequired\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*[\'\"]([a-zA-Z0-9_]+)[\'\"]',
    ]

    referenced_keys = set()
    for cf, content in file_contents.items():
        for pat in engine_read_patterns:
            for match in re.finditer(pat, content):
                referenced_keys.add(match.group(1))

    # Fallbacks internes ou synonymes algorithmiques documentes
    INTERNAL_ENGINE_FALLBACKS = {
        "min_gc", "max_gc", "salt_conc", "mg_conc",
        "primer3_tm_formula", "primer3_salt_corrections",
        "fixed_primer_margin"
    }

    orphan_parameters = []
    for k in sorted(referenced_keys):
        if k not in all_declared and k not in INTERNAL_ENGINE_FALLBACKS:
            orphan_parameters.append(k)

    # -------------------------------------------------------------------------
    # PARTIE 3 : Verification de l'absence des parametres supprimes
    # PART 3: Verification of absence of deprecated/removed parameters
    # -------------------------------------------------------------------------
    forbidden_parameters = [
        "min_signatures_for_success",
        "max_dist_outer_middle",
        "max_dist_middle_inner",
        "penalty_slope",
        "penalty_plateau"
    ]
    forbidden_found = {}
    for cf, content in file_contents.items():
        for f_param in forbidden_parameters:
            if re.search(r'\b' + re.escape(f_param) + r'\b', content):
                forbidden_found.setdefault(f_param, []).append(cf)

    # -------------------------------------------------------------------------
    # BILAN ET VERDICT DU TEST
    # SUMMARY AND TEST VERDICT
    # -------------------------------------------------------------------------
    success = True

    print("\n[2] Resultat de l'audit des parametres morts :")
    if dead_parameters:
        print(f"    FAIL: {len(dead_parameters)} parametre(s) declare(s) mais jamais lu(s) :")
        for dp, src in dead_parameters:
            print(f"      - {dp} (dans {src})")
        success = False
    else:
        print("    PASS: Aucun parametre mort detecte.")

    print("\n[3] Resultat de l'audit des parametres orphelins :")
    if orphan_parameters:
        print(f"    FAIL: {len(orphan_parameters)} parametre(s) lu(s) mais jamais declare(s) :")
        for op in orphan_parameters:
            print(f"      - {op}")
        success = False
    else:
        print("    PASS: Aucun parametre orphelin detecte.")

    print("\n[4] Resultat de l'audit des parametres supprimes (interdits) :")
    if forbidden_found:
        print("    FAIL: Des parametres interdits/supprimes ont ete trouves dans le code :")
        for fp, files in forbidden_found.items():
            print(f"      - {fp} dans : {', '.join(files)}")
        success = False
    else:
        print("    PASS: Aucun parametre interdit detecte dans le code executable.")

    print("=" * 70)
    if success:
        print("VERDICT FINAL : PASS (Tous les audits sont valides)")
        print("=" * 70)
        return 0
    else:
        print("VERDICT FINAL : FAIL (Anomalies detectees)")
        print("=" * 70)
        return 1

if __name__ == "__main__":
    sys.exit(run_audit())
