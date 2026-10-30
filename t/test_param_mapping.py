#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Audit bidirectionnel des parametres du moteur LAVA.
Bidirectional audit of LAVA engine parameters.

Detecte quatre classes d'anomalies :
1. Parametre declare mais jamais lu (parametre mort).
2. Parametre lu mais jamais declare dans GetOptions (parametre orphelin).
3. Parametre lu via optionWithDefault dans une variable $VAR qui n'est jamais
   reutilisee dans le fichier (lecture inerte), avec gestion stricte de la
   liste des defauts connus (t/param_audit_known.txt).
4. Verification de l'absence totale des anciens parametres supprimes.
"""

import os
import re
import sys
import glob
import subprocess


def get_repo_root():
    """
    Determine la racine du depot Git via git rev-parse --show-toplevel.
    Determines the Git repository root via git rev-parse --show-toplevel.
    """
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        if out and os.path.isdir(out):
            return out
    except Exception:
        pass
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def get_code_files(repo_root):
    """
    Retourne la liste des fichiers de code executables du moteur et des librairies.
    Returns the list of executable code files of the engine and libraries.
    """
    files = ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    for p in glob.glob(os.path.join(repo_root, "lib/LLNL/LAVA/**/*.pm"), recursive=True):
        files.append(os.path.relpath(p, repo_root))
    for p in glob.glob(os.path.join(repo_root, "lib/Bio/**/*.pm"), recursive=True):
        files.append(os.path.relpath(p, repo_root))
    return sorted(files)


def load_known_defects(repo_root):
    """
    Charge la liste des defauts connus depuis t/param_audit_known.txt.
    Loads the list of known defects from t/param_audit_known.txt.
    Format : nom_option | YYYY-MM-DD | justification
    """
    known_path = os.path.join(repo_root, "t", "param_audit_known.txt")
    known = {}
    if not os.path.isfile(known_path):
        return known

    with open(known_path, "r", encoding="utf-8", errors="replace") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("|")]
            opt_name = parts[0]
            date_str = parts[1] if len(parts) > 1 else ""
            reason = parts[2] if len(parts) > 2 else ""
            known[opt_name] = {
                "date": date_str,
                "reason": reason,
                "line": line_no,
            }
    return known


def extract_option_map(repo_root, relpath):
    """
    Extrait les cles et variables cibles du hash %optionMap passe a GetOptions.
    Extracts keys and target variables from %optionMap passed to GetOptions.
    """
    filepath = os.path.join(repo_root, relpath)
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    match = re.search(r"my\s+%optionMap\s*=\s*\((.*?)\n\s*\);", content, re.DOTALL)
    if not match:
        raise ValueError(f"Impossible de localiser %optionMap dans {relpath}")

    block = match.group(1)
    declared_options = {}

    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Format Getopt::Long: "cle=type" => \$options{"cible"} ou "cle|alias=type" => ...
        m_target = re.search(
            r'\"([^\"]+)\"\s*=>\s*.*?(?:\$options|\$opt)\{[\'\"]?([a-zA-Z0-9_]+)[\'\"]?\}',
            line,
        )
        if m_target:
            spec = m_target.group(1)
            target = m_target.group(2)
            base = spec.split("=")[0].split("!")[0].split(":")[0]
            aliases = [a.strip() for a in base.split("|") if a.strip()]
            for alias in aliases:
                declared_options[alias] = {
                    "target": target,
                    "spec": spec,
                    "file": relpath,
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
                        "file": relpath,
                    }

    return declared_options


def find_inert_option_reads(repo_root, scripts):
    """
    Detecte les affectations de la forme :
      my $VAR = optionWithDefault($options_r, "nom", ...)
    pour lesquelles $VAR n'apparait qu'une seule fois (son affectation) dans le fichier.
    Detects assignments of the form:
      my $VAR = optionWithDefault($options_r, "name", ...)
    where $VAR appears only once (its assignment) in the file.
    """
    assign_pattern = re.compile(
        r"my\s+\$([a-zA-Z_][a-zA-Z0-9_]*)\s*=\s*(?:int\s*\(\s*)?optionWithDefault\s*\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*[\"\']([a-zA-Z0-9_]+)[\"\']",
        re.MULTILINE,
    )

    inert_reads = []
    for relpath in scripts:
        filepath = os.path.join(repo_root, relpath)
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()

        # Exclure les lignes purement commentees pour compter les occurrences actives du nom de variable
        # Exclude pure comment lines when counting active occurrences of the variable name
        active_lines = [
            re.sub(r"#.*$", "", line)
            for line in content.splitlines()
            if not line.lstrip().startswith("#")
        ]
        active_text = "\n".join(active_lines)

        for match in assign_pattern.finditer(content):
            line_start = content.rfind("\n", 0, match.start()) + 1
            if content[line_start:match.start()].lstrip().startswith("#"):
                continue

            var_name = match.group(1)
            opt_name = match.group(2)
            line_no = content.count("\n", 0, match.start()) + 1

            # Compter toutes les occurrences du nom de variable ($VAR ou ${VAR})
            # Count all occurrences of the variable name ($VAR or ${VAR})
            occ_pattern = re.compile(
                r"\$(?:\{" + re.escape(var_name) + r"\}|" + re.escape(var_name) + r"\b)"
            )
            occ_count = len(occ_pattern.findall(active_text))

            if occ_count <= 1:
                inert_reads.append(
                    {
                        "option": opt_name,
                        "var": var_name,
                        "file": relpath,
                        "line": line_no,
                        "count": occ_count,
                    }
                )

    return inert_reads


def run_audit():
    repo_root = get_repo_root()
    print("=" * 70)
    print("AUDIT BIDIRECTIONNEL DES PARAMETRES DU MOTEUR LAVA")
    print("=" * 70)

    loop_options = extract_option_map(repo_root, "lava_loop_primer.pl")
    stem_options = extract_option_map(repo_root, "lava_stem_primer.pl")

    all_declared = {}
    all_declared.update(loop_options)
    all_declared.update(stem_options)

    print(f"[1] Options declarees dans GetOptions :")
    print(f"    - lava_loop_primer.pl : {len(loop_options)} options")
    print(f"    - lava_stem_primer.pl : {len(stem_options)} options")
    print(f"    - Total unique        : {len(all_declared)} options")

    code_files = get_code_files(repo_root)

    # Contenus sources hors declaration %optionMap
    file_contents = {}
    for cf in code_files:
        full_cf = os.path.join(repo_root, cf)
        with open(full_cf, "r", encoding="utf-8", errors="replace") as f:
            c = f.read()
        if cf in ["lava_loop_primer.pl", "lava_stem_primer.pl"]:
            c = re.sub(r"my\s+%optionMap\s*=\s*\(.*?\n\s*\);", "", c, flags=re.DOTALL)
        file_contents[cf] = c

    # -------------------------------------------------------------------------
    # PARTIE 1 : Detection des parametres morts (declares mais jamais lus)
    # PART 1: Detection of dead parameters (declared but never read)
    # -------------------------------------------------------------------------
    dead_parameters = []

    KNOWN_TRANSITIONAL_PARAMS = {"output_dir", "assembly_batch_size"}

    for opt_name, info in all_declared.items():
        if opt_name in KNOWN_TRANSITIONAL_PARAMS:
            continue

        target = info["target"]
        found = False

        for cf, content in file_contents.items():
            patterns = [
                r"options(?:_r|_ref)?->\{[\'\"]?" + re.escape(opt_name) + r"[\'\"]?\}",
                r"options(?:_r|_ref)?->\{[\'\"]?" + re.escape(target) + r"[\'\"]?\}",
                r"options\{[\'\"]?" + re.escape(opt_name) + r"[\'\"]?\}",
                r"options\{[\'\"]?" + re.escape(target) + r"[\'\"]?\}",
                r"optionWithDefault\([^,]+,\s*[\'\"]" + re.escape(opt_name) + r"[\'\"]",
                r"optionWithDefault\([^,]+,\s*[\'\"]" + re.escape(target) + r"[\'\"]",
                r"optionRequired\([^,]+,\s*[\'\"]" + re.escape(opt_name) + r"[\'\"]",
                r"optionRequired\([^,]+,\s*[\'\"]" + re.escape(target) + r"[\'\"]",
                r"optionDefaults\{[\'\"]?" + re.escape(opt_name) + r"[\'\"]?\}",
                r"optionDefaults\{[\'\"]?" + re.escape(target) + r"[\'\"]?\}",
                r"getOption\([\'\"]" + re.escape(opt_name) + r"[\'\"]\)",
                r"getOption\([\'\"]" + re.escape(target) + r"[\'\"]\)",
                r"optionTitle\s*=\s*[\'\"]" + re.escape(target) + r"[\'\"]",
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
        r"options_r->\{[\'\"]([a-zA-Z0-9_]+)[\'\"]\}",
        r"options_ref->\{[\'\"]([a-zA-Z0-9_]+)[\'\"]\}",
        r"optionWithDefault\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*[\'\"]([a-zA-Z0-9_]+)[\'\"]",
        r"optionRequired\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*[\'\"]([a-zA-Z0-9_]+)[\'\"]",
    ]

    referenced_keys = set()
    for cf, content in file_contents.items():
        for pat in engine_read_patterns:
            for match in re.finditer(pat, content):
                referenced_keys.add(match.group(1))

    INTERNAL_ENGINE_FALLBACKS = {
        "min_gc",
        "max_gc",
        "salt_conc",
        "mg_conc",
        "primer3_tm_formula",
        "primer3_salt_corrections",
        "fixed_primer_margin",
    }

    orphan_parameters = []
    for k in sorted(referenced_keys):
        if k not in all_declared and k not in INTERNAL_ENGINE_FALLBACKS:
            orphan_parameters.append(k)

    # -------------------------------------------------------------------------
    # PARTIE 3 : Detection des lectures inertes (my $VAR = optionWithDefault
    #            dont $VAR n'apparait qu'une seule fois dans le fichier)
    # PART 3: Detection of inert reads ($VAR assigned via optionWithDefault
    #         and never reused in the file)
    # -------------------------------------------------------------------------
    known_defects = load_known_defects(repo_root)
    inert_reads = find_inert_option_reads(
        repo_root, ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    )

    unlisted_inert = []
    listed_inert = []
    reproduced_known_options = set()

    for item in inert_reads:
        opt = item["option"]
        if opt in known_defects:
            listed_inert.append(item)
            reproduced_known_options.add(opt)
        else:
            unlisted_inert.append(item)

    stale_known = []
    for k_opt, k_info in known_defects.items():
        if k_opt not in reproduced_known_options:
            stale_known.append((k_opt, k_info))

    # -------------------------------------------------------------------------
    # PARTIE 4 : Verification de l'absence des parametres supprimes
    # PART 4: Verification of absence of deprecated/removed parameters
    # -------------------------------------------------------------------------
    forbidden_parameters = [
        "min_signatures_for_success",
        "max_dist_outer_middle",
        "max_dist_middle_inner",
        "penalty_slope",
        "penalty_plateau",
        "outer_pair_target_length",
        "middle_pair_target_length",
        "inner_pair_target_length",
    ]
    forbidden_found = {}
    for cf, content in file_contents.items():
        for f_param in forbidden_parameters:
            if re.search(r"\b" + re.escape(f_param) + r"\b", content):
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

    print("\n[4] Resultat de l'audit des lectures inertes (variables assignees mais jamais utilisees) :")
    if listed_inert:
        for item in listed_inert:
            k_meta = known_defects[item["option"]]
            print(
                f"    KNOWN: option '{item['option']}' (${item['var']}) dans {item['file']}:{item['line']} "
                f"[{k_meta['date']} - {k_meta['reason']}]"
            )
    if unlisted_inert:
        print(f"    FAIL: {len(unlisted_inert)} lecture(s) inerte(s) non listee(s) detectee(s) :")
        for item in unlisted_inert:
            print(
                f"      - option '{item['option']}' (${item['var']}) dans {item['file']} ligne {item['line']} "
                f"(occurrences de ${item['var']} = {item['count']})"
            )
        success = False
    if stale_known:
        print(f"    FAIL: {len(stale_known)} entree(s) obsolete(s) dans t/param_audit_known.txt :")
        for k_opt, k_info in stale_known:
            print(
                f"      - '{k_opt}' (ligne {k_info['line']}) : entree obsolete, le defaut est corrige, retire la ligne"
            )
        success = False
    if not unlisted_inert and not stale_known:
        print("    PASS: Aucune lecture inerte non listee ni entree obsolete dans t/param_audit_known.txt.")

    print("\n[5] Resultat de l'audit des parametres supprimes (interdits) :")
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
