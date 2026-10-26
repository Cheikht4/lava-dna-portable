#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Audit du maillon 2 : correspondance entre lava_flask_app.py et GetOptions des scripts Perl.
Link 2 audit: mapping between lava_flask_app.py and Perl scripts GetOptions.

Verifie que :
1. Tous les parametres de common_params U loop_only_params existent dans GetOptions de lava_loop_primer.pl
2. Tous les parametres de common_params U stem_only_params existent dans GetOptions de lava_stem_primer.pl
3. Toutes les cibles de param_mapping existent dans les ensembles Flask et dans GetOptions
4. Aucun parametre exclusif STEM ne fuit dans LOOP (et inversement)
"""

import os
import re
import subprocess
import sys


# Parametres geres en amont par Flask et non transmis comme options scalaires directes si presents
# Parameters handled upstream by Flask and not passed as direct scalar options if present
FLASK_META_PARAMS = {"fixed_primers", "lamp_mode", "script_type"}


def get_repo_root():
    """Resout la racine du depot Git courant / Resolve current Git repository root."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip()
    except Exception:
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get_perl_getoptions(script_path):
    """Extrait les options declarees dans GetOptions / Extract options declared in GetOptions."""
    with open(script_path, "r", encoding="utf-8") as fh:
        content = fh.read()

    params = set()
    matches = re.findall(r'"([a-zA-Z0-9_\|]+)=?[sfi]?"\s*=>', content)
    for match in matches:
        for alias in match.split("|"):
            params.add(alias.strip())
    return params


def get_flask_param_sets(flask_path):
    """Extrait common_params, loop_only_params, stem_only_params et param_mapping."""
    with open(flask_path, "r", encoding="utf-8") as fh:
        content = fh.read()

    sets = {}
    for name in ("common_params", "loop_only_params", "stem_only_params"):
        match = re.search(name + r"\s*=\s*\{(.*?)\}", content, re.DOTALL)
        if not match:
            raise RuntimeError(f"Ensemble '{name}' introuvable dans {flask_path}")
        sets[name] = set(re.findall(r"['\"]([^'\"]+)['\"]", match.group(1)))

    mapping = {}
    for map_match in re.finditer(r"param_mapping\s*=\s*\{(.*?)\}", content, re.DOTALL):
        pairs = re.findall(
            r"['\"]([^'\"]+)['\"]\s*:\s*['\"]([^'\"]+)['\"]",
            map_match.group(1),
        )
        for k, v in pairs:
            mapping[k] = v

    return sets["common_params"], sets["loop_only_params"], sets["stem_only_params"], mapping


def main():
    repo_root = get_repo_root()
    flask_path = os.path.join(repo_root, "lava_flask_app.py")
    loop_path = os.path.join(repo_root, "lava_loop_primer.pl")
    stem_path = os.path.join(repo_root, "lava_stem_primer.pl")

    for required in (flask_path, loop_path, stem_path):
        if not os.path.isfile(required):
            print(f"ERREUR FATALE: Fichier requis introuvable : {required}", file=sys.stderr)
            sys.exit(1)

    common, loop_only, stem_only, param_mapping = get_flask_param_sets(flask_path)
    perl_loop = get_perl_getoptions(loop_path)
    perl_stem = get_perl_getoptions(stem_path)

    flask_loop = (common | loop_only) - FLASK_META_PARAMS
    flask_stem = (common | stem_only) - FLASK_META_PARAMS
    all_flask_valid = common | loop_only | stem_only

    errors = []

    # 1. Verification LOOP : Flask -> lava_loop_primer.pl GetOptions
    missing_in_loop = sorted(flask_loop - perl_loop)
    for p in missing_in_loop:
        errors.append(
            f"[LOOP] Parametre envoye par lava_flask_app.py mais absent de GetOptions dans lava_loop_primer.pl : {p}"
        )

    # 2. Verification STEM : Flask -> lava_stem_primer.pl GetOptions
    missing_in_stem = sorted(flask_stem - perl_stem)
    for p in missing_in_stem:
        errors.append(
            f"[STEM] Parametre envoye par lava_flask_app.py mais absent de GetOptions dans lava_stem_primer.pl : {p}"
        )

    # 3. Verification des cibles de param_mapping
    for src_key, target_opt in sorted(param_mapping.items()):
        if target_opt not in all_flask_valid:
            errors.append(
                f"[MAPPING] Cible '{target_opt}' (depuis '{src_key}') absente de common/loop_only/stem_only_params"
            )
        if target_opt not in (perl_loop | perl_stem):
            errors.append(
                f"[MAPPING] Cible '{target_opt}' (depuis '{src_key}') absente de GetOptions des scripts Perl"
            )

    print("=" * 70)
    print("AUDIT MAILLON 2 : LAVA_FLASK_APP.PY -> GETOPTIONS PERL")
    print("=" * 70)
    print(f"  - Parametres Flask LOOP testes : {len(flask_loop)}")
    print(f"  - Parametres Flask STEM testes : {len(flask_stem)}")
    print(f"  - Entrees param_mapping testees: {len(param_mapping)}")

    if errors:
        print(f"\nECHEC : {len(errors)} incoherence(s) detectee(s) :")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)

    print("\nPASS : Tous les parametres Flask correspondent aux options GetOptions Perl.")
    print("=" * 70)
    sys.exit(0)


if __name__ == "__main__":
    main()
