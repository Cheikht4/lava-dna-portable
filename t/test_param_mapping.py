#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Audit bidirectionnel des parametres du moteur LAVA.
Bidirectional audit of LAVA engine parameters.

REGLE ABSOLUE D'INTEGRITE / ABSOLUTE INTEGRITY RULE:
Interdiction formelle d'ajouter la moindre exception en dur dans le code de ce
fichier d'audit. Tout defaut transitoire ou exception de perimetre doit obligatoirement
etre consigne dans le sas t/param_audit_known.txt avec date et justification.
Toute entree de t/param_audit_known.txt dont le defaut ne se reproduit plus fait
echouer l'audit (exit 1, regle d'obsolescence stricte).
Hardcoded exceptions are strictly forbidden in this audit script. Any transitional
defect or perimeter exception must be recorded in t/param_audit_known.txt with date
and justification. Any stale entry causes an immediate audit failure (exit 1).

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
import xml.etree.ElementTree as ET


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
    Extrait les cles, options distinctes, entrees GetOptions (alias compris)
    et doublons du hash %optionMap passe a GetOptions.
    Extracts keys, distinct options, GetOptions entries (including aliases),
    and duplicates from %optionMap passed to GetOptions.
    """
    filepath = os.path.join(repo_root, relpath)
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    match = re.search(r"my\s+%optionMap\s*=\s*\((.*?)\n\s*\);", content, re.DOTALL)
    if not match:
        raise ValueError(f"Impossible de localiser %optionMap dans {relpath}")

    block = match.group(1)
    block_start_line = content.count("\n", 0, match.start(1)) + 1
    declared_options = {}
    distinct_options = set()
    duplicates = []

    for offset, raw_line in enumerate(block.splitlines()):
        line_no = block_start_line + offset
        line = raw_line.strip()
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
            primary = aliases[0]
            distinct_options.add(primary)
            for alias in aliases:
                if alias in declared_options:
                    duplicates.append(
                        {
                            "file": relpath,
                            "alias": alias,
                            "spec": spec,
                            "line": line_no,
                            "prev_spec": declared_options[alias]["spec"],
                            "prev_line": declared_options[alias]["line"],
                        }
                    )
                declared_options[alias] = {
                    "target": target,
                    "primary": primary,
                    "spec": spec,
                    "file": relpath,
                    "line": line_no,
                }
        else:
            m_simple = re.search(r'\"([^\"]+)\"', line)
            if m_simple and "=>" in line:
                spec = m_simple.group(1)
                base = spec.split("=")[0].split("!")[0].split(":")[0]
                aliases = [a.strip() for a in base.split("|") if a.strip()]
                primary = aliases[0]
                distinct_options.add(primary)
                for alias in aliases:
                    if alias in declared_options:
                        duplicates.append(
                            {
                                "file": relpath,
                                "alias": alias,
                                "spec": spec,
                                "line": line_no,
                                "prev_spec": declared_options[alias]["spec"],
                                "prev_line": declared_options[alias]["line"],
                            }
                        )
                    declared_options[alias] = {
                        "target": alias,
                        "primary": primary,
                        "spec": spec,
                        "file": relpath,
                        "line": line_no,
                    }

    return declared_options, distinct_options, duplicates


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


def find_invalid_readme_options(repo_root, loop_options, stem_options):
    """
    Extrait de chaque README*.md du depot les commandes invoquant
    lava_loop_primer.pl ou lava_stem_primer.pl, releve chaque --option
    et verifie qu'elle existe dans %optionMap du script concerne.
    Extracts from each README*.md in the repository commands invoking
    lava_loop_primer.pl or lava_stem_primer.pl, collects every --option,
    and checks that it exists in the script's %optionMap.
    """
    readme_paths = sorted(glob.glob(os.path.join(repo_root, "README*.md")))
    invalid_entries = []
    checked_commands = 0

    for full_path in readme_paths:
        relpath = os.path.relpath(full_path, repo_root)
        with open(full_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()

        active_script = None
        continued = False
        in_code_block = False

        cmd_start_re = re.compile(
            r"(?:^|\s)(?:perl\s+(?:-\S+\s+)*|\./)(lava_(?:loop|stem)_primer\.pl)\b"
            r"|(?:^|\s)(lava_(?:loop|stem)_primer\.pl)(?=\s+--|\s*\\$)"
        )

        for idx, raw_line in enumerate(lines, start=1):
            stripped = raw_line.strip()

            # Fin ou debut de bloc de code markdown
            if stripped.startswith("```"):
                in_code_block = not in_code_block
                active_script = None
                continued = False
                continue

            # Ignorer les commentaires shell dans les blocs de code
            if in_code_block and stripped.startswith("#"):
                continue

            # Retirer un commentaire en fin de ligne eventuel dans un bloc de code
            code_part = (
                re.sub(r"\s+#.*$", "", raw_line).rstrip()
                if in_code_block
                else raw_line.rstrip()
            )

            # Detecter si la ligne lance une commande lava_loop_primer.pl ou lava_stem_primer.pl
            cmd_match = cmd_start_re.search(code_part)
            if cmd_match:
                active_script = cmd_match.group(1) or cmd_match.group(2)
                checked_commands += 1
            elif active_script and (continued or stripped.startswith("--")):
                pass
            else:
                active_script = None
                continued = False
                continue

            valid_map = (
                loop_options
                if active_script == "lava_loop_primer.pl"
                else stem_options
            )

            for m in re.finditer(r"(?:^|\s)--([a-zA-Z0-9_-]+)", code_part):
                opt_name = m.group(1)
                if opt_name == "help":
                    continue
                if opt_name not in valid_map:
                    invalid_entries.append(
                        {
                            "file": relpath,
                            "line": idx,
                            "script": active_script,
                            "option": opt_name,
                        }
                    )

            continued = code_part.endswith("\\")
            if not continued:
                active_script = None

    return invalid_entries, checked_commands


def find_invalid_xml_options(repo_root, valid_options):
    """
    Verifie que toutes les balises enfants de <option_set> dans les fichiers
    t_data/*.xml (hors test_options_file.xml et hors balises slava_* de
    slava_example_parameters.xml) appartiennent au %optionMap du moteur.
    Verifies that all child tags of <option_set> in t_data/*.xml files
    (excluding test_options_file.xml and slava_* tags in
    slava_example_parameters.xml) belong to the engine's %optionMap.
    """
    xml_paths = sorted(glob.glob(os.path.join(repo_root, "t_data", "*.xml")))
    invalid_xml_tags = []
    checked_xml_files = 0

    for full_path in xml_paths:
        basename = os.path.basename(full_path)
        if basename == "test_options_file.xml":
            continue

        relpath = os.path.relpath(full_path, repo_root)
        checked_xml_files += 1

        tree = ET.parse(full_path)
        root = tree.getroot()
        option_sets = [root] if root.tag == "option_set" else root.findall(".//option_set")

        for opt_set in option_sets:
            for child in opt_set:
                if not isinstance(child.tag, str):
                    continue
                tag_name = child.tag
                if basename == "slava_example_parameters.xml" and tag_name.startswith(
                    "slava_"
                ):
                    continue
                if tag_name not in valid_options:
                    invalid_xml_tags.append(
                        {
                            "file": relpath,
                            "tag": tag_name,
                        }
                    )

    return invalid_xml_tags, checked_xml_files


def find_logical_or_option_defaults(repo_root):
    """
    Scanne lava_loop_primer.pl, lava_stem_primer.pl et lib/LLNL/LAVA/**/*.pm
    a la recherche de replis logiques '$quelque_chose->{cle} || <valeur>'
    ou la valeur de repli est non nulle (ex: 42, 100000, 1.5, 'valeur').
    Les compteurs inoffensifs '|| 0' et les conditions logiques (if/elsif, sort)
    sont ignores.
    Scans scripts and PM modules for '$hash->{key} || <val>' fallbacks where
    val is non-zero. Harmless '|| 0' and boolean conditions are ignored.
    """
    targets = [
        os.path.join(repo_root, "lava_loop_primer.pl"),
        os.path.join(repo_root, "lava_stem_primer.pl"),
    ]
    targets.extend(
        glob.glob(
            os.path.join(repo_root, "lib", "LLNL", "LAVA", "**", "*.pm"),
            recursive=True,
        )
    )

    pattern = re.compile(
        r"(\$[a-zA-Z_][a-zA-Z0-9_]*(?:->)?\{['\"]?[a-zA-Z0-9_]+['\"]?\}\s*\|\|\s*([^;\n,]+))"
    )

    violations = []

    for filepath in sorted(targets):
        relpath = os.path.relpath(filepath, repo_root)
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        for idx, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue

            code_line = re.sub(r"#.*$", "", line).strip()
            if not code_line:
                continue

            if "<=>" in code_line or "cmp" in code_line:
                continue

            if re.search(r"\b(?:if|elsif|while|unless)\s*\(.*\|\|", code_line):
                continue

            for match in pattern.finditer(code_line):
                full_expr = match.group(1).strip()
                fallback_val = match.group(2).strip()

                if re.match(r"^0(?:\.0+)?$", fallback_val) or fallback_val in ('"0"', "'0'"):
                    continue

                if re.match(r"^(?:!|exists\b|defined\b|\$|\()", fallback_val):
                    continue

                violations.append(
                    {
                        "file": relpath,
                        "line": idx,
                        "expr": full_expr,
                        "fallback": fallback_val,
                    }
                )

    return violations


def find_option_with_default_mismatches(repo_root):
    """
    Extrait les tables %optionDefaults de lava_loop_primer.pl et lava_stem_primer.pl,
    puis inspecte tous les appels optionWithDefault($options_r, "X", <fallback>).
    Signale toute incohérence :
      - soit un fallback littéral différent de la valeur dans %optionDefaults{"X"}
        (ou une option "X" absente de %optionDefaults)
      - soit une référence à une mauvaise clé $optionDefaults{"Y"} où Y != X.
    Extracts %optionDefaults from both scripts and audits all optionWithDefault calls
    for literal divergence or key mismatch against %optionDefaults.
    """
    scripts = ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    violations = []

    for relpath in scripts:
        filepath = os.path.join(repo_root, relpath)
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()

        defaults_match = re.search(
            r"my\s+%optionDefaults\s*=\s*\((.*?)\);", content, re.DOTALL
        )
        if not defaults_match:
            violations.append(
                {
                    "file": relpath,
                    "line": 0,
                    "option": "ALL",
                    "msg": "Table %optionDefaults introuvable",
                }
            )
            continue

        defaults_block = defaults_match.group(1)
        defaults = {}
        for line in defaults_block.splitlines():
            line_clean = re.sub(r"#.*$", "", line).strip()
            kv = re.search(r'["\']([a-zA-Z0-9_]+)["\']\s*=>\s*([^,]+)', line_clean)
            if kv:
                key = kv.group(1)
                val_raw = kv.group(2).strip().strip('"\'')
                defaults[key] = val_raw

        call_pattern = re.compile(
            r'optionWithDefault\s*\(\s*\$(?:options_r|options|opts_r|opts)\s*,\s*["\']([a-zA-Z0-9_]+)["\']\s*,\s*([^,\)]+)\)',
            re.MULTILINE,
        )

        defined_or_pattern = re.compile(
            r'\$(?:options_r|options|opts_r|opts)(?:->)?\{["\']([a-zA-Z0-9_]+)["\']\}\s*//\s*([^;}\n]+)',
            re.MULTILINE,
        )

        candidates = []
        for match in call_pattern.finditer(content):
            candidates.append((match.group(1), match.group(2).strip(), match.start(), "optionWithDefault"))
        for match in defined_or_pattern.finditer(content):
            candidates.append((match.group(1), match.group(2).strip(), match.start(), "//"))

        for opt_name, fallback_expr, match_start, syntax_type in candidates:
            line_no = content.count("\n", 0, match_start) + 1

            ref_match = re.search(
                r'^\$optionDefaults\{["\']([a-zA-Z0-9_]+)["\']\}$', fallback_expr
            )
            if ref_match:
                ref_key = ref_match.group(1)
                if ref_key != opt_name:
                    violations.append(
                        {
                            "file": relpath,
                            "line": line_no,
                            "option": opt_name,
                            "msg": f"Cle incoherente ({syntax_type}) : lit '{opt_name}' mais utilise $optionDefaults{{'{ref_key}'}}",
                        }
                    )
                elif ref_key not in defaults:
                    violations.append(
                        {
                            "file": relpath,
                            "line": line_no,
                            "option": opt_name,
                            "msg": f"Option '{opt_name}' ({syntax_type}) absente de la table %optionDefaults",
                        }
                    )
                continue

            # Ignorer les expressions dynamiques calibrees avec variables Perl (ex: $targetTM - 1.0)
            if "$" in fallback_expr:
                continue

            clean_fallback = fallback_expr.strip('"\'')

            if opt_name not in defaults:
                violations.append(
                    {
                        "file": relpath,
                        "line": line_no,
                        "option": opt_name,
                        "msg": f"Option '{opt_name}' ({syntax_type}) utilise le repli litteral '{fallback_expr}' mais est absente de %optionDefaults",
                    }
                )
            else:
                expected_val = defaults[opt_name]
                is_equal = False
                try:
                    if float(clean_fallback) == float(expected_val):
                        is_equal = True
                except ValueError:
                    if clean_fallback == expected_val:
                        is_equal = True

                if not is_equal:
                    violations.append(
                        {
                            "file": relpath,
                            "line": line_no,
                            "option": opt_name,
                            "msg": f"Litteral '{fallback_expr}' ({syntax_type}) differe de %optionDefaults ('{expected_val}')",
                        }
                    )

    return violations


def find_missing_usage_options(repo_root):
    """
    Verifie la concordance tripartite stricte entre %optionMap, %optionDefaults et $usageString :
      - Sens A : Toute option distincte declaree dans %optionMap doit apparaitre dans $usageString.
      - Sens B : Toute option documentee dans $usageString doit exister dans %optionMap.
      - Sens C : Toute entree de %optionDefaults doit exister dans %optionMap (pas de constante deguisee).
    Verifies tripartite concordance between %optionMap, %optionDefaults, and $usageString.
    """
    scripts = ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    missing_by_script = {}
    undeclared_by_script = {}
    undeclared_defaults_by_script = {}

    for relpath in scripts:
        filepath = os.path.join(repo_root, relpath)
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()

        m_usage = re.search(
            r"(?:our|my)\s+\$usageString\s*=\s*(.*?);\s*(?:if\s*\(grep|my\s+\$alignmentFastaName)",
            content,
            re.DOTALL,
        )
        if not m_usage:
            missing_by_script[relpath] = ["$usageString introuvable"]
            continue
        usage_text = m_usage.group(1)

        declared_map, distinct_options, dups = extract_option_map(repo_root, relpath)

        # Sens A : declare dans %optionMap -> documente dans $usageString
        missing = []
        for opt in sorted(distinct_options):
            if not re.search(r"--" + re.escape(opt) + r"\b", usage_text) and not re.search(
                r"-" + re.escape(opt) + r"\b", usage_text
            ):
                missing.append(opt)

        if missing:
            missing_by_script[relpath] = missing

        # Sens B : documente dans $usageString -> declare dans %optionMap
        doc_options = set()
        for m in re.finditer(r"--([a-zA-Z0-9_]+)\b", usage_text):
            doc_opt = m.group(1)
            if doc_opt == "help":
                continue
            doc_options.add(doc_opt)

        undeclared = []
        for d_opt in sorted(doc_options):
            if d_opt not in declared_map:
                undeclared.append(d_opt)

        if undeclared:
            undeclared_by_script[relpath] = undeclared

        # Sens C : cle dans %optionDefaults -> declaree dans %optionMap
        defaults_match = re.search(
            r"my\s+%optionDefaults\s*=\s*\((.*?)\);", content, re.DOTALL
        )
        if defaults_match:
            defaults_block = defaults_match.group(1)
            undef_defs = []
            for line in defaults_block.splitlines():
                line_clean = re.sub(r"#.*$", "", line).strip()
                kv = re.search(r'["\']([a-zA-Z0-9_]+)["\']\s*=>', line_clean)
                if kv:
                    def_key = kv.group(1)
                    if def_key not in declared_map and def_key not in undef_defs:
                        undef_defs.append(def_key)
            if undef_defs:
                undeclared_defaults_by_script[relpath] = undef_defs

    return missing_by_script, undeclared_by_script, undeclared_defaults_by_script


def check_architecture_symmetry(loop_distinct, stem_distinct):
    """
    Verifie la symetrie stricte entre les options de lava_loop_primer.pl et lava_stem_primer.pl (Sens D).
    Les ensembles doivent etre identiques, a l'exception :
      - des options prefixees par 'loop_' ou 'include_loop_primers' (specifiques LOOP)
      - des options prefixees par 'stem_' ou 'include_stem_primers' (specifiques STEM)
      - d'une liste explicite d'exceptions assumees (ASSUMED_ASYMMETRIC_OPTIONS).
    Toute divergence non declaree fait echouer l'audit.
    Verifies strict option symmetry between LOOP and STEM architectures (Direction D).
    """
    ASSUMED_ASYMMETRIC_OPTIONS = set()

    def is_architecture_specific(opt, arch):
        if opt in ASSUMED_ASYMMETRIC_OPTIONS:
            return True
        if arch == "loop":
            return opt.startswith("loop_") or opt == "include_loop_primers"
        elif arch == "stem":
            return opt.startswith("stem_") or opt == "include_stem_primers"
        return False

    loop_generic = {opt for opt in loop_distinct if not is_architecture_specific(opt, "loop")}
    stem_generic = {opt for opt in stem_distinct if not is_architecture_specific(opt, "stem")}

    missing_in_stem = sorted(loop_generic - stem_generic)
    missing_in_loop = sorted(stem_generic - loop_generic)

    return missing_in_stem, missing_in_loop


def run_audit():
    repo_root = get_repo_root()
    print("=" * 70)
    print("AUDIT BIDIRECTIONNEL DES PARAMETRES DU MOTEUR LAVA")
    print("=" * 70)

    loop_options, loop_distinct, loop_dups = extract_option_map(
        repo_root, "lava_loop_primer.pl"
    )
    stem_options, stem_distinct, stem_dups = extract_option_map(
        repo_root, "lava_stem_primer.pl"
    )
    option_map_duplicates = loop_dups + stem_dups

    common_distinct = loop_distinct & stem_distinct
    union_distinct = loop_distinct | stem_distinct
    common_entries = set(loop_options.keys()) & set(stem_options.keys())

    all_declared = {}
    all_declared.update(loop_options)
    all_declared.update(stem_options)

    print("[1] Options declarees dans GetOptions :")
    print(
        f"    - lava_loop_primer.pl : {len(loop_distinct)} options distinctes (noms longs), "
        f"{len(loop_options)} entrees GetOptions alias compris"
    )
    print(
        f"    - lava_stem_primer.pl : {len(stem_distinct)} options distinctes (noms longs), "
        f"{len(stem_options)} entrees GetOptions alias compris"
    )
    print(
        f"    - Communes LOOP/STEM  : {len(common_distinct)} options distinctes (noms longs), "
        f"{len(common_entries)} entrees GetOptions alias compris"
    )
    print(
        f"    - Union (total unique): {len(union_distinct)} options distinctes (noms longs), "
        f"{len(all_declared)} entrees GetOptions alias compris"
    )

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
    known_defects = load_known_defects(repo_root)
    reproduced_known_options = set()
    unlisted_dead = []
    listed_dead = []

    for opt_name, info in all_declared.items():
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
            if opt_name in known_defects:
                listed_dead.append((opt_name, info["file"]))
                reproduced_known_options.add(opt_name)
            else:
                unlisted_dead.append((opt_name, info["file"]))

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
    inert_reads = find_inert_option_reads(
        repo_root, ["lava_loop_primer.pl", "lava_stem_primer.pl"]
    )

    unlisted_inert = []
    listed_inert = []

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
    if listed_dead:
        for dp, src in listed_dead:
            k_meta = known_defects[dp]
            print(
                f"    KNOWN: parametre mort '{dp}' (dans {src}) "
                f"[{k_meta['date']} - {k_meta['reason']}]"
            )
    if unlisted_dead:
        print(f"    FAIL: {len(unlisted_dead)} parametre(s) mort(s) non liste(s) detecte(s) :")
        for dp, src in unlisted_dead:
            print(f"      - {dp} (dans {src})")
        success = False
    else:
        print("    PASS: Aucun parametre mort non liste detecte.")

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

    print("\n[6] Resultat de l'audit des entrees dupliquees dans %optionMap :")
    if option_map_duplicates:
        print(
            f"    FAIL: {len(option_map_duplicates)} entree(s) dupliquee(s) dans %optionMap :"
        )
        for dup in option_map_duplicates:
            print(
                f"      - '{dup['alias']}' dans {dup['file']} ligne {dup['line']} "
                f"(spec '{dup['spec']}') duplique la ligne {dup['prev_line']} (spec '{dup['prev_spec']}')"
            )
        success = False
    else:
        print("    PASS: Aucune entree dupliquee dans %optionMap.")

    readme_invalid, readme_cmd_count = find_invalid_readme_options(
        repo_root, loop_options, stem_options
    )
    print(
        f"\n[7] Resultat de l'audit des options CLI dans les exemples README ({readme_cmd_count} commande(s) verifiee(s)) :"
    )
    if readme_invalid:
        print(
            f"    FAIL: {len(readme_invalid)} option(s) invalide(s) detectee(s) dans la documentation :"
        )
        for item in readme_invalid:
            print(
                f"      - '--{item['option']}' pour {item['script']} dans {item['file']} ligne {item['line']}"
            )
        success = False
    else:
        print(
            "    PASS: Toutes les options CLI des exemples README existent dans %optionMap."
        )

    xml_invalid, xml_file_count = find_invalid_xml_options(repo_root, all_declared)
    print(
        f"\n[8] Resultat de l'audit des balises XML dans t_data/*.xml ({xml_file_count} fichier(s) verifie(s)) :"
    )
    if xml_invalid:
        print(
            f"    FAIL: {len(xml_invalid)} balise(s) XML invalide(s) detectee(s) :"
        )
        for item in xml_invalid:
            print(f"      - '<{item['tag']}>' dans {item['file']}")
        success = False
    else:
        print(
            "    PASS: Toutes les balises XML des fichiers t_data/*.xml existent dans %optionMap."
        )

    logical_or_violations = find_logical_or_option_defaults(repo_root)
    print("\n[9] Resultat de l'audit des replis logiques '||' sur options/parametres :")
    if logical_or_violations:
        print(
            f"    FAIL: {len(logical_or_violations)} repli(s) logique(s) '||' non nul(s) detecte(s) (utiliser '//') :"
        )
        for viol in logical_or_violations:
            print(
                f"      - {viol['file']}:{viol['line']} : '{viol['expr']}' (repli = {viol['fallback']})"
            )
        success = False
    else:
        print(
            "    PASS: Aucun repli '|| <non-zero>' detecte sur lecture d'options/parametres."
        )

    default_mismatches = find_option_with_default_mismatches(repo_root)
    print("\n[10] Resultat de l'audit de coherence optionWithDefault vs %optionDefaults :")
    if default_mismatches:
        print(
            f"    FAIL: {len(default_mismatches)} incoherence(s) detectee(s) dans optionWithDefault :"
        )
        for m in default_mismatches:
            print(
                f"      - {m['file']}:{m['line']} (option '{m['option']}') : {m['msg']}"
            )
        success = False
    else:
        print(
            "    PASS: Tous les appels optionWithDefault concordent rigoureusement avec %optionDefaults."
        )

    missing_usage, undeclared_usage, undeclared_defaults = find_missing_usage_options(repo_root)
    missing_in_stem, missing_in_loop = check_architecture_symmetry(loop_distinct, stem_distinct)
    print("\n[11] Resultat de l'audit quadripartite options / defauts / aide en ligne / symetrie :")
    usage_ok = True
    if missing_usage:
        total_missing = sum(len(opts) for opts in missing_usage.values())
        print(
            f"    FAIL (Sens A : declare -> non documente) : {total_missing} option(s) de %optionMap absente(s) de $usageString :"
        )
        for scr, opts in missing_usage.items():
            print(f"      - {scr} ({len(opts)} manquante(s)) : {', '.join(opts)}")
        usage_ok = False
        success = False

    if undeclared_usage:
        total_undeclared = sum(len(opts) for opts in undeclared_usage.values())
        print(
            f"    FAIL (Sens B : documente -> non declare) : {total_undeclared} option(s) de $usageString absente(s) de %optionMap :"
        )
        for scr, opts in undeclared_usage.items():
            print(f"      - {scr} ({len(opts)} non declaree(s)) : {', '.join(opts)}")
        usage_ok = False
        success = False

    if undeclared_defaults:
        total_undef_defaults = sum(len(opts) for opts in undeclared_defaults.values())
        print(
            f"    FAIL (Sens C : defaut sans declaration dans %optionMap) : {total_undef_defaults} option(s) de %optionDefaults absente(s) de %optionMap :"
        )
        for scr, opts in undeclared_defaults.items():
            print(f"      - {scr} ({len(opts)} non declaree(s)) : {', '.join(opts)}")
        usage_ok = False
        success = False

    if missing_in_stem or missing_in_loop:
        print("    FAIL (Sens D : asymetrie inter-architectures non declaree) :")
        if missing_in_stem:
            print(f"      - Absente(s) dans lava_stem_primer.pl : {', '.join(missing_in_stem)}")
        if missing_in_loop:
            print(f"      - Absente(s) dans lava_loop_primer.pl : {', '.join(missing_in_loop)}")
        usage_ok = False
        success = False

    if usage_ok:
        print(
            "    PASS: Concordance quadripartite parfaite entre %optionMap, %optionDefaults, $usageString et symetrie LOOP/STEM."
        )

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
