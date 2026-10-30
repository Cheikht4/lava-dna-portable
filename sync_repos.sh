#!/bin/bash
set -eo pipefail

DRY_RUN=1
if [[ " $* " =~ " --apply " ]]; then
    DRY_RUN=0
fi

# Define paths from repository root / Definir les chemins depuis la racine du depot
REPO_ROOT="$(git rev-parse --show-toplevel)"
LAVA_ROOT="$(dirname "$REPO_ROOT")"
if [ ! -d "$LAVA_ROOT/.agents" ] && [ -d "$(dirname "$LAVA_ROOT")/.agents" ]; then
    LAVA_ROOT="$(dirname "$LAVA_ROOT")"
fi
cd "$REPO_ROOT"

# Ensure we are in lava-dna-portable / Verifier qu'on est dans lava-dna-portable
if [[ ! -f "lava_loop_primer.pl" || ! -d "lib/LLNL/LAVA" ]]; then
    echo "ERREUR: Ce script doit etre lance depuis lava-dna-portable." >&2
    exit 1
fi

# Try to read target repos from arguments or use default / Lire les cibles ou utiliser les valeurs par defaut
TARGETS=()
for arg in "$@"; do
    if [[ "$arg" != "--apply" && "$arg" != "--dry-run" ]]; then
        TARGETS+=("$arg")
    fi
done

if [ ${#TARGETS[@]} -eq 0 ]; then
    TARGETS=("$LAVA_ROOT/lava-virus-public" "$LAVA_ROOT/lava-interface-public")
fi

# Determine expected origin URL for a target / Determiner l'URL d'origin attendue pour une cible
get_expected_origin() {
    local target="$1"
    local index="$2"

    if [ -f "$target/.engine-only" ]; then
        echo "https://github.com/Cheikht4/LAVA-Virus.git"
    elif [ -f "$target/lava_flask_app.py" ] && [ ! -f "$target/sync_repos.sh" ]; then
        echo "https://github.com/Cheikht4/LAVA_Virus-Interface.git"
    elif [ "$index" -eq 0 ]; then
        echo "https://github.com/Cheikht4/LAVA-Virus.git"
    else
        echo "https://github.com/Cheikht4/LAVA_Virus-Interface.git"
    fi
}

# Strict guard for each target repository / Garde stricte pour chaque depot cible
validate_targets() {
    local idx=0
    for target in "${TARGETS[@]}"; do
        local expected_origin
        expected_origin="$(get_expected_origin "$target" "$idx")"

        # 1. Directory exists / Le repertoire existe
        if [ ! -d "$target" ]; then
            echo "ERREUR FATALE: La cible '$target' est introuvable (condition non remplie : le repertoire n'existe pas)." >&2
            exit 1
        fi

        # 2. Valid git repository / C'est un depot git valide
        if [ ! -d "$target/.git" ] && ! git -C "$target" rev-parse --git-dir >/dev/null 2>&1; then
            echo "ERREUR FATALE: La cible '$target' est invalide (condition non remplie : ce n'est pas un depot git)." >&2
            exit 1
        fi

        # 3. Expected origin URL / L'origin correspond a l'URL attendue
        local actual_origin
        actual_origin="$(git -C "$target" remote get-url origin 2>/dev/null || echo "AUCUN_ORIGIN")"
        if [ "$actual_origin" != "$expected_origin" ]; then
            echo "ERREUR FATALE: La cible '$target' a un origin inattendu : '$actual_origin' (condition non remplie : origin attendu '$expected_origin')." >&2
            exit 1
        fi

        # 4. Clean working tree / L'arbre de travail est propre
        local dirty
        dirty="$(git -C "$target" status --porcelain 2>/dev/null)"
        if [ -n "$dirty" ]; then
            echo "ERREUR FATALE: La cible '$target' n'est pas propre (condition non remplie : modifications locales ou fichiers non suivis presents)." >&2
            exit 1
        fi

        echo "Cible validee : $target (origin: $actual_origin, arbre propre)"
        idx=$((idx + 1))
    done
}

WHITELIST=(
    "lava_loop_primer.pl"
    "lava_stem_primer.pl"
    "lava_flask_app.py"
    "templates"
    "lib"
    "t/baseline"
    "t/canary_regression.t"
    "t/fixtures"
    "t/test_form_names.py"
    "t/test_flask_param_mapping.py"
    "t/test_param_mapping.py"
    "t/param_audit_known.txt"
    ".github/workflows/canary.yml"
)

# Configuration du hash selon l'OS (Linux md5sum vs macOS md5)
if command -v md5sum >/dev/null 2>&1; then
    HASH_CMD="md5sum"
elif command -v md5 >/dev/null 2>&1; then
    HASH_CMD="md5"
else
    echo "ERREUR : ni md5sum ni md5 disponible." >&2
    exit 1
fi

compute_md5() {
    local file="$1"
    if [ -f "$file" ]; then
        if [ "$HASH_CMD" = "md5sum" ]; then
            md5sum "$file" 2>/dev/null | awk '{print $1}' || echo "ERROR"
        else
            md5 -q "$file" 2>/dev/null || echo "ERROR"
        fi
    else
        echo "MISSING_FILE"
    fi
}

# Auto-test de protection
echo "Verification du systeme de hachage..."
TEST_HASH1=$(compute_md5 "lava_loop_primer.pl")
TEST_HASH2=$(compute_md5 "t/canary_regression.t")
if [ "$TEST_HASH1" = "ERROR" ] || [ "$TEST_HASH2" = "ERROR" ] || [ "$TEST_HASH1" = "MISSING_FILE" ]; then
    echo "ERREUR CRITIQUE : Impossible de calculer le hash des fichiers sources." >&2
    exit 1
fi
if [ "$TEST_HASH1" = "$TEST_HASH2" ]; then
    echo "ERREUR CRITIQUE : Le systeme de hachage renvoie des valeurs identiques pour des fichiers differents. Annulation." >&2
    exit 1
fi
echo "Hachage OK."

sync_item() {
    local src_item="$1"
    local target_repo="$2"
    local dry_run="$3"

    # If it's a file
    if [ -f "$src_item" ]; then
        if [ -f "$target_repo/.engine-only" ]; then
            if [[ "$src_item" == *"lava_flask_app.py"* ]] || [[ "$src_item" == *"templates"* ]] || [[ "$src_item" == *"deployment"* ]] || [[ "$src_item" == *"test_form_names.py"* ]] || [[ "$src_item" == *"test_flask_param_mapping.py"* ]]; then
                return
            fi
        fi

        local target_item="$target_repo/$src_item"
        local src_md5=$(compute_md5 "$src_item")
        
        if [ "$src_md5" = "MISSING_FILE" ] || [ "$src_md5" = "ERROR" ]; then
            echo "ERREUR FATALE: Le fichier source $src_item est manquant ou illisible. Le script est interrompu." >&2
            exit 1
        fi
        
        local target_md5=$(compute_md5 "$target_item")

        if [[ "$src_md5" != "$target_md5" ]]; then
            echo "  - [MODIFIE] $src_item"
            if [ "$dry_run" -eq 0 ]; then
                mkdir -p "$(dirname "$target_item")"
                cp "$src_item" "$target_item"
            fi
        fi
    # If it's a directory
    elif [ -d "$src_item" ]; then
        if [ -f "$target_repo/.engine-only" ]; then
            if [[ "$src_item" == *"templates"* ]] || [[ "$src_item" == *"deployment"* ]]; then
                return
            fi
        fi

        # Recursively process files
        find "$src_item" -type f | while read -r src_file; do
            # Ignore certain generated or system files / Ignorer les fichiers systeme
            if [[ "$src_file" == *".DS_Store" ]]; then
                continue
            fi
            # Ignore amplified_noms.txt EXCEPT in baseline directory / Ignorer _amplified_noms.txt SAUF dans t/baseline/
            if [[ "$src_file" == *"_amplified_noms.txt" ]] && [[ ! "$src_file" =~ "t/baseline/" ]]; then
                continue
            fi
            # Ignore amplified fasta EXCEPT in baseline directory / Ignorer _amplified.fasta SAUF dans t/baseline/
            if [[ "$src_file" == *"_amplified.fasta" ]] && [[ ! "$src_file" =~ "t/baseline/" ]]; then
                continue
            fi

            local target_item="$target_repo/$src_file"
            local src_md5=$(compute_md5 "$src_file")
            
            if [ "$src_md5" = "MISSING_FILE" ] || [ "$src_md5" = "ERROR" ]; then
                echo "ERREUR FATALE: Le fichier source $src_file est manquant ou illisible. Le script est interrompu." >&2
                exit 1
            fi
            
            local target_md5=$(compute_md5 "$target_item")

            if [[ "$src_md5" != "$target_md5" ]]; then
                echo "  - [MODIFIE] $src_file"
                if [ "$dry_run" -eq 0 ]; then
                    mkdir -p "$(dirname "$target_item")"
                    cp "$src_file" "$target_item"
                fi
            fi
        done
    fi
}

# Strict guard on baseline file counts in git index / Garde stricte sur le compte des baselines dans l'index git
verify_baseline_counts() {
    local dry_run="$1"
    if [ "$dry_run" -eq 1 ]; then
        echo "Garde baselines : non evaluable en --dry-run (necessite copie et indexation git dans la cible)."
        return 0
    fi

    local expected_count obtained_count disk_count
    expected_count="$(git -C "$REPO_ROOT" ls-files t/baseline/ | wc -l | tr -d ' ')"

    for target in "${TARGETS[@]}"; do
        # Indexer t/baseline/ dans la cible afin que git applique son .gitignore
        # Stage t/baseline/ in target so git applies its .gitignore rules
        git -C "$target" add -A t/baseline/ 2>/dev/null || true

        obtained_count="$(git -C "$target" ls-files t/baseline/ | wc -l | tr -d ' ')"
        disk_count=0
        if [ -d "$target/t/baseline" ]; then
            disk_count="$(find "$target/t/baseline" -type f ! -name ".DS_Store" | wc -l | tr -d ' ')"
        fi

        if [ "$obtained_count" -ne "$expected_count" ]; then
            git -C "$target" reset -q HEAD -- t/baseline/ 2>/dev/null || true
            echo "ERREUR FATALE: Garde baselines en echec pour '$target' (attendu : $expected_count, obtenu : $obtained_count, sur disque : $disk_count). Index git de la cible restaure (git reset)." >&2
            exit 1
        fi
        echo "Garde baselines OK : $target (attendu : $expected_count, obtenu : $obtained_count, sur disque : $disk_count)"
    done
}

# Full MD5 coherence check across all synchronized items / Verification complete de coherence MD5 sur tout le perimetre
verify_checksum_coherence() {
    declare -a FILES_TO_CHECK=(
        "lava_loop_primer.pl"
        "lava_stem_primer.pl"
        "lava_flask_app.py"
        "t/canary_regression.t"
        "t/test_form_names.py"
        "t/test_flask_param_mapping.py"
        "t/test_param_mapping.py"
        "t/param_audit_known.txt"
        ".github/workflows/canary.yml"
    )
    while IFS= read -r f; do
        FILES_TO_CHECK+=("$f")
    done < <(find lib -type f -name "*.pm" | sort)
    while IFS= read -r f; do
        FILES_TO_CHECK+=("$f")
    done < <(find t/baseline -type f ! -name ".DS_Store" | sort)
    while IFS= read -r f; do
        FILES_TO_CHECK+=("$f")
    done < <(find t/fixtures -type f ! -name ".DS_Store" | sort)
    while IFS= read -r f; do
        FILES_TO_CHECK+=("$f")
    done < <(find templates -type f ! -name ".DS_Store" | sort)

    local coherence_ok=1
    local ref_md5 tgt_md5
    for file in "${FILES_TO_CHECK[@]}"; do
        ref_md5="$(compute_md5 "$file")"
        for target in "${TARGETS[@]}"; do
            if [ ! -d "$target" ]; then continue; fi

            # Exclusions declaratives pour depot moteur pur (.engine-only) / Declarative UI exclusions via .engine-only
            if [ -f "$target/.engine-only" ]; then
                if [[ "$file" == "lava_flask_app.py" ]] || [[ "$file" == templates/* ]] || [[ "$file" == "t/test_form_names.py" ]] || [[ "$file" == "t/test_flask_param_mapping.py" ]]; then
                    continue
                fi
            fi

            tgt_md5="$(compute_md5 "$target/$file")"
            if [[ "$ref_md5" != "$tgt_md5" ]]; then
                echo "DIVERGENCE DETECTEE: $file differe entre principal ($ref_md5) et $target ($tgt_md5)" >&2
                coherence_ok=0
            fi
        done
    done

    if [ "$coherence_ok" -eq 1 ]; then
        echo "Tous les fichiers du perimetre (.pl, lib/**, t/baseline/**, tests, CI, templates) sont strictement identiques."
    else
        echo "ERREUR FATALE: Des divergences MD5 ont ete detectees apres la synchronisation !" >&2
        exit 1
    fi
}

echo "============================================================"
if [ "$DRY_RUN" -eq 1 ]; then
    echo "MODE: --dry-run (aucune modification ne sera apportee)"
else
    echo "MODE: --apply (copie des fichiers en cours...)"
fi
echo "============================================================"

# Verification stricte des cibles avant toute operation / Strict target validation before any operation
validate_targets
echo "------------------------------------------------------------"

for target in "${TARGETS[@]}"; do
    echo "Analyse des modifications pour : $target"
    
    # Capture output of sync to count and display later
    SYNC_OUT=$(
    for item in "${WHITELIST[@]}"; do
        if [ -e "$item" ]; then
            sync_item "$item" "$target" "$DRY_RUN"
        fi
    done
    )
    
    if [ -n "$SYNC_OUT" ]; then
        echo "$SYNC_OUT"
    else
        echo "  (Deja a jour)"
    fi
done

echo "------------------------------------------------------------"

if [ "$DRY_RUN" -eq 1 ]; then
    verify_baseline_counts 1
    echo ""
    echo "Relancez avec --apply pour appliquer les modifications et executer les tests."
    exit 0
fi

echo ""
echo "============================================================"
echo "GARDE SUR LE COMPTE DES BASELINES INDEXEES PAR GIT"
echo "============================================================"
verify_baseline_counts 0

echo ""
echo "============================================================"
echo "TESTS CANARY SUR LES DEPOTS CIBLES"
echo "============================================================"

for target in "${TARGETS[@]}"; do
    if [ ! -d "$target" ]; then
        continue
    fi
    echo "-> Lancement du canary sur $target"
    (
        cd "$target"
        export PERL5LIB="$PWD/lib"
        rm -rf t/canary_*_signatures_individuelles
        if prove t/canary_regression.t > /dev/null 2>&1; then
            echo "   [PASS] Canary 15/15 OK"
            exit 0
        else
            echo "   [FAIL] ERREUR LORS DU CANARY !"
            prove t/canary_regression.t # Run again to show output before exit
            exit 1
        fi
    )
    if [ $? -ne 0 ]; then
        echo "ERREUR CRITIQUE: Le depot $target a echoue au canary test."
        exit 1
    fi
done

echo ""
echo "============================================================"
echo "VERIFICATION DE COHERENCE (CHECKSUMS)"
echo "============================================================"
verify_checksum_coherence

echo ""
echo "============================================================"
echo "COMMANDES GIT SUGGEREES"
echo "============================================================"
for target in "${TARGETS[@]}"; do
    if [ ! -d "$target" ]; then continue; fi
    echo "Pour $target :"
    echo "  cd \"$target\""
    echo "  git status"
    echo "  git add ."
    echo "  git commit -m \"Sync: mise a jour du moteur scientifique et des baselines depuis lava-dna-portable\""
    echo "  git push origin main"
    echo ""
done
