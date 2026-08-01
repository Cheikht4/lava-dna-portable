#!/bin/bash
set -eo pipefail

DRY_RUN=1
if [[ " $* " =~ " --apply " ]]; then
    DRY_RUN=0
fi

# Define paths
SOURCE_DIR=$(pwd)
# Ensure we are in lava-dna-portable
if [[ ! -f "lava_loop_primer.pl" || ! -d "lib/LLNL/LAVA" ]]; then
    echo "ERREUR: Ce script doit etre lance depuis la racine de lava-dna-portable."
    exit 1
fi

# Try to read target repos from arguments or use default
TARGETS=()
for arg in "$@"; do
    if [[ "$arg" != "--apply" && "$arg" != "--dry-run" ]]; then
        TARGETS+=("$arg")
    fi
done

if [ ${#TARGETS[@]} -eq 0 ]; then
    TARGETS=("../../lava-virus-public" "../../lava-interface-public")
fi

WHITELIST=(
    "lava_loop_primer.pl"
    "lava_stem_primer.pl"
    "lava_flask_app.py"
    "templates"
    "lib"
    "t/baseline"
    "t/canary_regression.t"
    "t/fixtures"
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
        # Recursively process files
        find "$src_item" -type f | while read -r src_file; do
            # Ignore certain generated or system files
            if [[ "$src_file" == *".DS_Store" ]] || [[ "$src_file" == *"_amplified_noms.txt" ]]; then
                continue
            fi
            # Ignore amplified fasta EXCEPT in baseline directory
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

echo "============================================================"
if [ "$DRY_RUN" -eq 1 ]; then
    echo "MODE: --dry-run (aucune modification ne sera apportee)"
else
    echo "MODE: --apply (copie des fichiers en cours...)"
fi
echo "============================================================"

for target in "${TARGETS[@]}"; do
    if [ ! -d "$target" ]; then
        echo "ATTENTION: Le depot cible '$target' n'existe pas."
        continue
    fi
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

if [ "$DRY_RUN" -eq 1 ]; then
    echo ""
    echo "Relancez avec --apply pour appliquer les modifications et executer les tests."
    exit 0
fi

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
# Compare .pl, app and lib/LLNL/LAVA
declare -a FILES_TO_CHECK=("lava_loop_primer.pl" "lava_stem_primer.pl" "lava_flask_app.py")
while IFS= read -r f; do
    FILES_TO_CHECK+=("$f")
done < <(find lib/LLNL/LAVA -type f -name "*.pm")

COHERENCE_OK=1
for file in "${FILES_TO_CHECK[@]}"; do
    ref_md5=$(compute_md5 "$file")
    for target in "${TARGETS[@]}"; do
        if [ ! -d "$target" ]; then continue; fi
        tgt_md5=$(compute_md5 "$target/$file")
        if [[ "$ref_md5" != "$tgt_md5" ]]; then
            echo "DIVERGENCE DETECTEE: $file differe entre principal et $target"
            COHERENCE_OK=0
        fi
    done
done

if [ "$COHERENCE_OK" -eq 1 ]; then
    echo "Tous les fichiers critiques (.pl, modules LAVA) sont strictement identiques."
else
    echo "ATTENTION: Des divergences ont ete detectees apres la synchronisation !"
    exit 1
fi

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
