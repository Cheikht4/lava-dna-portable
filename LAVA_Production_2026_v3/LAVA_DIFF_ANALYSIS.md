# 🧬 Analyse Exhaustive des Modifications LAVA-DNA
## Comparaison entre `pseudogene/lava-dna` (original) et le fork modifié (2026)

**Date de l'analyse** : 24 avril 2026  
**Dépôt original** : https://github.com/pseudogene/lava-dna.git  
**Auteur de l'analyse** : Analyse bioinformatique automatisée  

---

## 📊 Résumé Statistique Global

| Catégorie | Original | Fork modifié | Delta |
|-----------|----------|-------------|-------|
| **Fichiers Perl principaux** | 1 (`lava.pl`, 2818 lignes) | 2 (`lava_stem_primer.pl` 3934L, `lava_loop_primer.pl` 3452L) | +1 fichier, +4568 lignes |
| **Modules Perl (lib/)** | 14 fichiers | 20 fichiers (+6 nouveaux) | +6 modules |
| **Application Web (Python)** | 0 | 2 fichiers (1741 lignes) | +1741 lignes |
| **Templates HTML** | 0 | 5 fichiers (1542 lignes) | +1542 lignes |
| **Déploiement** | 0 | 4 fichiers (394 lignes) | +394 lignes |
| **Documentation** | 1 (`README`) | 6 fichiers | +5 documents |
| **Script utilitaire (`slava.pl`)** | 1 (9808 octets) | 0 (supprimé) | -1 fichier |

---

## 1. 🔴 FICHIERS SUPPRIMÉS

### 1.1 `lava.pl` → Remplacé par deux scripts spécialisés
- **Rôle original** : Script monolithique unique de design LAMP (2818 lignes)
- **Raison de la suppression** : Séparation architecturale en deux pipelines distincts (STEM et LOOP) pour permettre un design d'amorces adapté à chaque topologie LAMP
- **Impact biologique** : Permet de traiter indépendamment les amorces STEM (positionnées entre F1c et B1c) et les amorces LOOP (positionnées entre F1c/F2 et B1c/B2), avec des contraintes géométriques spécifiques à chaque architecture

### 1.2 `slava.pl` → Supprimé
- **Rôle original** : Script auxiliaire de lancement LAVA (9808 octets)
- **Raison** : Fonctionnalité absorbée par la nouvelle interface Flask et le script `launch_lava_smart_kill.py`

---

## 2. 🟢 FICHIERS AJOUTÉS — Scripts Perl Principaux

### 2.1 `lava_stem_primer.pl` [NOUVEAU] — 3934 lignes (153 Ko)
**Nature** : Version spécialisée dérivée de `lava.pl` pour l'architecture STEM

**Architecture LAMP STEM** :
```
5' ←── Signature LAMP (≤ 400 nt) ──→ 3'
F3 ──(≥1)── F2 ──(≥1)── F1──(≥1)── FL════BL ──(≥1)── B1 ──(≥1)── B2 ──(≥1)── B3
                         └── Inner+STEM Pair (≈50 nt) ──┘
             └──────── Middle Pair (≈180 nt) ────────┘
└──────────────────── Outer Pair (≈250 nt) ──────────────────────┘
```

**Modifications majeures par rapport à `lava.pl`** :

| Fonctionnalité | `lava.pl` (original) | `lava_stem_primer.pl` (fork) |
|---|---|---|
| Tolérance mismatches | Aucune (conservation 100%) | `getOligosWithMismatchTolerance()` — analyse IUPAC position par position |
| Gestion dégénérescence | Inexistante | Paramètres `max_total_degenerate_bases`, `max_consecutive_degenerate_bases`, `max_3prime_degenerate_bases` |
| Validation de signature | Aucune | `calculateSignatureIntersection()` — intersection de couverture de tous les primers |
| Export de résultats | Fichier texte simple | Fichiers par signature individuelle avec rapport de validation (VALID/REJECT) + FASTA amplifiés/exclus |
| Combinatoire | Aucune | `analyzeSignatureCombinations()` — recherche de la meilleure combinaison de signatures pour couverture maximale |
| Paramètres STEM | Inexistants | `stem_primer_target_length/min/max`, `stem_primer_target_tm/min/max` |
| Modules importés | Aucun module spécifique | `LLNL::LAVA::Core`, `LLNL::LAVA::Validator` |
| Filtrage entropie | Basique | Paramètre `entropy_threshold` configurable |

**Nouvelles fonctions ajoutées** :
- `getOligosWithMismatchTolerance()` — Pipeline de tolérance aux mismatches avec codes IUPAC
- `calculateSignatureIntersection()` — Calcul d'intersection de couverture entre tous les primers (6 ou 8)
- `analyzeSignatureCombinations()` — Analyse combinatoire pour maximiser la couverture génomique
- `generateCombinations()` — Génération récursive de combinaisons
- `createPerSignatureFiles()` — Export détaillé par signature avec statut de validation
- `createAmplificationFiles()` — Export FASTA des séquences amplifiées et exclues

### 2.2 `lava_loop_primer.pl` [NOUVEAU] — 3452 lignes (140 Ko)
**Nature** : Version spécialisée dérivée de `lava.pl` pour l'architecture LOOP

**Différences clés par rapport à `lava_stem_primer.pl`** :

| Aspect | STEM | LOOP |
|--------|------|------|
| Primers supplémentaires | FSTEM/BSTEM (entre F1c et B1c) | FLOOP/BLOOP (entre F1c/F2 et B1c/B2) |
| Paramètre de couverture | `include_stem_primers` | `include_loop_primers` avec seuil `min_signature_coverage` (défaut 70%) |
| Fonction de pénalité géométrique | `calculate_proportional_geometry` | `generateSigmoidPenalty` |
| Résolution chevauchement | Non | `max_overlap_percent`, `resolve_overlap_by` (penalty/coverage) |
| Validation signature | Retourne 2 valeurs | Retourne 3 valeurs (+ `validation_status`) |
| Tags stockés dans signature | Couverture basique | `signature_intersection_ids`, `signature_coverage_percent`, `signature_target_count`, `validation_status`, `primer_coverage_details` |

---

## 3. 🟢 FICHIERS AJOUTÉS — Nouveaux Modules Perl

### 3.1 `lib/LLNL/LAVA/Core.pm` [NOUVEAU] — 115 lignes
**Rôle** : Noyau de fonctions utilitaires partagées entre STEM et LOOP

**Fonctions exportées** :
- **`calculate_proportional_geometry($L)`** — Calcule les distances cibles inter-primers basées sur des ratios proportionnels (F3-F2: 12%, F2-F1: 18%, F1-B1: 40%)
- **`generateSigmoidPenalty($actual, $target, $plateau_ratio, $k_slope)`** — Fonction de pénalité sigmoïde remplaçant la parabole originale. Zone de confort (plateau 0) à ±25% de la cible, montée logistique douce au-delà
- **`generateDistancePenalties($maxDistance, $targetLength)`** — Génère un tableau de pénalités pour toutes les distances possibles
- **`countDegenerateBases($sequence)`** — Compte les bases non-standard (IUPAC dégénérées) dans une séquence

**Justification biologique** : La courbe sigmoïde (vs. parabole) modélise plus fidèlement le comportement thermodynamique réel de la Bst polymérase à 65°C : une zone de tolérance large suivie d'une pénalité croissante fluide, évitant les rejets abruptes de candidats potentiellement viables.

### 3.2 `lib/LLNL/LAVA/Validator.pm` [NOUVEAU] — 449 lignes
**Rôle** : Module centralisé de validation des amorces et signatures LAMP

**Fonctions exportées** :
- **`checkPrimerMismatchTolerance()`** — Algorithme en 4 phases :
  1. Extraction des régions cibles (gap-aware)
  2. Test d'orientation SENSE vs ANTISENSE avec reverse-complement automatique
  3. Analyse position par position avec génération de codes IUPAC sous contraintes (max total, max consécutifs, max en 3')
  4. Validation finale avec protection de la zone 3' et tolérance aux mismatches
- **`isIUPACCompatible($base, $iupac_code)`** — Vérifie la compatibilité IUPAC base-à-code
- **`rev_comp($seq)`** — Complément inversé avec support IUPAC complet
- **`generateIUPACCode($bases_ref)`** — Génère le code IUPAC pour un ensemble de bases
- **`getPrimerTargetedSequences()`** — Identifie les séquences ciblées par un primer donné
- **`validateCompleteSignatureSpacing()`** — Valide l'espacement et le non-chevauchement de tous les primers d'une signature

**Impact biologique** : Ce module permet de designer des amorces LAMP sur des virus hautement variables (ex: Dengue) en tolérant une diversité génomique contrôlée via les codes IUPAC, tout en protégeant l'extrémité 3' critique pour l'extension enzymatique.

### 3.3 `lib/Lava/Core.pm` [NOUVEAU] — 1001 lignes
**Rôle** : Noyau de refactoring orienté objet encapsulant l'ancien pipeline `lava.pl`

**Fonctions principales** :
- `run_lava_loop()` / `run_lava_stem()` — Points d'entrée programmatiques (vs. CLI)
- `buildReversePrimers()` — Construction des amorces reverse par complément inversé
- `analyzeAll()` — Analyse batch via PrimerAnalyzer
- `enumeratePairs()` — Énumération de paires Forward/Reverse compatibles
- `reducePairInfosByPenalty()` — Sélection des meilleures paires par score
- `reducePrimersByOverlap()` — Filtrage par chevauchement maximal
- `reduceSignaturesByOverlap()` — Filtrage de signatures LAMP complètes
- `validateF1B1Spacing()` — Validation de l'espacement F1/B1

### 3.4 `lib/Lava/Enumerator/StemConserved.pm` [NOUVEAU] — 81 lignes
**Rôle** : Énumérateur spécialisé pour les amorces STEM, héritant de `Primer3Conserved`

**Fonctionnement** : Convertit les paramètres nommés `stem_primer_*` vers les tags Primer3 internes (`PRIMER_INTERNAL_*`), permettant une interface utilisateur claire.

---

## 4. 🟡 FICHIERS MODIFIÉS — Modules Perl Existants

### 4.1 `lib/Bio/Tools/Run/Primer3.pm` — 2 modifications

| Modification | Détail |
|---|---|
| Ajout de paramètres Primer3 | `PRIMER_TM_FORMULA`, `PRIMER_SALT_CORRECTIONS`, `PRIMER_THERMODYNAMIC_ALIGNMENT`, `SEQUENCE_EXCLUDED_REGION` |
| Fix parsing des résultats | `split '='` → `split('=', $_, 2)` pour gérer les valeurs contenant `=` |

**Justification** : Les nouveaux paramètres permettent d'utiliser la formule de Tm SantaLucia (plus précise pour LAMP à 65°C) et la correction de sel de SantaLucia, essentiels pour un design thermodynamique fiable.

### 4.2 `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm` — Modifications majeures

| Modification | Lignes | Description |
|---|---|---|
| **Fenêtre glissante d'entropie** | +50 lignes | Lissage par fenêtre de la taille du primer, marquage des régions à haute entropie comme `SEQUENCE_EXCLUDED_REGION` pour Primer3 |
| **Diagnostic zéro résultats** | +27 lignes | Affichage des régions exclues et du `PRIMER_*_EXPLAIN` quand aucun oligo n'est trouvé |
| **Suppression du filtre de conservation stricte** | -25 lignes | L'ancien filtre exigeant 100% d'identité sur toutes les séquences du MSA est supprimé |
| **Filtre homopolymères post-Primer3** | +20 lignes | Filtrage manuel des runs de bases identiques (AAAA, CCCC...) car Primer3 ne les filtre pas toujours |
| **Helper IUPAC** | +17 lignes | Fonction `_getIUPAC()` pour la conversion bases → code IUPAC |
| **Code IUPAC désactivé (commenté)** | +23 lignes (commentées) | Génération de séquences dégénérées directement dans l'énumérateur — désactivée au profit du Validator |

**Justification biologique** : La suppression du filtre de conservation 100% est le changement le plus important. L'original exigeait que chaque amorce soit parfaitement identique dans TOUTES les séquences du MSA — rendant impossible le design sur des virus RNA hautement variables. Le nouveau système délègue la gestion de la diversité au `Validator.pm` via les codes IUPAC.

### 4.3 `lib/LLNL/LAVA/PrimerSet/LAMP.pm` — 1 ajout

| Modification | Description |
|---|---|
| **`getStemLocationSummary()`** | +68 lignes — Nouvelle méthode retournant les positions génomiques des primers STEM (FSTEM/BSTEM) avec gestion strand sense/antisense |

### 4.4 `lib/LLNL/LAVA/PrimerSetAnalyzer/PCRPair.pm` — 1 modification

| Modification | Description |
|---|---|
| **Import explicite** | +1 ligne — `use LLNL::LAVA::PrimerSetInfo::PCRPair;` ajouté pour garantir le chargement du module |

---

## 5. 🟢 FICHIERS AJOUTÉS — Interface Web

### 5.1 `lava_flask_app.py` [NOUVEAU] — 1433 lignes
**Rôle** : Application web Flask complète pour le design LAMP

**Fonctionnalités** :
- Upload de fichiers FASTA (jusqu'à 1 Go)
- Configuration interactive de ~40 paramètres LAVA (outer, middle, inner, stem/loop primers)
- Lancement asynchrone des scripts Perl via `subprocess`
- Monitoring en temps réel des exécutions
- Téléchargement des résultats (signatures, FASTA amplifiés/exclus)
- Interface bilingue (français/anglais) via dictionnaire `TRANSLATIONS`
- Authentification par session
- Historique des exécutions passées

### 5.2 `launch_lava_smart_kill.py` [NOUVEAU] — 308 lignes
**Rôle** : Lanceur intelligent avec arrêt automatique à la fermeture du navigateur (environnement WSL)

### 5.3 Templates HTML (5 fichiers, 1542 lignes)
| Fichier | Lignes | Rôle |
|---|---|---|
| `templates/base.html` | 128 | Template de base Jinja2 |
| `templates/index.html` | 592 | Page principale — formulaire de configuration |
| `templates/monitor.html` | 510 | Monitoring temps réel des exécutions |
| `templates/executions.html` | 199 | Historique des exécutions |
| `templates/login.html` | 113 | Page d'authentification |

---

## 6. 🟢 FICHIERS AJOUTÉS — Déploiement Production

### 6.1 `deployment/` [NOUVEAU] — 4 fichiers (394 lignes)
| Fichier | Rôle |
|---|---|
| `deploy.sh` | Script de déploiement automatisé (193 lignes) |
| `gunicorn_config.py` | Configuration Gunicorn (workers, timeouts, bind) |
| `nginx_lava.conf` | Configuration Nginx (reverse proxy, SSL, limites upload) |
| `lava-dna.service` | Unité systemd pour service Linux |

---

## 7. 🟢 FICHIERS AJOUTÉS — Documentation

| Fichier | Taille | Rôle |
|---|---|---|
| `README.md` | 10 Ko | Documentation principale du fork (Markdown) |
| `README_Interface.md` | 3.4 Ko | Guide de l'interface web |
| `DOCUMENTATION_LAVA.txt` | 13 Ko | Documentation technique détaillée |
| `LAVA_PARAMETERS_REFERENCE.txt` | 9 Ko | Référence complète des ~40 paramètres |
| `LAVA_EVOLUTION_JOURNAL.md` | 50 Ko | Journal d'évolution du projet |
| `Makefile` | 30 Ko | Makefile de build complet (vs. `Makefile.PL` seul dans l'original) |

---

## 8. 🟡 FICHIERS MODIFIÉS — Configuration

### 8.1 `.gitignore` — Réécriture complète
- **Original** : Liste de fichiers de build Perl (blib/, Makefile, etc.)
- **Fork** : Gitignore complet pour projet Python+Perl (venv, __pycache__, .DS_Store, logs, IDE, résultats)

---

## 9. 📋 Synthèse des Innovations Biologiques et Algorithmiques

### 9.1 Gestion de la Diversité Génomique
| Aspect | Original | Fork |
|---|---|---|
| Stratégie | Conservation stricte 100% | Tolérance IUPAC avec seuils configurables |
| Virus cibles | Conservés uniquement | Hautement variables (Dengue, etc.) |
| Couverture génomique | Tout ou rien | Pourcentage configurable (défaut 70%) |
| Protection 3' | Non | Zone 3' protégée contre les mismatches |
| Filtrage bruit | Non | `min_base_frequency` pour ignorer les variants rares |

### 9.2 Modèle Thermodynamique
| Aspect | Original | Fork |
|---|---|---|
| Pénalité d'espacement | Parabole | Sigmoïde généralisée (plateau + montée logistique) |
| Formule Tm | Basique | SantaLucia (via `PRIMER_TM_FORMULA=1`) |
| Correction de sel | Basique | SantaLucia (via `PRIMER_SALT_CORRECTIONS=1`) |
| Géométrie | Fixe | Proportionnelle à la longueur totale de signature |

### 9.3 Architecture Logicielle
| Aspect | Original | Fork |
|---|---|---|
| Script principal | 1 monolithique | 2 spécialisés (STEM + LOOP) |
| Modules utilitaires | 0 | 4 nouveaux (`Core.pm`, `Validator.pm`, `Lava::Core`, `StemConserved.pm`) |
| Interface utilisateur | CLI uniquement | CLI + Interface Web Flask |
| Déploiement | Manuel | Automatisé (Nginx + Gunicorn + systemd) |
| Export des résultats | Fichier texte brut | Rapports par signature + FASTA amplifiés/exclus |
| Validation post-design | Aucune | Intersection de couverture + combinatoire |

---

## 10. 📁 Arborescence Comparative

```
pseudogene/lava-dna (ORIGINAL)          Fork LAVA (MODIFIÉ)
├── .gitignore                          ├── .gitignore [MODIFIÉ]
├── .travis.yml                         ├── .travis.yml
├── MANIFEST                            ├── MANIFEST
├── Makefile.PL                         ├── Makefile.PL
├── README                              ├── README
├── environment.yml                     ├── environment.yml
├── lava.pl [SUPPRIMÉ]                  ├── lava_stem_primer.pl [NOUVEAU]
├── slava.pl [SUPPRIMÉ]                 ├── lava_loop_primer.pl [NOUVEAU]
│                                       ├── lava_flask_app.py [NOUVEAU]
│                                       ├── launch_lava_smart_kill.py [NOUVEAU]
│                                       ├── Makefile [NOUVEAU]
│                                       ├── README.md [NOUVEAU]
│                                       ├── README_Interface.md [NOUVEAU]
│                                       ├── DOCUMENTATION_LAVA.txt [NOUVEAU]
│                                       ├── LAVA_PARAMETERS_REFERENCE.txt [NOUVEAU]
│                                       ├── LAVA_EVOLUTION_JOURNAL.md [NOUVEAU]
│                                       ├── .streamlit/ [NOUVEAU]
│                                       ├── deployment/ [NOUVEAU]
│                                       │   ├── deploy.sh
│                                       │   ├── gunicorn_config.py
│                                       │   ├── nginx_lava.conf
│                                       │   └── lava-dna.service
│                                       ├── templates/ [NOUVEAU]
│                                       │   ├── base.html
│                                       │   ├── index.html
│                                       │   ├── monitor.html
│                                       │   ├── executions.html
│                                       │   └── login.html
│                                       ├── static/ [NOUVEAU]
├── lib/                                ├── lib/
│   ├── Bio/Tools/Run/Primer3.pm        │   ├── Bio/Tools/Run/Primer3.pm [MODIFIÉ]
│   └── LLNL/LAVA/                      │   ├── LLNL/LAVA/
│       ├── Constants.pm                │   │   ├── Constants.pm
│       ├── Oligo.pm                    │   │   ├── Oligo.pm
│       ├── OligoEnumerator.pm          │   │   ├── OligoEnumerator.pm
│       ├── OligoEnumerator/            │   │   ├── OligoEnumerator/
│       │   └── Primer3Conserved.pm     │   │   │   └── Primer3Conserved.pm [MODIFIÉ]
│       ├── Options.pm                  │   │   ├── Options.pm
│       ├── PrimerAnalyzer.pm           │   │   ├── PrimerAnalyzer.pm
│       ├── PrimerAnalyzer/PCRPrimer.pm │   │   ├── PrimerAnalyzer/PCRPrimer.pm
│       ├── PrimerInfo.pm               │   │   ├── PrimerInfo.pm
│       ├── PrimerSet.pm                │   │   ├── PrimerSet.pm
│       ├── PrimerSet/                  │   │   ├── PrimerSet/
│       │   ├── LAMP.pm                 │   │   │   ├── LAMP.pm [MODIFIÉ]
│       │   └── PCRPair.pm              │   │   │   └── PCRPair.pm
│       ├── PrimerSetAnalyzer/          │   │   ├── PrimerSetAnalyzer/
│       │   └── PCRPair.pm              │   │   │   └── PCRPair.pm [MODIFIÉ]
│       ├── PrimerSetInfo.pm            │   │   ├── PrimerSetInfo.pm
│       ├── PrimerSetInfo/PCRPair.pm    │   │   ├── PrimerSetInfo/PCRPair.pm
│       └── TagHolder.pm               │   │   ├── TagHolder.pm
│                                       │   │   ├── Core.pm [NOUVEAU]
│                                       │   │   └── Validator.pm [NOUVEAU]
│                                       │   └── Lava/ [NOUVEAU]
│                                       │       ├── Core.pm
│                                       │       └── Enumerator/
│                                       │           └── StemConserved.pm
├── t/                                  ├── t/
└── t_data/                             └── t_data/
```

---

*Document généré le 24/04/2026 — Analyse bioinformatique comparative LAVA-DNA*
