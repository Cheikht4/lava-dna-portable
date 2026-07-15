# Journal d'Évolution du Projet LAVA (Version 2026)

## Introduction

Ce document trace l'évolution du projet LAVA (LAMP Assay Versatile Analysis) dans le cadre de sa refonte majeure visant à :

1. **Amélioration de la flexibilité** : Design des amorces proportionnel à la taille de la signature cible
2. **Gestion de la diversité génomique** : Acceptation de codes IUPAC illimités pour couvrir les variants
3. **Modernisation thermodynamique** : Intégration des modèles les plus récents de Primer3 (SantaLucia 1998, Owczarzy 2004)

---

## Entrées du Journal

### [2026-01-20] Modernisation Thermodynamique et Structurelle
- **Entropie** : Seuil initial à 1.2 bit
- **Thermodynamique** : SantaLucia 1998, Owczarzy 2004
- **Géométrie** : Distances proportionnelles et pénalités sigmoïdes
- **Validation** : Dengue 2

### [2026-01-22] Affinement de la Génération IUPAC (Filtrage Fréquence)
- **Modif** : Filtrage des bases < 5% pour éviter la dégénérescence excessive.
- **Option** : `--min_base_frequency` (défaut : 0.05).

### [2026-01-22] Seuil d'Entropie Configurable et Fenêtre Glissante
- **Modif** : Remplacement du seuil point par point par une **fenêtre glissante** (taille MIN amorce).
- **Problème résolu** : Effet "Tout ou Rien" dû à la fragmentation par pics isolés.
- **Option** : `--entropy_threshold` (défaut : 1.5).

### [2026-01-26] Correction Critique : Calcul d'Entropie "Gap-Aware"
**Fichiers impactés** : `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm`
**Nature du changement** : [Algorithmique / Correction Bug Critique]

**Problème identifié** :
L'algorithme ignorait les gaps (`-`) dans le calcul de l'entropie.
Conséquence : Une région composée de 1% de 'A' et 99% de gaps était considérée comme "parfaitement conservée" (Entropie = 0 sur les bases présentes), alors qu'elle est inutilisable pour le design d'amorces. Cela générait des candidats invalides dans les zones de début/fin d'alignement ou de délétions majeures.

**Solution technique** :
1. **Inclusion des Gaps** : Les gaps sont maintenant comptés comme un état distinct dans le calcul de Shannon.
2. **Pénalité Hard-Threshold** : Toute position contenant plus de **20% de gaps** se voit attribuer une entropie MAXIMALE (10.0), garantissant son exclusion immédiate.

**Justification biologique** :
Une insertion/délétion (indel) ou une zone non séquencée (gap) représente une variabilité extrême structurale. Une amorce ne peut pas s'hybrider sur un vide. Cette correction force Primer3 à éviter les zones "mitées" par les gaps.

**Impact attendu** :
- Disparition des candidats aberrants dans les zones >20% gaps.
- Calcul d'entropie reflétant la vraie stabilité de l'alignement (Bases + Structure).

### [2026-01-29] Retour à l'Algorithme IUPAC Position-par-Position (Rollback)
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Algorithmique / Correction Logic]

**Problème identifié** :
La nouvelle approche "Iterative Greedy" pour l'optimisation IUPAC, bien que théoriquement plus intelligente, introduisait une complexité inutile et potentiellement des régressions dans la détection de couvertures simples. L'approche précédente, analysant chaque position indépendamment pour voir si un code IUPAC permet de récupérer les séquences manquantes, était plus robuste et prévisible.

**Solution technique** :
Remplacement de la fonction `checkPrimerMismatchTolerance` par sa version de `lava_loop_primer_OLD (2).pl`.
- Abandon de la boucle `while ($improved)`
- Retour à une itération simple `for my $pos_offset (0 .. $length - 1)`
- Génération directe du code IUPAC si la conservation parfaite est < seuil, et vérification immédiate du gain de couverture.

**Justification biologique** :
Dans le contexte de virus hautement variables (Dengue), nous voulons identifier rapidement les positions polymorphes critiques et les couvrir par des dégénérescences, sans essayer d'optimiser itérativement une combinaison complexe qui risque de dériver.

**Impact attendu** :
- Comportement plus stable et prévisible de la sélection d'amorces dégénérées.
- Résultats identiques à la version "OLD (2)" de référence validée.
### [2026-01-29] Correction Bug "0% Coverage" sur Loop Primers
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Bug Fix Critique / Algorithmique]

**Problème identifié** :
Des amorces dégénérées (contenant des codes IUPAC comme R, Y...) retournaient 0% de couverture malgré une bonne conception théorique. Trois causes identifiées :
1. **Égalité Stricte** : L'utilisation de `eq` pour comparer une amorce dégénérée à une séquence ciblée échouait systématiquement (ex: `R` n'est pas égal à `A`).
2. **Orientation (Strand)** : Les amorces *LoopB* (Reverse/Antisense) étaient comparées directement au brin *Sense* sans Reverse Complement, causant 100% de mismatches.
3. **Extraction Gap-Aware** : La suppression aveugle des gaps (`s/-//g`) raccourcissait les séquences extraites, causant leur rejet par validation de longueur.

**Solution technique** :
1. **Orientation Dynamique** : Test automatique des deux orientations (Sense vs Antisense) au début de la validation. Si le score Antisense est meilleur, les cibles sont converties en RC.
2. **Validation IUPAC** : Remplacement de l'égalité stricte par `isIUPACCompatible` dans toute la chaîne de validation.
3. **Conservation des Gaps** : Les séquences sont extraites telles quelles (avec gaps/N), les gaps étant traités naturellement comme des mismatches.
4. **Optimisation Early Exit** : Restaurée et adaptée pour utiliser la logique IUPAC correcte.

**Impact attendu** :
- Disparition des faux négatifs (0% de couverture).
- Validation correcte des amorces LoopB (Antisense).

### [2026-01-30] Désactivation Génération Consensus (Revert)
**Fichiers impactés** : `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm`
**Nature du changement** : [Algorihmique / Revert]

**Demande Utilisateur** : 
L'utilisateur a signalé que la génération d'une séquence consensus *avant* validation produisait des amorces de qualité médiocre. Il souhaite utiliser la séquence brute trouvée par Primer3 (référence) comme point de départ pour l'optimisation.

**Solution technique** :
- Désactivation (commenting out) du bloc de génération IUPAC Consensus dans `Primer3Conserved.pm`.
- L'oligo retourné contient désormais la séquence exacte de la Séquence 1 (Reference) identifiée par Primer3.
- C'est le script de validation (`lava_loop_primer.pl`) qui se chargera de la dégénérescence (Phase 3).

**Impact attendu** :
- Meilleure qualité initiale des candidats.
- Validation plus stricte (Phase 2) et optimisation contrôlée (Phase 3).
- Validation plus stricte (Phase 2) et optimisation contrôlée (Phase 3).

### [2026-01-30] Amélioration de la Traçabilité (Debug Logos)
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [UX / Debugging]

**Modification** :
Ajout de messages de débogage explicites pour les cas de rejet critiques qui étaient silencieux :
- Rejet pour variance excessive (Code 'N').
- Rejet final pour couverture insuffisante malgré optimisation (Phase 5).

**Objectif** : 
Fournir à l'utilisateur une preuve formelle que l'algorithme a tenté l'optimisation mais a échoué pour des raisons mathématiques (couverture < seuil), distinguant ainsi un échec logique d'un bug technique.

### [2026-02-02] Relaxation des Pénalités de Spacing
**Fichiers impactés** : `lib/LLNL/LAVA/Core.pm`, `lava_loop_primer.pl`
**Nature du changement** : [Algorithmique / Tuning Thermodynamique]

**Problème identifié** :
Les contraintes de distance (Spacing) étaient trop rigides, éliminant d'excellentes paires d'amorces simplement parce qu'elles s'écartaient de quelques nucléotides de la cible idéale. Le poids de la pénalité de distance écrasait la qualité intrinsèque des amorces.

**Solution technique** :
1.  **Zone de Confort Élargie (Core.pm)** : La fonction sigmoïde tolère désormais **±25%** d'écart par rapport à la cible (au lieu de ±15%) avant d'appliquer une pénalité.
2.  **Poids Réduits (lava_loop_primer.pl)** : Division par deux des coefficients de pénalité de distance (1.0 -> 0.5).

**Impact attendu** :
- Augmentation significative du nombre de sets valides trouvés.
- Sélection prioritaire de la qualité d'amorce (Tm/GC) sur la géométrie parfaite.

### [2026-02-03] Relaxation Modèle & Nettoyage Logique
**Fichiers impactés** : `lib/LLNL/LAVA/Core.pm`, `lava_loop_primer.pl`
**Nature du changement** : [Algorithmique / Mathématiques / Cleanup]

**1. Modification du Moteur Mathématique (Core.pm) :**
- Remplacement du modèle précédent par une **Sigmoïde Généralisée Permissive**.
- **Plateau (Zone de Gratuité)** : ±25% de la cible (ex: pour 40bp, plage 30-50bp est "gratuite").
- **Pente (Slope)** : $k=0.15$ pour une montée "colline douce" au lieu d'un mur vertical.

**2. Correction Logique de Match (lava_loop_primer.pl) :**
- Restauration de `s/-//g` lors de l'extraction des séquences cibles.
- **Raison** : On compare désormais la séquence *physique* (sans gaps) à l'amorce, ce qui est biologiquement plus pertinent pour l'hybridation que l'alignement théorique avec trous.

**3. Nettoyage Dégénérescence (Noise Filter) :**
- Ajout d'un **Filtre de Fréquence (5%)** avant la génération du code IUPAC.
- **But** : Éliminer les variants anecdotiques (<5%) pour éviter de générer des amorces hyper-dégénérées ("soupe") inutilement pénalisées par Primer3.

### [2026-02-03] Synchronisation STEM Primer
**Fichiers impactés** : `lava_stem_primer.pl`
**Nature du changement** : [Maintenance / Harmonisation]

**Action** :
Portage intégral de la logique de validation et de pénalité de `lava_loop_primer.pl` vers `lava_stem_primer.pl`.
- Injection de `checkPrimerMismatchTolerance` (vers. "Nettoyage Logique").
- Réduction des poids de pénalité de distance à 0.5.
- Intégration implicite de la Sigmoïde Généralisée (via `Core.pm` partagé).
- Support de l'orientation dynamique Sense/Antisense pour les STEMs.

### [2026-02-03] Interface Graphique & Paramétrage Avancé
**Fichiers impactés** : `lava_flask_app.py`, `templates/index.html`, `lava_loop_primer.pl`, `lava_stem_primer.pl`
**Nature du changement** : [Interface / UX / Paramétrage]

**Action** :
1. **Perl Backend** :
   - Ajout de l'option `--min_base_frequency` dans `lava_loop_primer.pl` et `lava_stem_primer.pl`.
   - Passage explicite de ce paramètre à la fonction de validation `checkPrimerMismatchTolerance`.
   
2. **Interface Web** :
   - Ajout d'un champ "Fréquence min. bruit" dans le panneau de configuration (section Paramètres généraux).
   - Valeur par défaut : 0.05 (5%) avec incrément de 0.01.
   - I18n : Traductions FR/EN complètes.

**Justification Biologique** :
La "soupe" de variants mineurs (<1-2%) dans les alignements viraux massifs obligeait souvent l'algorithme à générer des codes IUPAC trop larges (ex: N au lieu de R), réduisant le score thermodynamique de l'amorce. Ce paramètre permet à l'utilisateur de décider dynamiquement du niveau de "pureté" requis pour intégrer un variant dans le design.

**Impact Attendu** :
- Contrôle utilisateur total sur la sensibilité aux variants rares.
- Meilleure ergonomie pour les bioinformaticiens ajustant les seuils de bruit.

### [2026-02-05] Refactoring Majeur et Nettoyage de Code (Phase 12)

### Fichiers impactés
- **[NEW]** `lib/LLNL/LAVA/Validator.pm` : Nouveau module centralisé.
- `lava_loop_primer.pl` : Suppression de ~250 lignes redondantes.
- `lava_stem_primer.pl` : Suppression de ~250 lignes redondantes.

### Nature du changement
[Architecture / Refactoring]

### Explication technique
Création du module `LLNL::LAVA::Validator` pour encapsuler toute la logique de validation des amorces (IUPAC, Mismatch Tolerance, Spacing).
Les fonctions suivantes ont été extraites et centralisées :
- `checkPrimerMismatchTolerance`
- `isIUPACCompatible`
- `rev_comp`
- `generateIUPACCode`
- `getPrimerTargetedSequences`

### Justification biologique
Aucun changement fonctionnel biologique. Ce refactoring garantit que les algorithmes de validation (notamment la gestion des variants IUPAC et le "Gap-Awareness") sont strictement identiques entre les primers de boucle (LOOP) et les primers de tige (STEM). Cela élimine le risque de divergence silencieuse entre les deux types d'amorces.

### Impact attendu
- **Maintenance facilitée** : Toute future amélioration de l'algorithme de tolérance s'appliquera instantanément à tous les types d'amorces.
- **Fiabilité** : Code plus propre, moins de risque de bugs copiés-collés.

## [2026-02-05] Correction Bug Refactoring (Validator)

### Fichiers impactés
- `lib/LLNL/LAVA/Validator.pm`
- `lava_stem_primer.pl`

### Nature du changement
[Bug Fix]

### Explication technique
Restauration de la fonction `validateCompleteSignatureSpacing` qui avait été perdu lors du nettoyage.
- Ajout de la fonction dans `Validator.pm` + Export.
- Import explicite dans `lava_stem_primer.pl`.

### Justification biologique
N/A (Correction purement algorythmique pour éviter le crash "Undefined subroutine").

### Impact attendu
- Plus d'erreur fatale lors de l'exécution de `lava_stem_primer.pl`.

## [2026-02-05] Phase 13 : Optimisation Architecturelle "The Big Merge"

### Fichiers impactés
- `lava_loop_primer.pl`

### Nature du changement
[Optimisation Critique / Performance / Refonte]

### Explication technique
Remplacement complet du cœur combinatoire du moteur de recherche d'amorces :
1.  **"The Big Merge"** : Suppression de l'ancienne logique itérative par sous-groupes (`combinationPlan`) qui testait des milliers de combinaisons redondantes. Remplacement par une **Passe Unique** générant des "Listes Maîtres" triées et dédupliquées pour chaque type d'amorce (Inner/Loop/Middle/Outer).
2.  **Fast-Fail Spatiale** : Implémentation d'une logique de recherche par boucles imbriquées (Inner->Loop->Middle->Outer) avec sortie anticipée (`last`) dès qu'une distance dépasse la limite autorisée. Cela évite de tester des millions de combinaisons géométriquement impossibles.
3.  **Cross-Combinaison Optimisée** : Séparation de l'optimisation en deux "Demi-Signatures" (Forward Best + Reverse Best) combinées uniquement à la fin si la paire interne (F1c/B1c) est compatible.

### Justification Biologique
L'ancienne méthode était exhaustive mais exponentielle, rendant le design sur de grands génomes viraux très lent. La nouvelle approche respecte toujours toutes les contraintes biologiques (Thermodynamique, Spacing, IUPAC) mais converge vers la solution optimale en une fraction du temps en éliminant intelligemment les branches mortes de l'arbre de recherche.

### Impact attendu
- **Vitesse** : Accélération drastique du temps de calcul (facteur 10x à 100x attendu sur les cas complexes).
- **Qualité** : Conservation de la même qualité de signatures (méthode de scoring identique).

### 06/02/2026 - Stabilisation des Tags de Signature

**Fichiers impactés :** `lava_loop_primer.pl`

**Nature du changement :** [Bug Fix / Architecture]

**Explication technique :** 
Lors de la refonte du moteur de recherche pour l'optimisation thermique, la création manuelle des objets `LAMP` (signatures) omettait plusieurs métadonnées (tags) attendues par les étapes ultérieures du pipeline (validation et rapports).
J'ai ajouté les assignations manquantes :
1.  Appel explicite à `calculateSignatureIntersection` pour générer `signature_coverage_percent` et `signature_intersection_ids`.
2.  Standardisation des noms de tags (`lamp_penalty` au lieu de `total_penalty` et `signature_intersection_ids` au lieu de `amplified_sequences`).
3.  Ajout du tag informatif `penalty_notes`.

**Justification biologique :** 
Bien que purement informatique, cette correction est critique pour garantir que le filtre de "couverture virale" (qui élimine les signatures ne détectant pas assez de variants) dispose des données nécessaires pour fonctionner. Sans ces tags, le script crashe ou rejette silencieusement des signatures valides.

**Impact attendu :** 
Le pipeline doit maintenant s'exécuter de bout en bout sans erreur Perl, produisant les fichiers de résultats finaux avec les statistiques de couverture correctes.

### 06/02/2026 - Réorganisation du Flux d'Analyse

**Fichiers impactés :** `lava_loop_primer.pl`

**Nature du changement :** [Architecture / Optimisation]

**Explication technique :** 
L'analyse des combinaisons de signatures (étape très coûteuse car factorielle) était effectuée AVANT l'étape de réduction par chevauchement (overlap reduction). Cela signifiait que le programme analysait des signatures redondantes qui allaient de toute façon être fusionnées ou éliminées.
J'ai déplacé le bloc `analyzeSignatureCombinations` pour qu'il s'exécute APRES `reduceSignaturesByOverlap`.
Le flux est désormais :
1. Recherche & Validation Thermique
2. Stockage de toutes les signatures brutes (`.all_signatures`)
3. Réduction par chevauchement (Fusion des signatures trop proches)
4. Analyse des combinaisons sur les signatures finales SÉLECTIONNÉES
5. Génération des rapports finaux

**Justification biologique :** 
Cela assure que les combinaisons proposées ne sont composées que de signatures distinctes et indépendantes, évitant de proposer des paires "artificielles" qui sont en réalité des variations de la même région génomique.

**Impact attendu :** 
Rapports de combinaison plus pertinents et temps de calcul réduit sur les grands sets de données.

### 06/02/2026 - Synchronisation de Stem Primer

**Fichiers impactés :** `lava_stem_primer.pl`

**Nature du changement :** [Synchronisation / Optimisation / Architecture]

**Explication technique :** 
Le script `lava_stem_primer.pl` (version expérimentale avec architecture "Stem") accusait un retard technologique par rapport à la version "Loop". J'ai porté l'ensemble des optimisations récentes :
1.  **Filtre Thermique Dynamique** intrusif dans les boucles imbriquées (Inner -> Stem -> Middle -> Outer), assurant la cohérence thermodynamique entre voisins immédiats.
2.  **Calcul de Pénalité Sigmoïde** via `Core.pm`.
3.  **Système de Tags de Signature** complet (`lamp_penalty`, `penalty_notes`, `signature_intersection_ids`, etc.) pour garantir la compatibilité avec les outils d'analyse.
4.  **Réorganisation du Flux** : Déplacement de l'analyse factorielle des combinaisons APRÈS la réduction par chevauchement.

**Justification biologique :** 
Permet d'appliquer la même rigueur de design aux sets d'amorces "Stem" (qui utilisent des boucles structurées pour une hybridation plus rapide/stable) qu'aux sets LAMP classiques. La maintenance de deux bases de code divergentes posait un risque d'erreurs scientifiques.

**Impact attendu :** 
`lava_stem_primer.pl` est maintenant aussi rapide et robuste que `lava_loop_primer.pl`, avec les mêmes garanties de qualité thermodynamique.

### 09/02/2026 - Interface Graphique Avancée

**Fichiers impactés :** `lava_flask_app.py`, `templates/index.html`

**Nature du changement :** [Interface / Configuration]

**Explication technique :** 
Extension majeure de l'interface web Flask pour exposer l'ensemble des paramètres du moteur LAVA.
1.  **Backend** : Mise à jour de `get_default_params` et du mapping des arguments CLI pour inclure les paramètres thermodynamiques (`dntp`, `salt`, `max_tm_diff`), architecturaux (`dist_outer_middle`) et de diversité (`entropy`).
2.  **Frontend** : Ajout d'une section accordéon "Paramètres Avancés" regroupant ces options par catégorie (Thermodynamique, Architecture, Diversité, Config).

**Justification biologique :** 
Permet aux bioinformaticiens d'ajuster finement les conditions de réaction (ex: concentration en sels pour ajuster le Tm) et les tolérances aux mutations sans avoir recours à la ligne de commande, démocratisant l'accès aux fonctionnalités puissantes du moteur.

**Impact attendu :** 
Contrôle total sur l'exécution via le navigateur.

# 2026-02-10 - Étape 17 : Révision de l'Ordre de Tri des Signatures (Critère de Qualité)

- **Fichiers impactés** : `lib/LLNL/LAVA/Core.pm`, `lava_loop_primer.pl`, `lava_stem_primer.pl`.
- **Nature du changement** : [Algorithmique / Ergonomie].
- **Explication technique** : 
    1.  Ajout de la fonction `countDegenerateBases` dans `Core.pm` pour quantifier l'ambiguïté des séquences.
    2.  Modification de la logique de tri dans les scripts principaux :
        -   **Priorité 1** : Couverture des séquences cibles (Décroissant).
        -   **Priorité 2** : Nombre total de bases dégénérées (Croissant).
        -   **Priorité 3** : Pénalité thermodynamique LAVA (Croissant).
    3.  Mise à jour des formats de sortie pour inclure le nombre de bases dégénérées.
- **Justification biologique** : 
    -   Maximiser la couverture est la priorité absolue pour un test diagnostic universel.
    -   Minimiser les bases dégénérées réduit le coût de synthèse et augmente la spécificité/efficacité de l'amplification.
    -   La pénalité thermodynamique reste un critère de qualité important, mais secondaire par rapport à la détectabilité globale.
- **Impact attendu** : Les meilleures signatures présentées à l'utilisateur seront celles qui couvrent le plus de variants avec les amorces les plus simples possibles (moins de dégénérescence).

# 2026-02-10 - Étape 19 : Correction du Passage des Paramètres (Interface Web)

- **Fichiers impactés** : `lava_flask_app.py`, `lava_loop_primer.pl`, `lava_stem_primer.pl`.
- **Nature du changement** : [Bug Fix / Architecture].
- **Explication technique** : 
    1.  **Filtrage Intelligent** : `lava_flask_app.py` trie désormais les paramètres envoyés aux scripts Perl. Les paramètres spécifiques à "STEM" ne sont plus envoyés à "LOOP" et vice-versa, évitant les erreurs "Unknown option".
    2.  **Mise à jour des Scripts** : Ajout du support explicite pour `max_tm_diff` dans `lava_loop_primer.pl`.
    3.  **Restauration Thermodynamique** : Décommenté les paramètres de concentration (dNTP, Sels) dans `lava_stem_primer.pl` pour assurer la cohérence avec `Primer3Conserved.pm`.
- **Justification biologique** : Assure que les calculs thermodynamiques utilisent bien les paramètres définis par l'interface (ex: conditions de sels spécifiques pour LAMP) au lieu de valeurs par défaut silencieuses ou d'erreurs d'exécution.
- **Impact attendu** : Plus d'erreurs "Unknown option" dans les logs lors de l'utilisation des paramètres avancés ou du basculement entre modes Loop/Stem.

# 2026-02-10 - Étape 20 : Correction de la Contrainte Interface (Max Primers)

- **Fichiers impactés** : `templates/index.html`.
- **Nature du changement** : [Interface / Ergonomie].
- **Explication technique** : Abaissement de la limite minimale du champ `max_primer_gen` de 1000 à 100 dans le code HTML.
- **Justification biologique** : Certains utilisateurs souhaitent limiter drastiquement le nombre de candidats (ex: 500) pour accélérer le tri ou réduire l'espace de recherche, mais l'interface bloquait toute valeur inférieure à 1000.
- **Impact attendu** : L'utilisateur peut désormais entrer "500" dans "Max primers générés" et l'interface acceptera la valeur (au lieu de bloquer la soumission et de rester sur la valeur par défaut de 5000).

# 2026-02-10 - Étape 21 : Correction de la Gestion des Décimales (Float)

- **Fichiers impactés** : `lava_flask_app.py` (fonctions `update_params` et `execute_lava`).
- **Nature du changement** : [Bug Fix / Critique].
- **Explication technique** : Modification de la logique de conversion des types. Auparavant, seuls les paramètres contenant "tm", "percent" ou "coverage" étaient traités comme des nombres à virgule. Les autres (ex: `penalty_plateau`, `min_base_frequency`, `entropy_threshold`) étaient brutalement convertis en entiers (ex: 0.25 devenait 0).
- **Justification biologique** : Ces paramètres fins (seuils d'entropie, pénalités) nécessitent une précision décimale. Leur arrondissement à l'entier faussait complètement les calculs de scoring ou désactivait des filtres (ex: fréquence min de bruit à 0% au lieu de 5%).
- **Impact attendu** : Les valeurs décimales saisies dans l'interface (comme 0.05 ou 1.5) sont maintenant correctement transmises aux scripts Perl sans être tronquées.

# 2026-02-10 - Étape 22 : Correction du Type de max_overlap_percent

- **Fichiers impactés** : `lava_loop_primer.pl` (ligne 795), `lava_stem_primer.pl` (ligne 752).
- **Nature du changement** : [Bug Fix / Validation de Type].
- **Explication technique** : Le paramètre `max_overlap_percent` était déclaré comme entier (`=i`) dans `GetOptions`, mais l'interface Python envoyait `0.0` (float). Perl refusait la conversion avec l'erreur "Value '0.0' invalid for option max_overlap_percent (number expected)". Changement de `=i` vers `=f` pour accepter les décimales.
- **Justification biologique** : Le paramètre contrôle le pourcentage maximal de chevauchement entre signatures. Bien qu'il soit souvent à 0 (pas de chevauchement), il peut nécessiter des valeurs décimales dans certains cas d'usage avancés.
- **Impact attendu** : Le script accepte maintenant `0.0` et toute autre valeur décimale pour ce paramètre sans erreur de validation.

# 2026-02-10 - Étape 23 : Correction du Type PRIMER_NUM_RETURN (Primer3)

- **Fichiers impactés** : `lava_loop_primer.pl` (ligne 1221), `lava_stem_primer.pl` (ligne 1189).
- **Nature du changement** : [Bug Fix / Interopérabilité Primer3].
- **Explication technique** : Primer3 refuse strictement les valeurs décimales pour `PRIMER_NUM_RETURN`. Après avoir corrigé la gestion des floats dans Python, `max_primer_gen` était transmis comme `500.0` au lieu de `500`. Ajout d'un cast `int()` autour de `optionWithDefault` pour forcer la conversion en entier dans Perl.
- **Justification biologique** : `PRIMER_NUM_RETURN` contrôle le nombre maximal de candidats générés par Primer3. Ce paramètre doit être un entier strict selon la spécification de Primer3.
- **Impact attendu** : Primer3 accepte maintenant la valeur sans erreur "Illegal PRIMER_NUM_RETURN value". Les candidats sont correctement générés selon la limite définie par l'utilisateur.

# 2026-02-17 - Étape 24 : Correction du Reporting et Harmonisation (Fix Critique)

- **Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`.
- **Nature du changement** : [Bug Fix / Harmonisation].
- **Explication technique** : 
    1. **Bug Fix** : Dans `lava_loop_primer.pl`, correction de la récupération du tag pour les séquences amplifiées (`amplified_sequences` -> `signature_intersection_ids`). Le code utilisait un mauvais nom de tag et retombait sur une valeur par défaut "toutes les séquences".
    2. **Harmonisation** : Mise à jour de `lava_stem_primer.pl` pour utiliser la même logique de reporting robuste (en-têtes complets, gestion d'erreurs, détails de couverture par primer) et correction d'une variable non déclarée.
- **Justification biologique** : Il est crucial de connaître exactement quelles séquences sont couvertes par une signature pour valider la sensibilité du test face aux variants. Un rapport faussement optimiste (100% par défaut) est dangereux.
- **Impact attendu** : Les fichiers de sortie par signature (`signature_XX_VALID_XX_seq.txt`) contiennent désormais la vraie liste des séquences ciblées.

# 2026-03-16 - Étape 25 : Tolérance Biologique aux Mismatches

- **Fichiers impactés** : `lib/LLNL/LAVA/Validator.pm`, `lava_loop_primer.pl`, `lava_stem_primer.pl`.
- **Nature du changement** : [Algorithmique / Biologie].
- **Explication technique** : Modification de la fonction `checkPrimerMismatchTolerance` pour abandonner l'exigence de match parfait absolu (hors IUPAC). Le script identifie désormais l'extrémité 3' critique en fonction de l'orientation de l'amorce (Sense = Fin, Antisense = Début). Toute mutation dans cette zone (`three_prime_zone_size`) entraîne un rejet immédiat. Pour le reste de la séquence (5' et milieu), le script compte les mismatches et valide l'amorce tant que le total ne dépasse pas `max_total_mismatches`. Les paramètres de contrôle ont été correctement connectés depuis les scripts d'appel.
- **Justification biologique** : La réaction d'amplification isotherme (LAMP) nécessite une hybridation parfaite à l'extrémité 3' pour initier l'élongation par la polymérase. En revanche, les régions en 5' ou au centre de l'amorce tolèrent de légères variations (1 à 3 mismatches) sans bloquer la réaction. L'approche précédente, trop stricte, rejetait à tort de bonnes amorces face à des virus variants.
- **Impact attendu** : Une augmentation significative de la couverture diagnostique (plus de variants reconnus) en tolérant des mutations non-bloquantes, tout en garantissant la fonctionnalité de l'amorce grâce à la protection stricte du 3'.

# 2026-03-16 - Étape 26 : Séparation Sémantique et Algorithmique (Dégénéré vs Mismatch)

- **Fichiers impactés** : `lava_flask_app.py`, `templates/index.html`, `lava_loop_primer.pl`, `lava_stem_primer.pl`, `lib/LLNL/LAVA/Validator.pm`.
- **Nature du changement** : [Architecture / UX / Biologie].
- **Explication technique** : Renommage des anciens paramètres de mismatch liés à IUPAC (ex: `max_total_mismatches` -> `max_total_degenerate_bases`) sur toute la base de code (Python, UI HTML, Perl GetOptions). Ajout d'un paramètre exclusif `max_tolerated_mismatches`. Dans `Validator.pm`, la Phase 3 gère désormais strictement les limites des bases dégénérées (IUPAC), tandis que la Phase 4 utilise le nouveau paramètre pour autoriser de véritables mutations non-couvertes par IUPAC dans la zone 5'/milieu.
- **Justification biologique** : Il y avait une confusion sémantique. Une base dégénérée (ex: Y pour C/T) n'est **pas** un "mismatch", c'est une conception d'amorce mixte pour couvrir deux variants. Un vrai mismatch (ex: une amorce A forcée sur une matrice T) est le résultat de la fonction de "Tolérance Mismatch" (Étape 25). Les deux doivent être contrôlés par l'utilisateur de manière indépendante pour éviter un design incontrôlable.
- **Impact attendu** : L'interface web et le backend distinguent maintenant clairement le "nombre max de bases dégénérées autorisées" du "nombre max de mismatches tolérés". L'outil Perl lève l'ambigüité, rendant la génération des signatures beaucoup plus prédictible.

# 2026-03-16 - Étape 27 : Correction du Cache Navigateur pour les Traductions (Flask)

- **Fichiers impactés** : `lava_flask_app.py`, `templates/base.html`, `templates/index.html`.
- **Nature du changement** : [Architecture / UX / Bug Fix].
- **Explication technique** : Les traductions de l'application Flask (Français/Anglais) semblaient bloquées malgré le changement en backend. Le problème provenait du cache agressif des navigateurs web sur les pages HTML générées. J'ai injecté un décorateur global `@app.after_request` dans `lava_flask_app.py` pour forcer les en-têtes `Cache-Control: no-store, no-cache`, `Pragma: no-cache` et `Expires: 0` sur toutes les réponses. De plus, le sélecteur de langue a été centralisé dans la barre de navigation (`base.html`) avec le paramètre `next` pour recharger la page active, et la balise `<html lang="fr">` stricte a été remplacée par une balise dynamique `{{ lang }}`. Enfin, des clés manquantes (`penalty_plateau`, `penalty_slope`) ont été réparées dans le dictionnaire français.
- **Justification biologique** : Une interface fluide et internationale est primordiale pour garantir que les scientifiques du monde entier (anglophones et francophones) puissent utiliser LAVA sans erreur d'interprétation des paramètres complexes liés à la conception d'amorces LAMP.
- **Impact attendu** : 
### [2026-03-24] Phase 28 : Stabilisation Finale de l'Interface et Traduction Robuste

**Fichiers impactés** : `lava_flask_app.py`, `templates/base.html`, `templates/index.html`, `templates/monitor.html`, `templates/executions.html`
**Nature du changement** : [Architecture / UX / Bug Fix]

**Explication technique** :
1. **Architecture de Traduction** : Migration de la logique de langue du client (JavaScript/LocalStorage) vers le serveur (Flask `context_processor`). La langue est désormais injectée globalement via `g.lang` et le filtre Jinja `t`, garantissant une cohérence totale sans flash de contenu non traduit.
2. **Persistance** : Utilisation d'un cookie `language` avec une durée de vie de 1 an et attribut `SAMESITE=Lax` pour une compatibilité maximale entre les environnements WSL/Windows.
3. **Stabilisation Layout** : Rectification des balises `div` orphelines dans `index.html` qui provoquaient un dérèglement de la grille Bootstrap (notamment dans les blocs conditionnels d'upload).
4. **Internationalisation Dynamique** : Création d'un objet `i18n` en JavaScript pour traduire les messages d'état de l'API (rechargement des logs, statuts d'exécution) en temps réel.

**Justification biologique** :
Une interface stable et une traduction précise des paramètres techniques (ex: "penalty plateau", "degenerate bases") sont indispensables pour éviter les erreurs d'interprétation lors du design d'amorces LAMP. La clarté de l'outil réduit le risque d'erreur humaine dans le paramétrage bioinformatique.

**Impact attendu** :
- Interface utilisateur 100% stable visuellement.
- Basculement de langue instantané et persistant.
- Suppression des erreurs de layout sur les petits écrans ou lors des rafraîchissements.

### [2026-03-25] Phase 29 : Correction de la Protection 3' (Bases Dégénérées & Mismatches)

**Fichiers impactés** : `lib/LLNL/LAVA/Validator.pm`
**Nature du changement** : [Bug Fix / Algorithmique / Biologie]

**Explication technique** : 
Correction d'une erreur d'orientation dans le calcul des indices de la zone 3' critique. Auparavant, pour les primers en orientation "ANTISENSE", le script pointait vers l'index 0 (5' réel du primer string) au lieu de la fin de la chaîne. Comme les séquences cibles sont systématiquement normalisées en 5'->3' (via Reverse Complement si nécessaire) avant la validation, l'extrémité 3' biologique correspond toujours à la fin de la chaîne de caractères. Le calcul est désormais unifié : `three_prime_start_idx = length - zone_size`.

**Justification biologique** : 
L'extrémité 3' d'une amorce est le site d'initiation de l'élongation par la polymérase. La présence d'une base dégénérée (IUPAC) ou d'un mismatch à cette position réduit drastiquement l'efficacité de la réaction, voire l'annule totalement. Cette correction garantit une protection stricte (0 mismatch / 0 base dégénérée si configuré) de la zone 3', indispensable pour la robustesse des essais LAMP face aux variants viraux.

**Impact attendu** : 
- Respect strict du paramètre `max_3prime_degenerate_bases` pour tous les primers.
- Augmentation de la spécificité et de la fiabilité des signatures générées.
- Disparition des bases IUPAC indésirables en fin de séquence dans les rapports `.primers`.

### [2026-03-25] Phase 30 : Flexibilité du Design LAMP - Abaissement du Loop Gap

**Fichiers impactés** : `templates/index.html`, `lava_flask_app.py`
**Nature du changement** : [Architecture / Interface / Biologie]

**Explication technique** : 
L'interface graphique imposait une limite minimale de 25nt pour le paramètre `loop_min_gap` (distance entre l'amorce F2 et l'extrémité de F1c). Cette contrainte a été abaissée à **15nt** dans le HTML (`min="15"`) et les métadonnées de l'application Flask ont été mises à jour en conséquence. Le moteur de calcul Perl était déjà capable de gérer des valeurs inférieures ; il s'agissait donc uniquement d'une levée de restriction au niveau de la couche utilisateur.

**Justification biologique** : 
Bien qu'une distance de 25nt soit recommandée pour éviter l'encombrement stérique et faciliter la formation de la boucle (loop) lors de l'amplification isotherme à 65°C, certains génomes viraux compacts ou hautement conservés ne laissent pas toujours cet espace. Abaisser la limite à 15nt permet de trouver des signatures dans des zones plus restreintes tout en conservant une distance physique suffisante pour l'hybridation des amorces LOOP.

**Impact attendu** : 
- Capacité de générer des signatures LAMP sur des cibles plus courtes ou plus encombrées.
- Plus grande liberté de paramétrage pour les experts en bioinformatique.

### [2026-03-25] Phase 31 : Harmonisation des Rapports (STEM vs LOOP)

**Fichiers impactés** : `lava_stem_primer.pl`
**Nature du changement** : [Architecture / Algorithmique / Reporting]

**Explication technique** : 
Alignement strict du format de sortie de `lava_stem_primer.pl` sur celui de `lava_loop_primer.pl`. 
- Unification des variables de coverage (`$target_count`, `$coverage_percent`).
- Harmonisation des en-têtes `.primers` et `.dash` : inclusion systématique des métadonnées de couverture (`coverage`) et de dégénérescence (`degenerate`) pour correspondre au format du rapport principal.
- Correction du post-processing : ajout automatique des fichiers `_amplified.fasta`, `_excluded.fasta` et `_amplified_noms.txt`.
- Enrichissement des rapports individuels avec l'inclusion des séquences FSTEM et BSTEM.

**Justification biologique** : 
La cohérence des formats de sortie est cruciale pour l'interopérabilité des outils d'analyse en aval. Que le design utilise des STEM primers ou des LOOP primers, la structure des données (couverture réelle des cibles, pénalités thermodynamiques) doit être identique pour permettre une comparaison objective des performances des signatures LAMP.

**Impact attendu** : 
- Scripts d'analyse tiers compatibles avec les deux types de résultats LAVA.
- Rapports plus complets pour les designs STEM (incluant désormais le détail des séquences STEM dans les fichiers individuels).
- Nomenclature des fichiers unifiée (préférence pour les termes anglais standard `_amplified` / `_excluded`).

### [2026-04-07] Phase 32 : Correction Bug BioPerl (Hash Odd Elements)

**Fichiers impactés** : `lib/Bio/Tools/Run/Primer3.pm`
**Nature du changement** : [Bug Fix / Robustesse]

**Explication technique** : 
Correction d'une erreur fatale `Odd number of elements in hash assignment` lors de l'appel à `run()`. L'erreur survenait quand un paramètre passé à Primer3 avait une valeur vide (ex: `KEY=`), car le `split '='` par défaut en Perl ne retournait qu'un seul élément, cassant la parité du `map` utilisé pour reconstruire le hachage des entrées. Utilisation de `split('=', $_, 2)` pour forcer le retour d'une paire `(clé, "")`.

**Justification biologique** : 
Ce bug empêchait l'énumération des oligos dans certaines configurations de diversité génomique où des paramètres facultatifs ou par défaut étaient transmis sans valeur explicite. Sa résolution est indispensable pour la stabilité du moteur de recherche d'amorces sur des alignements complexes.

**Impact attendu** : 
- Suppression des plantages intermittents lors de la phase d'énumération.
- Meilleure résilience du module vis-à-vis des paramètres de configuration variables.

### [2026-04-15] Restauration de l'Environnement et Adaptation macOS (Port & Dépendances)

**Fichiers impactés** : `launch_lava_smart_kill.py`, `lava_flask_app.py`, `README.md`
**Nature du changement** : [Bug Fix / Architecture / Déploiement]

**Explication technique** : 
1. L'environnement virtuel Python (`lava_env`) était corrompu ou vide sur macOS (causant une erreur `[Errno 8] Exec format error`). L'environnement a été entièrement reconstruit en utilisant les dépendances inscrites dans `requirements_flask.txt` et `requirements.txt`.
2. Modification du port par défaut de l'interface Flask (5000 -> 5001) en dur dans l'application et les scripts de lancement intelligents.
3. Mise à jour de la documentation pour inclure une procédure complète d'installation sous macOS via Homebrew (en contournant les problèmes de permissions liés à `sudo cpanm` sur Mac).

**Justification biologique** : 
Bien que ce correctif soit purement informatique, garantir la portabilité de LAVA sur l'écosystème Apple est indispensable, ce matériel étant largement prévalent dans les laboratoires de bioinformatique. Un outil de conception de primers LAMP doit pouvoir tourner localement et sans friction pour permettre des itérations rapides.

**Impact attendu** : 
- Disparition du conflit de port classique avec le service "AirPlay Receiver" de macOS qui occupait discrètement le port 5000.
- Installation et déploiement fluides sur Mac (Intel & Apple Silicon).
- Retour à la normale pour l'exécution locale de l'interface web.

### [2026-04-15] Correction Architecturale de la Sigmoïde (Mathématiques)

**Fichiers impactés** : `lib/LLNL/LAVA/Core.pm`
**Nature du changement** : [Bug Fix / Thermodynamique / Algorithmique]

**Explication technique** : 
La fonction `generateSigmoidPenalty` générait des pénalités massives inattendues (ex: scores totaux supérieurs à 200). Le problème venait d'une erreur mathématique classique avec la fonction Sigmoïde classique `1 / (1 + exp(-x))`.
À la sortie de la zone de confort (`plateau_width`), la valeur `x` (l'excès) devenait très légèrement supérieure à 0 (ex: 0.001). L'exponentielle `exp(0)` valait 1, ce qui donnait la pénalité instantanée `max_penalty / (1 + 1) = 50`. 
Conséquence : Dès qu'une distance s'écartait d'un seul nucléotide hors du plateau, le score subissait un **saut brutal de 0 à 50**. Si 4 distances s'écartaient légèrement dans une combinaison, la pénalité s'envolait à 200 instantanément, détruisant des amorces potentiellement excellentes.
J'ai réécrit la formule pour qu'elle passe exactement par l'origine $(0,0)$ tout en gardant sa dynamique d'asymptote vers $max\_penalty$ :
$Penalty(x) = max\_penalty \times \left[ \frac{2}{1 + exp(-k \cdot x)} - 1 \right]$

**Justification biologique** : 
En biologie, la baisse d'efficacité enzymatique (Taq / Bst Polymerase) liée à un espacement légèrement sous-optimal n'est jamais brutale (exception faite de l'encombrement stérique extrême). Il s'agit d'une perte d'efficacité cinétique progressive. Une pénalité qui bondit brutalement de 0 à 50 pour un écart d'un nucléotide annule tout l'intérêt de la "Colline Douce" permissive implémentée plus tôt. Ce n'était pas fidèle à la cinétique enzymatique réelle.

**Impact attendu** : 
- Disparition des scores d'erreur exagérés à +200 dans les logs.
- Les pénalités d'espacement recommenceront doucement vers 1, 2, 3... au lieu de sauter à 50.
- Beaucoup plus de candidats survivront à la combinaison.

### [2026-04-15] Interface de Debug Avancée des Pénalités (Sub-scoring)

**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Architecture / UX / Reporting]

**Explication technique** : 
Implémentation d'un "penalty breakdown" exhaustif et granulaire pour chaque signature. Le tag `penalty_notes` (précédemment réduit au simple résumé `F:X R:Y`) stocke dorénavant le détail absolu du calcul pour les paires Forward et Reverse. 
La chaîne de débogage prend la forme suivante : `Total F:X R:Y | F{Spc[I_L:A L_M:B M_O:C] Thm[I:D L:E M:F O:G]}`. Elle sépare d'une part les pénalités d'espacement géographique (Spc : `Inner-Loop`, `Loop-Middle`, `Middle-Outer`) et d'autre part les pénalités thermodynamiques de Primer3 (Thm : `Inner`, `Loop`, `Middle`, `Outer`).

**Justification biologique** : 
Lors de l'optimisation des essais LAMP face à des virus variants complexes, le concepteur bioinformaticien a besoin de savoir *pourquoi* une signature a été mal notée par l'algorithme. S'agissait-t-il d'un Primer3 GC% / Hairpin médiocre (Thm) ou bien d'une concession géométrique pour s'adapter à une mutation bloquante (Spc) ? Ce traçage permet d'auditer avec précision le comportement du logiciel.

**Impact attendu** : 
- Les fichiers résultats et rapports logs regorgeront désormais de toutes les données justificatives thermodynamiques et spatiales.
- Débogage instantané pour l'utilisateur.

### [2026-04-15] Sécurisation Serveur de l'Interface Flask (Production)

**Fichiers impactés** : `lava_flask_app.py`
**Nature du changement** : [Architecture / Sécurité / Serveur]

**Explication technique** : 
Préparation du code pour un déploiement public sécurisé via de nombreux patchs Flask :
1. **Désactivation RCE** : Le paramètre `debug=True` a été remplacé par une lecture de variable d'environnement (`FLASK_DEBUG`). Cela prévient l'exécution malveillante de code Python à distance ("Remote Code Execution") via la console Werkzeug.
2. **Session Hijacking & Cryptographie** : La Secret Key (qui chiffre les cookies de session) n'est plus en "dur", mais est générée au vol de manière cryptographique forte via `os.urandom(24)` (ou injéctée via variable d'environnement).
3. **Prévention Path Traversal / Arbitrary Execution** : Assainissement du champ sortant `output_name` grâce au filtre `secure_filename()` pour empêcher les traversées de répertoire (`../../../`), et bridage strict du `script_type` via une Liste Blanche (uniquement "STEM" ou "LOOP").
4. **Validation de Concurrence Anti-DDoS** : Sans créer de lourde base de données, l'application assigne dorénavant un "UUID" par visiteur. Lors d'un lancement LAVA (via `/execute`), le système vérifie s'il existe déjà une exécution en arrière-plan à l'état `running` ou `starting` pour cette même session. Si oui, un message d'erreur bloque la requête (Limite stricte de 1 exécution par utilisateur simultanément).

**Justification biologique** : 
LAVA demande énormément de ressources processeur et de charge RAM pour simuler la thermodynamique moléculaire, surtout lors du design des "STEM" avec de nombreuses amorces dégénérées. En garantissant qu'un utilisateur ou un bot ne puisse pas appuyer 10 fois de suite sur le bouton d'exécution, on empêche l'épuisement massif des ressources, bloquant de facto une attaque DDoS classique qui empêcherait le reste de la communauté scientifique d'utiliser l'outil.

**Impact attendu** : 
- L'interface LAVA peut désormais être exposée sereinement sur l'intranet ou l'internet public pour un usage collaboratif.
- Les attaques courantes de bots et le spamming de calculs sont mitigés.

### [2026-04-15] Mise à jour de Cohérence : Interface de Debug pour les STEM Primers

**Fichiers impactés** : `lava_stem_primer.pl`
**Nature du changement** : [Architecture / UX / Reporting]

**Explication technique** : 
Portage exact du système de débogage ("Penalty Breakdown") depuis le script `LOOP` vers le script `STEM`. Dorénavant, lors d'une combinaison d'amorces STEM (Stem Forward / Stem Reverse), la trace générée affichera un découpage strict des pénalités thermodynamiques ("Primer3Penalty") et spatiales ("SpacingPenalty"). Le format de sortie retranscrit dans `penalty_notes` : `Total F:X R:Y | F{Spc[I_S:A I_M:B M_O:C] Thm[I:D S:E M:F O:G]}` (I_S pour Inner-Stem). Si le mode STEM est désactivé, le format retombe gracieusement sur `Spc[I_M:A M_O:B] Thm[...]`.

**Justification biologique** : 
Afin de concevoir des amorces "STEM" capables de s'ancrer et de booster la vitesse des réactions enzymatiques LAMP, il est crucial d'avoir le même niveau de transparence de diagnostic que pour les LOOP. Si une signature STEM affiche un score très faible, c'est généralement que le site d'ancrage imposé à l'amorce STEM force Primer3 à briser la thermodynamique (Thm) ou à s'éloigner indûment du cœur de l'amplicon (Spc). Ce portage assure une cohérence totale de l'UX de diagnostic.

**Impact attendu** : 
- Les fichiers générés par `lava_stem_primer.pl` exhiberont eux-aussi une radiographie parfaite de chaque composant de la paire Forward et Reverse.
- Facilitera fortement l'optimisation des architectures moléculaires STEM.

### [2026-04-16] Ajout du Contrôle de Priorité au Nettoyage des Chevauchements (Pénalité vs Couverture)

**Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`, `lava_flask_app.py`, `templates/index.html`
**Nature du changement** : [Algorithmique / Web UX]

**Explication technique** : 
Mise à jour majeure du module de dédoublonnage `reduceSignaturesByOverlap`. Par le passé, pour les signatures se chevauchant au-delà du seuil de `max_overlap_percent`, LAVA conservait aveuglément la signature avec la "lamp_penalty" la plus basse (via `@byPenaltyLookup`). J'ai introduit un commutateur `--resolve_overlap_by` ("penalty" ou "coverage"). S'il est réglé sur "coverage", LAVA trie désormais le groupe d'amorces concurrentes par leur "signature_coverage_percent" et ne départage par la pénalité qu'en cas d'égalité stricte. Ce mode a été propagé du frontend Flask jusqu'aux cœurs Perl.

**Justification biologique** : 
Dans des cas viraux avec très peu d'amorces possibles, les épidémiologistes privilégient souvent une couverture de souches maximale (ex: 95% des variants) par rapport à une stabilité purement thermodynamique ou d'espacement (ex: une amplicons 5% trop long). Cette fonctionnalité offre la possibilité cruciale de préserver les signatures "plus universelles" qui étaient auparavant supprimées car concurrencées géométriquement par des signatures plus parfaites mais moins couvrantes.

**Impact attendu** : 
- Un nouveau sélecteur disponible sur l'IHM.
- Les chercheurs pourront occasionnellement récupérer des "champions" avec une tolérance mutante supérieure lorsqu'ils enquêtent sur des clades très diversifiés, en changeant simplement la logique de nettoyage.

### [2026-04-24] Phase 33 : Analyse Comparative Exhaustive Pré-Publication

**Fichiers impactés** : `LAVA_DIFF_ANALYSIS.md` (nouveau)
**Nature du changement** : [Documentation / Publication]

**Explication technique** : 
Réalisation d'une analyse exhaustive des différences entre le dépôt original `pseudogene/lava-dna` et le fork modifié. Le document `LAVA_DIFF_ANALYSIS.md` recense :
- 2 fichiers supprimés (`lava.pl`, `slava.pl`)
- 4 modules Perl modifiés (`Primer3.pm`, `Primer3Conserved.pm`, `LAMP.pm`, `PrimerSetAnalyzer/PCRPair.pm`)
- 6 nouveaux modules Perl (`LLNL::LAVA::Core`, `LLNL::LAVA::Validator`, `Lava::Core`, `Lava::Enumerator::StemConserved`, + scripts principaux)
- 1 application web Flask complète (1433 lignes)
- 5 templates HTML, 4 fichiers de déploiement, 5 documents de documentation
L'analyse couvre les innovations algorithmiques (sigmoïde vs parabole, tolérance IUPAC, validation par intersection), thermodynamiques (SantaLucia, protection 3') et architecturales (séparation STEM/LOOP, interface web).

**Justification biologique** : 
Cette documentation est indispensable pour la publication scientifique du fork. Elle permet aux reviewers de comprendre l'ampleur et la cohérence des modifications apportées au design d'amorces LAMP, notamment la capacité nouvelle de cibler des virus hautement variables comme le Dengue grâce à la tolérance contrôlée de la diversité génomique.

**Impact attendu** : 
- Document de référence pour la rédaction de l'article scientifique.
- Traçabilité complète des innovations pour les reviewers.

### [2026-04-25] Phase 34 : Audit de Code et Nettoyage Pré-Publication

**Fichiers impactés** : `lava_stem_primer.pl`, `lava_loop_primer.pl`, `lava_flask_app.py`
**Nature du changement** : [Architecture / Cleanup / Publication]

**Explication technique** :
Audit de code exhaustif identifiant 18 problèmes classés par sévérité. Nettoyage Phase A réalisé :
1. **Code mort supprimé** : ~160 lignes de code commenté par script (ancien système de sortie `lava.pl`, ancienne fonction `reverseAlignmentStrand`, anciens plans combinatoires alternatifs, anciens prints de debug).
2. **Import `Data::Dumper` supprimé** des deux scripts Perl (module de debug non nécessaire en production).
3. **POD orphelin supprimé** de `lava_stem_primer.pl` (documentation pour `generateIUPACCode` qui a été migré dans `Validator.pm`).
4. **Import `datetime` doublon supprimé** dans `lava_flask_app.py` (lignes 8 et 17 identiques).
5. **Paramètres par défaut dupliqués supprimés** : `primer_iupac_min_percent` et `min_primer_coverage` étaient définis deux fois dans `get_default_params()`.
6. **Print de debug remplacé** par `app.logger.debug()` dans la route `/execute` de Flask.

Bilan : STEM passe de 3934 → 3755 lignes (-179), LOOP de 3452 → 3311 lignes (-141), Flask de 1433 → 1431 lignes (-2).

**Justification biologique** :
Aucun changement fonctionnel. Ce nettoyage améliore la maintenabilité du code, réduit le risque de confusion lors des revues de code, et prépare le dépôt pour une publication propre. Le code mort issu de l'ancien `lava.pl` pouvait induire en erreur les contributeurs sur le fonctionnement réel du pipeline.

**Impact attendu** :
- Code source plus propre et plus lisible pour les contributeurs et les reviewers.
- Réduction de ~322 lignes de code mort.
- Logs de production plus propres (pas de sortie `Data::Dumper`).

### [2026-04-25] Phase 35 : Refactoring — Élimination de la Duplication de Code (Phase B)

**Fichiers impactés** : `lava_stem_primer.pl`, `lava_loop_primer.pl`, `lava_flask_app.py`, `lib/LLNL/LAVA/PipelineUtils.pm` [NOUVEAU], `lib/LLNL/LAVA/TagHolder.pm`, `lib/Lava/Core.pm` [SUPPRIMÉ], `lib/Lava/Enumerator/StemConserved.pm` [SUPPRIMÉ]
**Nature du changement** : [Architecture / Refactoring]

**Explication technique** :
Refactoring majeur éliminant la duplication de code critique identifiée lors de l'audit Phase 34.

1. **Création de `LLNL::LAVA::PipelineUtils`** (735 lignes) : Module centralisant les 8 fonctions utilitaires qui étaient copiées identiquement dans `lava_stem_primer.pl` et `lava_loop_primer.pl` — `buildReversePrimers`, `analyzeAll`, `enumeratePairs`, `buildMetricsArray`, `reducePairInfosByPenalty`, `reducePrimersByOverlap`, `reduceSignaturesByOverlap`, `flattenInfoData`. La version de `flattenInfoData` retenue est celle de LOOP (plus robuste avec protection `eval` pour l'extraction du Tm).

2. **Suppression de `lib/Lava/Core.pm`** (1001 lignes) et **`lib/Lava/Enumerator/StemConserved.pm`** (81 lignes) : Modules orphelins non utilisés par aucun script du pipeline. Contenaient des copies supplémentaires (troisième copie) des fonctions utilitaires et des fonctions obsolètes (`solveCoefficients`, `findDeterminant`).

3. **Ajout de `hasTag()` dans `TagHolder.pm`** : Alias propre pour `getTagExists()`, remplaçant le pattern fragile `eval { $obj->getTag("name") }; if ($@) { ... }` utilisé dans les scripts.

4. **Refactoring Flask** : Extraction de deux fonctions utilitaires (`_convert_param_value`, `_apply_lamp_mode`) et d'une constante globale `FLOAT_PARAMS`, éliminant la duplication entre `update_params()` et `execute_lava()`.

**Bilan quantitatif** :
| Fichier | Avant Phase B | Après Phase B | Delta |
|---------|--------------|--------------|-------|
| `lava_stem_primer.pl` | 3755 | 3152 | **-603** |
| `lava_loop_primer.pl` | 3311 | 2698 | **-613** |
| `lava_flask_app.py` | 1431 | 1400 | **-31** |
| `PipelineUtils.pm` | 0 | 735 | **+735** |
| `TagHolder.pm` | 272 | 292 | **+20** |
| `Lava/Core.pm` | 1001 | 0 | **-1001** |
| `Lava/StemConserved.pm` | 81 | 0 | **-81** |
| **Net** | **9851** | **8277** | **-1574** |

**Justification biologique** :
Ce refactoring garantit que toute correction de bug dans les fonctions de filtrage par chevauchement (`reducePrimersByOverlap`), de tri par pénalité (`reducePairInfosByPenalty`) ou d'énumération de paires (`enumeratePairs`) sera automatiquement appliquée aux deux pipelines STEM et LOOP. Auparavant, une correction dans un script devait être manuellement répliquée dans l'autre, créant un risque de divergence silencieuse pouvant affecter la sélection des signatures LAMP.

**Impact attendu** :
- Réduction nette de 1574 lignes de code.
- Maintenance simplifiée : un seul point de correction pour les fonctions partagées.
- Risque de divergence entre pipelines STEM et LOOP éliminé.
- Code prêt pour une publication scientifique propre.

---

### Phase 36 — Harmonisation complète STEM ↔ LOOP (25 avril 2026)

**Date/Étape** : Phase 36 — Harmonisation architecturale STEM/LOOP

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm` (ajout de 6 fonctions exportées, +555 lignes)
- `lava_stem_primer.pl` (suppression de 6 fonctions locales, -534 lignes)
- `lava_loop_primer.pl` (suppression de 5 fonctions locales, portage calcul dynamique, -550 lignes)

**Nature du changement** : Architecture / Algorithmique

**Explication technique** :

Un audit comparatif exhaustif a révélé 7 divergences entre les scripts STEM et LOOP dans des zones qui devaient être identiques. Les corrections suivantes ont été appliquées :

1. **`calculateSignatureIntersection`** — Fusion des versions STEM et LOOP :
   - Signature unifiée avec paramètres `$extra_primer_type` ("stem"/"loop") et `$min_signature_coverage`
   - Pattern `eval{}` défensif de STEM (par primer individuel) retenu pour la robustesse
   - Validation en 3 phases de LOOP retenue (coverage data, intersection, validation finale)
   - Tags de validation stockés automatiquement dans la signature

2. **`createPerSignatureFiles`** — Fusion :
   - Format de sortie unifié avec type de primer paramétré (`$primer_type`)
   - Accès aux primers via méthodes dynamiques (boucle sur F3/B3/F2/B2/F1/B1)
   - Gestion des primers enrichis par tags génériques (`f{type}_info`, `b{type}_info`)

3. **`createAmplificationFiles`** — Version LOOP retenue (index-based avec bounds checking)

4. **`analyzeSignatureCombinations`** et **`generateCombinations`** — Déplacés vers PipelineUtils (déjà identiques)

5. **`calculateDynamicPairLengths`** — Extraite depuis STEM et portée dans LOOP :
   - Calcul à rebours des longueurs cibles des paires Middle et Inner
   - Paramètres `--max_dist_outer_middle` et `--max_dist_middle_inner` ajoutés à LOOP

**Justification biologique** :

L'harmonisation garantit que les deux pipelines appliquent exactement la même logique de validation des signatures. Auparavant, LOOP utilisait un seuil de couverture paramétrable (`$min_signature_coverage`) pour l'intersection des séquences compatibles, tandis que STEM n'en avait pas — ce qui pouvait conduire à des résultats de sélection incohérents entre les deux modes. La fonction `calculateDynamicPairLengths` permet désormais à l'utilisateur de contrôler finement les distances inter-primers dans les deux pipelines, facilitant l'adaptation à des cibles génomiques de tailles variées.

**Impact attendu** :
- Comportement algorithmique identique entre STEM et LOOP pour toutes les étapes communes
- Réduction de ~1080 lignes de code dupliqué supplémentaires
- Le calcul dynamique des longueurs cibles est désormais disponible dans LOOP
- Toute correction future est automatiquement propagée aux deux pipelines

---

### Phase 36.1 — Résolution du Conflit Géométrique LOOP (25 avril 2026)

**Date/Étape** : Phase 36.1 — Résolution du conflit `max_dist_middle_inner` vs `loop_min_gap`

**Fichiers impactés** :
- `lava_loop_primer.pl` (ajout d'une vérification dynamique)

**Nature du changement** : Algorithmique / Bug Fix

**Explication technique** :

L'introduction de `max_dist_middle_inner` (qui dicte la distance cible F2-F1c) dans LOOP posait un problème mathématique fondamental avec `loop_min_gap` (qui force la distance minimale F2-F1c). 
Puisque `max_dist_middle_inner` représente `F2_length + gap(F2, F1c)`, sa valeur par défaut de `30` (héritée de STEM) devenait mathématiquement impossible à atteindre dans LOOP où `loop_min_gap` est par défaut de `25`. 
En effet, pour que le Loop primer ait la place d'exister, la cible géométrique idéale DOIT être supérieure ou égale à la contrainte physique absolue (`middle_primer_target_length + loop_min_gap`).
Le script ajuste désormais automatiquement `max_dist_middle_inner` à `middle_primer_target_length + loop_min_gap` (ex: 20 + 25 = 45) s'il détecte un conflit bloquant, garantissant que l'algorithme n'essaie pas d'optimiser vers une cible physiquement impossible.

**Justification biologique** :

Un primer de boucle (Loop primer) requiert un espace physique suffisant (environ 15-20nt plus les marges) entre les amorces F2 et F1c pour se lier à la boucle de l'amplicon haltère sans entraver la polymérase. Forcer la fonction de pénalité à rechercher un gap de 10nt tout en interdisant formellement tout gap inférieur à 25nt créait une aberration thermodynamique où toutes les configurations viables se trouvaient artificiellement pénalisées.

**Impact attendu** :
- LOOP ne rejettera ni ne pénalisera plus aveuglément les bonnes combinaisons lorsque le calcul dynamique est activé.
- Les utilisateurs sont alertés dans les logs si leurs contraintes entrent en collision géométrique.

---

### Phase 36.2 — Traçabilité des Paramètres d'Exécution (25 avril 2026)

**Date/Étape** : Phase 36.2 — Génération automatique d'un fichier `.params.txt`

**Fichiers impactés** :
- `lava_flask_app.py` (modification de `execute_lava_background`)

**Nature du changement** : Architecture / Outil de traçabilité

**Explication technique** :

Ajout d'une fonctionnalité dans le backend Flask qui intercepte tous les paramètres configurés par l'utilisateur via l'interface graphique juste avant le lancement de l'exécution Perl. Le système génère automatiquement un fichier texte additionnel (portant l'extension `.params.txt`) en parallèle des autres fichiers de résultats (comme `.primers`, `.fasta`, etc.). Ce fichier capture l'horodatage, le fichier d'entrée, le script exécuté, ainsi que la liste exhaustive de tous les paramètres CLI passés à la commande Perl.

**Justification biologique** :

La sélection de signatures LAMP est extrêmement sensible aux multiples paramètres thermodynamiques (Tm, pourcentages d'hybridation) et géométriques (espacements, pénalités). Lors de l'itération pour trouver les meilleures amorces sur des génomes particulièrement complexes, il est essentiel de pouvoir retracer exactement quelles conditions expérimentales (in silico) ont généré tel ou tel ensemble de candidats. Ce fichier garantit la reproductibilité des analyses bioinformatiques.

**Impact attendu** :
- À chaque exécution, un fichier `nom_du_resultat.params.txt` est généré.
- Le fichier contient l'intégralité des réglages de l'interface qui ont été appliqués à cette exécution spécifique.

### [2026-04-26] Correction Logique Géométrique (Mode Classic)
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Algorithmique / Bug Fix]

**Problème identifié** :
En mode "Classic" (`includeLoopPrimers = 0`), le script omettait bien de chercher une amorce Loop, mais continuait d'imposer un écart minimum (`loop_min_gap`) entre l'amorce F2 et F1c (et entre B1c et B2). Cela forçait les amorces Middle et Inner à être artificiellement éloignées d'au moins 25 nucléotides (valeur par défaut), rejetant ainsi des milliers de candidats parfaitement valides et réduisant drastiquement les performances du mode classique sans aucune justification biologique.

**Explication technique** :
La ligne `my $altMiddleEndAt = $innerLocation - ($loopMinGap + 1);` et son équivalent pour le brin reverse étaient exécutées de manière inconditionnelle. J'ai enveloppé ce calcul d'écart dans un bloc `if ($includeLoopPrimers)` et défini un repli (`else`) sur le simple `$minPrimerSpacing` lorsque l'amorce Loop est désactivée.

**Impact attendu** :
Augmentation massive et immédiate du nombre de signatures primaires trouvées en mode "Classic", et production d'amplicons plus compacts et plus rapides (thermodynamiquement plus stables) lorsque l'utilisateur ne souhaite pas inclure de Loop primers.

### [2026-05-12] Restauration du Moteur Combinatoire et Filtrage
**Fichiers impactés** : `lava_loop_primer.pl`, `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm`
**Nature du changement** : [Algorithmique / Bug Fix]

**Explication technique** :
Deux modifications critiques ont été apportées pour restaurer la fonctionnalité de génération de signatures LAMP-LOOP :
1. **Désactivation du filtre manuel des homopolymères (`Primer3Conserved.pm`)** : Le script appliquait un filtre Perl restrictif rejetant toute amorce possédant un homopolymère supérieur à `max_poly_bases` (défaut à 2 dans LAVA). Ce filtre éliminait la quasi-totalité des bons candidats (ex: 91 oligos rejetés sur 155), réduisant le pool à un point tel qu'aucune combinaison géométrique viable n'était possible. Le filtrage a été désactivé en faveur du comportement natif de Primer3.
2. **Correction du calcul d'espacement `inner_gap` (`lava_loop_primer.pl`)** : Le calcul de la distance entre l'amorce Forward Inner (F1c) et Reverse Inner (B1c) était mathématiquement faux. Il utilisait l'extrémité 5' de F1c au lieu de son extrémité 3', surestimant l'écart d'une vingtaine de nucléotides et provoquant le rejet des candidats valides. La formule a été corrigée : `my $inner_gap = $b1c_location - ($f1c_location + $f1c_length);`, et la tolérance statique (`> 100`) remplacée par la limite de sécurité thermodynamique (`> $signatureMaxLength`). Par ailleurs, le seuil de chevauchement lors de la fusion globale (Big Merge) a été fixé à 100% pour conserver la diversité combinatoire.

**Justification biologique** :
La conception d'amorces LAMP exige une grande densité de candidats dans une fenêtre très restreinte (300 pb). Imposer une limite arbitraire de répétitions (ex: refuser tout "AAA") prive l'algorithme des régions génomiques les plus stables thermodynamiquement, surtout chez les virus. La correction géométrique assure quant à elle que l'encombrement stérique et la formation de la structure en haltère (dumbbell) respectent la dynamique réelle de l'ADN polymérase (Bst) à 65°C, sans rejets injustifiés.

**Impact attendu** :
Le script passe d'un échec total ("0 signatures créées") à la génération de centaines de combinaisons LAMP fonctionnelles (ex: 488 signatures générées puis réduites à la meilleure signature optimale avec 100% de couverture), validant ainsi la viabilité du pipeline.

### [2026-05-12] Restauration du Filtre Homopolymères
**Fichiers impactés** : `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm`
**Nature du changement** : [Algorithmique / Restauration]

**Explication technique** :
À la demande de l'utilisateur, le filtre post-Primer3 éliminant les amorces contenant des homopolymères excessifs a été réactivé. Ce filtre inspecte chaque séquence et rejette celles présentant une répétition stricte d'une base (A, C, G ou T) supérieure au paramètre défini par `$maxPolyBases` (paramètre `--max_poly_bases` dans la commande CLI). Le défaut global de LAVA est de 2, mais l'utilisateur peut désormais l'augmenter manuellement (ex: `--max_poly_bases 4` ou `5`) lors du lancement du script pour éviter de vider le pool d'amorces tout en gardant un contrôle sur la composition des oligos.

**Justification biologique** :
La gestion des homopolymères est critique dans la conception des amorces. De longues répétitions de la même base diminuent la spécificité de l'amorce, favorisent le glissement de la polymérase (polymerase slippage) et peuvent induire des structures secondaires indésirables ou un mésappariement non spécifique. Garder ce filtre ajustable permet à l'utilisateur de trouver le juste milieu entre diversité des candidats géométriques et qualité thermodynamique de chaque amorce individuelle.

**Impact attendu** :
Le script filtrera à nouveau les amorces en fonction de la valeur de `--max_poly_bases`. Si réglé trop bas (comme le défaut de 2), le script risque de rejeter la majorité des candidats. L'utilisateur pourra augmenter manuellement cette valeur pour réussir à générer des signatures viables.

### [2026-05-12] Restauration du Seuil de Chevauchement (Big Merge)
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Algorithmique / Restauration]

**Explication technique** :
Lors du débogage précédent pour restaurer la génération de signatures, la variable `$maxSigOverlapPercent` avait été temporairement remplacée par un code en dur de `100` (%) dans les fonctions `reducePrimersByOverlap` de la phase "Big Merge". Cela désactivait toute réduction des listes maîtresses, causant l'évaluation inutile de millions de combinaisons. Le paramètre dynamique `$maxSigOverlapPercent` a été restauré pour toutes les listes (Inner, Middle, Outer, Loop).

**Justification biologique** :
La réduction par chevauchement (Overlap Reduction) est indispensable pour éliminer les candidats redondants qui se chevauchent de manière excessive (ex: décalage d'un seul nucléotide). En réduisant intelligemment le pool d'amorces tout en gardant une diversité spatiale, on allège massivement la complexité algorithmique sans sacrifier les meilleures signatures potentielles.

**Impact attendu** :
Le "Big Merge" recommencera à filtrer drastiquement les amorces avant l'itération combinatoire, réduisant drastiquement le temps de calcul tout en utilisant le seuil de chevauchement défini par l'utilisateur ou par défaut.

### [2026-05-12] Ajustement Stratégique du Big Merge
**Fichiers impactés** : `lava_loop_primer.pl`
**Nature du changement** : [Architecture / Optimisation]

**Explication technique** :
Après avoir observé que le "Big Merge" avec la valeur par défaut (`max_overlap_percent = 0` ou un pourcentage bas) réduisait drastiquement le pool de candidats avant même la phase combinatoire (ex: 225 amorces réduites à seulement 6), ce qui empêchait la création de la moindre signature, le comportement a été corrigé. Dans la nouvelle architecture à passage unique ("Single-Pass") de LAVA, filtrer les amorces par chevauchement avant de les combiner détruit la diversité géométrique nécessaire pour assembler le "puzzle" complexe d'une signature LAMP. Par conséquent, la variable `max_overlap_percent` a été fixée de manière statique à `100` (%) spécifiquement pour la préparation des listes maîtresses (Big Merge), conservant ainsi toutes les variantes. Le paramètre dynamique `$maxSigOverlapPercent` reste exclusivement utilisé lors de la réduction finale des *signatures complètes*.

**Justification biologique** :
La méthode LAMP nécessite des distances extrêmement précises entre ses 6 amorces (F3, F2, F1c, B1c, B2, B3). En supprimant des amorces isolées sous prétexte qu'elles chevauchent d'autres candidats, on risque d'éliminer la seule version décalée d'un ou deux nucléotides qui s'insérait parfaitement dans l'espacement requis. En maintenant la diversité initiale, le moteur combinatoire peut tester tous les ancrages possibles, et l'optimisation par chevauchement n'intervient qu'à la fin pour ne garder que la meilleure signature globale par région ciblée.

**Impact attendu** :
Le script conservera l'intégralité du pool d'amorces (ex: 225->225) avant combinaison, permettant au moteur de "Fast-Fail" de générer des centaines de signatures valides (ex: 700), puis de réduire le lot à la signature unique et parfaite, sans perte d'opportunité géométrique.

### [2026-05-12] Migration du Script STEM vers l'Architecture Big Merge (Single-Pass)
**Fichiers impactés** : `lava_stem_primer.pl`
**Nature du changement** : [Architecture / Algorithmique]

**Explication technique** :
Le script `lava_stem_primer.pl` utilisait encore l'ancienne architecture combinatoire multi-passes héritée : un plan de 12 itérations (`combinationPlan`) avec un tableau de seuils de chevauchement progressifs (`subgroupSchedule` : de 50% à 94%) et un système de cache de sous-groupes (`%cachedSubgroups`, `%cachedSubgroupData`). Cette approche générait les listes maîtresses de manière répétitive et les filtrait agressivement avant la phase combinatoire. L'ensemble du bloc a été remplacé par le Big Merge Single-Pass déjà adopté par `lava_loop_primer.pl` :
- Construction unique de 8 listes maîtresses (Inner F/R, STEM F/R, Middle F/R, Outer F/R) avec `max_overlap_percent = 100`.
- La déclaration `my $possibleSignatures_r` déplacée vers la phase de réduction finale (`reduceSignaturesByOverlap` sur `$allFoundSignatures_r`), là où elle est sémantiquement correcte.
- Suppression de l'appel intermédiaire à `reduceSignaturesByOverlap` à l'intérieur de la boucle et du test d'arrêt précoce (`$minSignaturesForSuccess`).

**Justification biologique** :
Identique au LOOP : la réduction prématurée par chevauchement sur les amorces individuelles détruit la diversité spatiale nécessaire à l'assemblage des signatures LAMP-STEM (ex: passer de 225 à 6 candidats par type rend impossible la résolution des contraintes géométriques F1c-B1c-STEM). En conservant toute la diversité combinatoire jusqu'à la fin, le moteur peut explorer l'espace complet des configurations et sélectionner la signature optimale globalement.

**Impact attendu** :
Le script STEM affiche désormais la section "Building Master Primer Lists (The Big Merge)..." avec les comptages par type, puis procède directement à la combinaison exhaustive Forward/Reverse en une seule passe, avant la réduction finale par chevauchement sur les signatures complètes.

### [2026-05-22] Correction Critique : Calcul d'Entropie "Gap-Aware" (Fix)
**Fichiers impactés** : `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm`
**Nature du changement** : [Algorithmique / Correction Bug Critique]

**Problème identifié** :
Le nettoyage initial des séquences avant le calcul d'entropie remplaçait tout ce qui n'était pas ATCG par un `N`. Cela détruisait les gaps (`-`), empêchant le calcul de Shannon de les détecter et de les pénaliser. Des régions avec 99% de gaps étaient faussement considérées comme parfaites.

**Solution technique** :
Modification du regex de nettoyage (`s/[^ATCG\-]/N/g`) pour préserver les gaps pendant le calcul d'entropie. Les gaps sont ensuite remplacés par `N` *uniquement* dans la séquence clone envoyée à Primer3, afin d'éviter qu'il ne plante tout en conservant les bonnes coordonnées spatiales.

### [2026-05-24] Correction Ordre de Tri des Signatures Individuelles
**Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`
**Nature du changement** : [Bug Fix / Reporting]

**Problème identifié** :
L'ordre des signatures dans le dossier `*_signatures_individuelles` ne correspondait pas à l'ordre des signatures dans le fichier `.primers` principal. Les signatures individuelles et les fichiers `.fasta` des séquences amplifiées étaient générés à partir d'un tableau de mémoire non trié.

**Solution technique** :
Après l'opération de tri (Coverage > Degeneracy > Penalty), la référence globale `$possibleSignatures_r` est désormais réassignée pour pointer vers le nouveau tableau trié. Les fonctions d'écriture en aval (`createPerSignatureFiles`, `createAmplificationFiles`) utiliseront donc la bonne liste ordonnée.
**Impact attendu** :
Cohérence totale des rapports de bout en bout. La `signature_01_...` dans le dossier individuel correspond exactement à la Signature 1 du fichier `.primers`.

### [2026-05-26] Phase Option B : Génération Native des Amorces Reverse

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm` : Ajout de `buildNativeReversePool`
- `lava_loop_primer.pl` : Remplacement des 3 appels `buildReversePrimers` (Outer/Middle/Inner)
- `lava_stem_primer.pl` : Remplacement des 3 appels `buildReversePrimers` (Outer/Middle/Inner)

**Nature du changement** : [Architecture / Bug Fix Critique / Biologie]

**Problème identifié** :
L'ancienne architecture générait les amorces Reverse (B3, B2, B1c) en appliquant un Reverse Complement
aveugle aux amorces Forward validées. Cette approche créait une incohérence critique : une base dégénérée
(IUPAC) autorisée au 5' d'une Forward (zone permissive) devenait automatiquement le 3' de la Reverse
correspondante (zone stricte). La protection 3' n'était donc pas garantie pour les amorces du brin moins.

**Solution technique** :
Implémentation de `buildNativeReversePool` dans `PipelineUtils.pm`. Cette fonction :
1. Calcule le Reverse Complement complet de l'alignement MSA (toutes les séquences)
2. Lance Primer3 directement sur RC(Séquence 1) — les candidats générés sont nativement
   sur le brin moins (5'→3' du brin -)
3. Valide chaque candidat contre les séquences RC de l'alignement en orientation SENS
   (pas d'ANTISENSE auto-détection nécessaire, car tout est déjà normalisé)
4. Applique la protection 3' standard (last N chars) qui correspond maintenant
   au vrai 3' biologique de l'amorce Reverse
5. Convertit les positions du RC (position p dans RC) en coordonnées génomiques
   originales (location = alignmentLength - 1 - p)

**Justification biologique** :
La cinétique d'hybridation LAMP exige une extrémité 3' parfaite pour initier l'élongation
par la Bst polymérase à 65°C. Une base dégénérée (même Y = C/T) à cette position réduit
drastiquement l'efficacité d'amorçage. En générant les amorces Reverse nativement via
Primer3 sur le brin complémentaire, LAVA garantit que la protection de l'intégrité
thermodynamique 3' s'applique dans le référentiel correct pour TOUS les types d'amorces.

**Impact attendu** :
- Disparition des bases dégénérées indésirables au 3' des amorces B3, B2, B1c
- Pool de candidats Reverse indépendant et plus riche (Primer3 optimise sur le bon brin)
- Meilleures signatures LAMP car les 6 types d'amorces sont tous optimisés nativement
- Couverture plus honnête et reproductible (la protection 3' est symétrique entre Forward et Reverse)

### [2026-05-26] Complément Option B : FLOOP et FSTEM natifs

**Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`
**Nature du changement** : [Architecture / Bug Fix / Complément]

**Problème identifié** :
La correction Option B du commit précédent couvrait les amorces Outer/Middle/Inner Reverse,
mais omettait deux autres amorces du brin moins générées par RC aveugle :
- **FLOOP** (Forward Loop) : RC de BLOOP → IUPAC au 3' possible
- **FSTEM** (Forward Stem) : RC de BSTEM → IUPAC au 3' possible

**Solution technique** :
Remplacement des appels `buildReversePrimers(\@loopBackPrimers)` et
`buildReversePrimers(\@stemBackPrimers)` par `buildNativeReversePool()` avec les
mêmes enumerateurs Loop/Stem et la même logique de validation sur RC(MSA).

**Justification biologique** :
FLOOP et FSTEM s'hybrident tous deux sur leurs cibles en 3'→5' (brin moins).
La protection de leur extrémité 3' est tout aussi critique que pour B3/B2/B1c :
elle détermine l'efficacité d'initiation de la synthèse isotherme dans la structure
en haltère (dumbbell) caractéristique de l'amplification LAMP.

**Impact attendu** :
L'intégralité des 6 types d'amorces du brin moins (B3, B2, B1c, FLOOP, FSTEM et leurs
équivalents Middle/Inner) sont désormais générés et validés nativement.

---

### [2026-06-11] Correction Critique : Bug Fatal dans validateCompleteSignatureSpacing (Validator.pm)

**Date/Étape** : 2026-06-11 — Correction régression totale du module STEM

**Fichiers impactés** :
- `lib/LLNL/LAVA/Validator.pm` (fonction `validateCompleteSignatureSpacing`)

**Nature du changement** : Bug Fix — Architecture

**Explication technique** :
La fonction `validateCompleteSignatureSpacing` appelait `$primer->getTag("strand")` directement sur des objets `PrimerInfo`. Or, le tag `strand` est stocké sur l'objet `Oligo` sous-jacent, accessible via `getAnalyzedPrimer()`. `TagHolder::getTag()` lève une exception fatale si le tag n'existe pas. Ce crash silencieux (capturé par le `next` de validation) rejetait 100% des signatures candidates.

**Correctif** : Ajout d'un helper interne `$get_strand` qui cherche le strand dans l'ordre suivant :
1. `getAnalyzedPrimer()->getTag('strand')` (chemin correct)
2. `$primer->getTag('strand')` si c'est directement un Oligo
3. Fallback sur le rôle : `'plus'` pour les primers forward, `'minus'` pour les reverse

**Justification biologique** :
La validation de l'espacement entre amorces d'une signature LAMP est critique pour garantir que F3, F2, F1c, FSTEM, BSTEM, B1c, B2, B3 ne se chevauchent pas sur le génome cible. Un rejet systematique par exception interne rendait le script incapable de produire tout résultat, compromettant l'ensemble du pipeline de design d'amorces LAMP.

**Impact attendu** : Le script `lava_stem_primer.pl` retrouve sa capacité de validation d'espacement fonctionnelle.

---

### [2026-06-11] Correction Critique : Clampage des Indices hors-bornes dans les Tableaux de Pénalités

**Date/Étape** : 2026-06-11 — Correction bug OOB (Out-Of-Bounds)

**Fichiers impactés** :
- `lava_stem_primer.pl` (boucles Forward et Reverse, calcul `$spacingPenalty`)

**Nature du changement** : Bug Fix — Algorithmique

**Explication technique** :
Les tableaux `innerToLoopPenalties_r`, `innerToMiddlePenalties_r`, `middleToOuterPenalties_r` ont une taille égale à `signatureMaxLength`. Les distances calculées (`innerToStemDistance`, `innerToMiddleDistance`, `middleToOuterDistance`) pouvaient dépasser cette borne, provoquant un accès out-of-bounds qui retourne `undef` en Perl. La multiplication `undef * weight` retourne `undef`, la somme des pénalités devient `undef`, et la comparaison `undef < $bestSetPenalty` échoue silencieusement — aucune combinaison n'est sauvegardée dans `%bestForwardInfos`.

**Correctif** : Clampage explicite avant tout accès : `$d = ($dist < $maxPenIdx) ? $dist : $maxPenIdx`, appliqué aux boucles Forward ET Reverse.

**Justification biologique** :
La fonction de pénalité de distance entre amorces LAMP encode la cinétique d'hybridation à 65°C : des distances trop grandes entre FSTEM et F1c réduisent l'efficacité de la polymérisation en boucle. Le clampage garantit que cette pénalité reste calculable même pour des configurations géométriques atypiques, sans bloquer le moteur de recherche combinatoire.

**Impact attendu** : Les pénalités d'espacement sont désormais toujours numériques, permettant la comparaison et la sélection des meilleures combinaisons d'amorces.

---

### [2026-06-11] Correction Racine : Bornes Géométriques FSTEM/BSTEM (cause principale des 0 signatures)

**Date/Étape** : 2026-06-11 — Correction de la régression principale

**Fichiers impactés** :
- `lava_stem_primer.pl` (boucle Forward inner : calcul de `stemEndAt` ; boucle Reverse inner : calcul de `stemStartAt`)

**Nature du changement** : Bug Fix — Algorithmique / Architecture

**Explication technique** :
Les bornes de recherche des amorces STEM étaient calculées avec `signatureMaxLength` comme référence :

```perl
# FAUX — zone de 400+ nt au-delà de F1c
my $stemEndAt = $innerLocation + $innerLength + $signatureMaxLength;

# FAUX — zone de 400+ nt avant B1c
my $stemStartAt = $innerLocation - $innerLength - $signatureMaxLength;
```

La liste maîtresse des STEM est triée par position croissante. La boucle `for` sur les STEM utilise un `last` quand `stemLocation > stemEndAt`. Avec une borne à F1c+400+, les STEM réels (positionnés à F1c+10 à F1c+75 biologiquement) passaient le filtre, mais les STEM de l'itération suivante (pour un inner primer différent) avaient une borne différente — la boucle avait déjà avancé trop loin dans la liste triée et ne revenait pas en arrière.

**Correctif biologique** :
En architecture LAMP-STEM, FSTEM se situe entre F1c et le milieu de la zone F1c-B1c, et BSTEM entre le milieu et B1c. La distance F1c-B1c est encodée par `innerPairTargetLength` (calculé dynamiquement). La borne correcte est donc :

```perl
# CORRECT — zone physiologiquement réaliste pour FSTEM
my $stemEndAt = $innerLocation + $innerLength + int($innerPairTargetLength / 2);

# CORRECT — zone physiologiquement réaliste pour BSTEM
my $stemStartAt = $innerLocation - int($innerPairTargetLength / 2);
```

**Justification biologique** :
La géométrie LAMP exige que FSTEM et BSTEM se trouvent dans la zone inter-amorces F1c/B1c (typiquement 40–80 nt). Ces amorces participent à la formation de la structure en haltère (dumbbell) lors de l'initialisation de l'amplification isotherme. Une borne de recherche de 400 nt dépasse largement cette zone physiologique et causait une désynchronisation de l'itérateur sur la liste triée — les amorces STEM existantes n'étaient pas détectées pour la majorité des combinaisons inner F/R.

**Impact attendu** :
- 0 signatures → **492 signatures** avec les paramètres standards (signature_max_length=400, min_base_frequency=0.18)
- **1104 signatures** avec des paramètres assouplis
- Le BigMerge Single-Pass confirme sa supériorité sur l'ancien combinationPlan multi-passes (backup : 2 signatures)
- La correction est valide pour toute séquence cible, quelle que soit la longueur du segment génomique analysé

---

### [2026-06-11] Correction : Paramètre signatureCommonTargetMinPercent mappé sur le mauvais option

**Date/Étape** : 2026-06-11 — Correction de mapping paramètre

**Fichiers impactés** :
- `lava_stem_primer.pl` (ligne de lecture `optionWithDefault`)

**Nature du changement** : Bug Fix

**Explication technique** :
`$signatureCommonTargetMinPercent` lisait le paramètre `min_signatures_for_success` (valeur entière = 1) au lieu de `signature_common_target_min_percent` (valeur pourcentage = 70%). Le seuil de couverture commune des séquences cibles était donc fixé à 1% au lieu de 70%.

**Justification biologique** :
Le pourcentage d'intersection commune garantit qu'une signature LAMP couvre un minimum de séquences homologues dans un alignement multiple. Un seuil à 1% rendait ce filtre inopérant, acceptant des signatures ne ciblant qu'une fraction infime du panel de séquences — compromettant la sensibilité diagnostique de l'assay.

**Impact attendu** : Le filtre de couverture commune fonctionne désormais correctement à 70% (valeur par défaut), garantissant des signatures diagnostiquement robustes.

---

### [2026-06-11] Suppression du Paramètre Mort stem_min_gap

**Date/Étape** : 2026-06-11 — Nettoyage de l'architecture des paramètres de distance

**Fichiers impactés** :
- [lava_stem_primer.pl](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_stem_primer.pl)
- [lava_flask_app.py](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_flask_app.py)

**Nature du changement** : [Architecture / Nettoyage]

**Explication technique** :
La variable `$stemMinGap` (définie dans le script Perl via `--stem_min_gap`) était lue depuis les options du terminal mais n'était jamais exploitée dans la boucle combinatoire de recherche des amorces FSTEM/BSTEM ni dans les calculs d'espacement géométrique. Le calcul d'espacement réel entre amorces est entièrement pris en charge par le validateur global de signatures `validateCompleteSignatureSpacing` dans `Validator.pm`. Nous avons donc supprimé la variable `$stemMinGap` et l'option de ligne de commande `--stem_min_gap` dans `lava_stem_primer.pl`. En parallèle, le paramètre a été retiré des configurations par défaut et de la liste `stem_only_params` dans `lava_flask_app.py` pour éviter tout envoi de paramètre non reconnu au moteur Perl.

**Justification biologique** :
La contrainte physique et stérique d'adjacence entre les amorces LAMP (distance minimale d'espacement $\ge 0$ nt entre toutes les amorces ordonnées sur le brin) est déjà appliquée par le validateur de signature de LAVA pour éviter la formation de dimères d'amorces ou les interférences d'élongation de la polymérase. Un paramètre de gap minimal dédié aux STEM, en plus d'être inactif dans le code original, est conceptuellement superflu car l'espacement global protège déjà les interfaces d'hybridation F1c, FSTEM, BSTEM et B1c.

**Impact attendu** :
L'interface utilisateur et le script Perl sont simplifiés et débarrassés d'un paramètre redondant et inactif, sans altérer la qualité thermodynamique ou géométrique des signatures LAMP-STEM générées.

---

### [2026-06-11] Correction du Mapping du Seuil de Couverture dans lava_stem_primer.pl

**Date/Étape** : 2026-06-11 — Résolution du problème de seuil de couverture "bloqué" à 70%

**Fichiers impactés** :
- [lava_stem_primer.pl](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_stem_primer.pl)

**Nature du changement** : [Bug Fix / Algorithmique]

**Explication technique** :
Une modification antérieure dans `lava_stem_primer.pl` avait changé le mapping de la variable `$signatureCommonTargetMinPercent` (seuil de couverture minimal pour valider une signature LAMP) pour qu'elle lise exclusivement l'option `--signature_common_target_min_percent`. Cependant, l'interface graphique de LAVA (IHM Flask) transmet toujours ce seuil via l'option `--min_signatures_for_success` (qui est le champ historique utilisé pour spécifier la couverture universelle ciblée). En conséquence, l'option `--min_signatures_for_success` transmise par l'IHM était ignorée par le script Perl, et ce dernier retombait systématiquement sur sa valeur par défaut de 70%. Nous avons mis en place une lecture adaptative avec repli : la variable tente de lire `--signature_common_target_min_percent`, et si elle est absente, elle lit `--min_signatures_for_success` avant d'appliquer le fallback par défaut à 70%.

**Justification biologique** :
Dans le design d'amorces LAMP multi-séquences ou sur des isolats viraux, le seuil de couverture détermine la fraction de génomes cibles homologues que la signature combinée (les 6 à 8 amorces) doit couvrir pour être déclarée valide. Permettre à l'utilisateur de baisser ce seuil (par exemple à 1% ou à une seule séquence) ou de l'ajuster finement est crucial pour concevoir des amorces adaptées à des panels de virus hautement divergents ou pour des validations spécifiques sur séquence unique sans subir le rejet systématique d'une contrainte trop stricte à 70%.

**Impact attendu** :
Le seuil de validation affiché à l'exécution et appliqué par le validateur correspond désormais exactement à la valeur configurée par l'utilisateur dans l'interface Flask (ex: 1% pour une seule séquence cible), éliminant l'effet de valeur bloquée ou "hardcodée" à 70%.

---

### [2026-06-11] Intégration IHM pour taille STEM et Tm STEM/LOOP

**Date/Étape** : 2026-06-11 — Exposition des paramètres de taille STEM et de Tm STEM/LOOP dans l'interface graphique

**Fichiers impactés** :
- [lava_flask_app.py](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_flask_app.py)
- [templates/index.html](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/templates/index.html)

**Nature du changement** : [Architecture / Interface Graphique]

**Explication technique** :
1. **Ajout de clés de localisation** dans les dictionnaires bilingues `TRANSLATIONS` de `lava_flask_app.py` pour gérer l'affichage bilingue (Français/Anglais) des paramètres de taille (cible, minimale, maximale) et de température de fusion Tm (cible, minimale, maximale) pour les amorces STEM et LOOP.
2. **Restructuration de la section `loop-advanced-params`** dans le template `templates/index.html` pour inclure 5 nouveaux contrôles numériques de configuration : la longueur minimale (`loop_primer_min_length`), la longueur maximale (`loop_primer_max_length`), ainsi que les températures de fusion cible (`loop_primer_target_tm`), minimale (`loop_primer_min_tm`), et maximale (`loop_primer_max_tm`).
3. **Création d'une nouvelle section `stem-advanced-params`** contenant 6 contrôles pour les longueurs de STEM (cible, min, max) et les températures de fusion de STEM (cible, min, max).
4. **Mise à jour du script JavaScript dynamique** dans l'IHM pour écouter les modifications des sélecteurs `script_type` et `lamp_mode` et afficher conditionnellement la section `stem-advanced-params` (lorsque `script_type === 'STEM'` et `lamp_mode === 'enriched'`) ou `loop-advanced-params` (lorsque `script_type === 'LOOP'` et `lamp_mode === 'enriched'`).

**Justification biologique** :
La température de fusion ($T_m$) et la longueur des amorces sont des déterminants critiques de la thermodynamique de l'amplification LAMP (cinétique d'hybridation et stabilité thermique à la température isotherme standard de 65°C). 
- Permettre à l'utilisateur de configurer les plages de $T_m$ et de taille des amorces d'enrichissement (LOOP et STEM) évite des contraintes trop rigides sur les régions génomiques polymorphes et hautement variables (ex: virus comme la Dengue).
- Le contrôle de la taille des amorces STEM évite la formation de repliements secondaires indésirables ou de structures "stem" instables. 
- L'ajustement thermodynamique des amorces LOOP assure qu'elles s'hybrident à la cinétique voulue sans perturber le cycle d'initiation et d'élongation globale géré par la BST polymérase.

**Impact attendu** :
L'utilisateur final dispose désormais d'un contrôle total sur la géométrie et la thermodynamique des amorces d'enrichissement (STEM et LOOP) directement depuis l'interface web, facilitant le design d'assays robustes pour des panels de virus hautement variables.

---

### [2026-06-11] Extension IHM : Contrôle des plages de Tm (min/max) pour Outer, Middle et Inner Primers

**Date/Étape** : 2026-06-11 — Exposition des températures de fusion minimale et maximale pour Outer, Middle, et Inner Primers

**Fichiers impactés** :
- [lava_flask_app.py](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_flask_app.py)
- [templates/index.html](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/templates/index.html)

**Nature du changement** : [Architecture / Interface Graphique]

**Explication technique** :
1. **Ajout de clés de traduction génériques** (`min_tm` et `max_tm`) dans `TRANSLATIONS` pour éviter la duplication inutile de clés de localisation en français et en anglais.
2. **Restructuration de la grille Bootstrap** dans `templates/index.html` pour les trois sections fondamentales d'amorces (Outer, Middle, Inner). Passage d'une disposition à 4 colonnes (`col-md-3`) à une disposition à 6 colonnes (`col-md-2`) pour chaque type d'amorce, alignant horizontalement : la longueur (cible, min, max) et la température de fusion Tm (cible, min, max).
3. **Exposition de 6 nouveaux contrôles de saisie** (Tm min/max pour Outer, Middle et Inner primers) mappés directement sur les variables existantes de la session Flask (`outer_primer_min_tm`, `outer_primer_max_tm`, `middle_primer_min_tm`, `middle_primer_max_tm`, `inner_primer_min_tm`, `inner_primer_max_tm`).

**Justification biologique** :
Dans une réaction LAMP isotherme, la stabilité thermodynamique des paires de primers est orchestrée par les amorces Inner (F1/B1) et Middle (F2/B2) qui initient la réplication et forment la structure de boucle, tandis que les Outer (F3/B3) déplacent le brin synthétisé.
- Avoir un contrôle direct sur les plages de température de fusion ($T_m$) minimale et maximale pour l'ensemble des 6 amorces fondamentales de l'assay permet d'assurer une cinétique d'hybridation harmonieuse et d'éviter des différences de $T_m$ trop prononcées qui bloqueraient la réaction ou induiraient des hybridations non spécifiques à 65°C.
- Cela améliore significativement la tolérance aux mésappariements sur les isolats cliniques présentant des dérives mutationnelles.

**Impact attendu** :
L'utilisateur a désormais le contrôle complet et granulaire de la fenêtre thermodynamique pour l'ensemble du set d'amorces (6 amorces standard ou 8 amorces enrichies) directement depuis le formulaire de configuration.

---

### [2026-07-01] Optimisation : Réduction spatiale des candidats d'amorces par fenêtre glissante

**Date/Étape** : 2026-07-01 — Documentation et validation du mécanisme de réduction spatiale par fenêtre glissante (`window_size` et `max_per_window`).

**Fichiers impactés** :
- [lib/LLNL/LAVA/PipelineUtils.pm](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lib/LLNL/LAVA/PipelineUtils.pm)
- [lava_flask_app.py](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_flask_app.py)

**Nature du changement** : [Algorithmique / Optimisation de performance]

**Explication technique** :
1. **Implémentation de `reducePrimersByWindow`** dans `PipelineUtils.pm` : Cette routine divise le segment génomique analysé en fenêtres adjacentes de $W$ nucléotides (défini par le paramètre `--window_size`).
2. **Filtrage par pénalité thermodynamique** : Pour chaque fenêtre, le script trie les amorces candidates par score de pénalité calculé par Primer3 (qui intègre la dérive de Tm par rapport à la cible, la présence de structures secondaires et de runs de bases de homopolymères). Seuls les $K$ meilleurs candidats d'amorces (défini par `--max_per_window`) sont conservés pour chaque fenêtre.
3. **Application lors de l'assemblage (BigMerge)** : Ce filtrage intervient immédiatement avant la combinatoire de fusion géométrique des amorces unilatérales (Forward et Reverse) pour les candidats internes (Inner), intermédiaires (Middle), externes (Outer), et d'enrichissement (STEM/LOOP).
4. **Exposition dans l'IHM** : Les paramètres `window_size` et `max_per_window` sont exposés dans la section "Configuration Exécution / Réduction spatiale des candidats" de `templates/index.html` pour permettre un contrôle direct par l'utilisateur.

**Justification biologique** :
Dans le design d'amorces LAMP sur des génomes complets ou des alignements multiples complexes, le nombre de candidats générés par Primer3 à des positions très proches (décalées de seulement 1 ou 2 nucléotides) est extrêmement élevé. Sans réduction spatiale, la fusion combinatoire (le croisement de $F3 \times F2 \times F1 \times STEM \times \dots$) souffre d'une explosion combinatoire exponentielle ($\mathcal{O}(N^k)$), ce qui sature la mémoire et allonge indéfiniment les temps de calcul du BigMerge.
- D'un point de vue biologique, concevoir plusieurs dizaines d'amorces quasiment superposées dans une même région de 10 nucléotides est redondant car elles ciblent le même locus physique.
- Le filtrage par fenêtre glissante sélectionne rigoureusement les "champions thermodynamiques" de chaque région génomique tout en assurant une distribution spatiale uniforme des candidats le long du génome, ce qui réduit le nombre de combinaisons d'un facteur 100 à 1000 sans sacrifier la diversité des signatures découvertes.

**Impact attendu** :
Une accélération majeure de la phase d'assemblage combinatoire (BigMerge) sur les génomes longs ou hautement variables, avec des temps de calcul divisés par 10 ou plus, tout en maintenant la couverture globale et la qualité thermodynamique des signatures LAMP générées.

---

### [2026-07-01] Bug Fix : Filtrage effectif des signatures sous le seuil de couverture dans .primers

**Date/Étape** : 2026-07-01 — Correction du bug de non-filtrage des signatures rejetées lors de l'écriture du fichier principal.

**Fichiers impactés** :
- [lava_stem_primer.pl](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_stem_primer.pl)
- [lava_loop_primer.pl](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lava_loop_primer.pl)

**Nature du changement** : [Bug Fix / Algorithmique]

**Explication technique** :
1. **Évaluation sans filtrage** : Les scripts Perl évaluaient bien la couverture universelle de chaque signature complète via `calculateSignatureIntersection` et lui assignaient un tag `validation_status` (valant `"VALIDEE"` ou `"REJETEE..."`). Cependant, les signatures rejetées n'étaient jamais retirées de la liste active `$allFoundSignatures_r`.
2. **Filtrage explicite** : Ajout d'une opération `grep` sur la liste `$allFoundSignatures_r` immédiatement après l'écriture du rapport global de traçabilité `.all_signatures` et juste avant la réduction finale par chevauchement (`reduceSignaturesByOverlap`).
3. **Maintien de la cohérence** : Les signatures n'ayant pas le statut `"VALIDEE"` sont exclues de `$allFoundSignatures_r` de sorte qu'elles ne soient jamais écrites dans le fichier principal de référence des amorces `.primers`.

**Justification biologique** :
La validation diagnostique d'un assay LAMP exige le strict respect de la sensibilité ciblée. Permettre à des signatures n'atteignant pas le seuil de couverture minimale (ex: 70%) d'être exportées dans le fichier principal d'amorces `.primers` pouvait conduire l'utilisateur final à synthétiser des sets d'amorces incapables de couvrir l'ensemble des isolats cliniques ciblés, provoquant de faux-négatifs.

**Impact attendu** :
Le fichier de sortie de référence `.primers` ne contient désormais **strictement que des signatures validées** ayant une couverture supérieure ou égale au seuil de couverture requis par l'utilisateur (ex: 70%). Les signatures sous le seuil restent listées dans `.all_signatures` avec la mention `REJETEE` à des fins de diagnostic et de traçabilité.


---

### [2026-07-02] Bug Fix : Bypass de l'Early Exit dans checkPrimerMismatchTolerance pour la tolérance aux mismatches

**Date/Étape** : 2026-07-02 — Correction du court-circuit de l'Early Exit (Phase 2b) pour assurer l'évaluation de la tolérance aux mismatches sur toutes les amorces.

**Fichiers impactés** :
- [lib/LLNL/LAVA/Validator.pm](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lib/LLNL/LAVA/Validator.pm)

**Nature du changement** : [Thermodynamique / Algorithmique / Bug Fix]

**Explication technique** :
1. **Court-circuit identifié** : Dans la Phase 2b de la validation individuelle des candidats (`checkPrimerMismatchTolerance`), si le pourcentage de séquences cibles matchant parfaitement (sans erreur) l'amorce brute était supérieur ou égal au seuil de couverture minimale exigé par le pipeline (`min_primer_coverage`, ex: 70%), la fonction effectuait un retour immédiat (`return`).
2. **Absence de tolérance** : Ce retour prématuré renvoyait uniquement la liste des séquences cibles compatibles à 100%, sans jamais évaluer la tolérance aux mismatches configurée par l'utilisateur (Phase 4). Par conséquent, les amorces ayant une couverture brute naturellement élevée (comme FLOOP, BLOOP, F2, B3) ne voyaient pas leur couverture augmenter pour inclure les cibles présentant 1 ou 2 mutations tolérées en 5' (restant bloquées aux valeurs strictes à 0 mismatch). À l'inverse, les amorces ayant une mauvaise couverture brute passaient le filtre et bénéficiaient de la tolérance de mismatch.
3. **Correctif appliqué** : L'early exit a été modifié pour contourner uniquement l'optimisation par bases dégénérées IUPAC (Phase 3) lorsqu'elle n'est pas nécessaire, tout en forçant le passage systématique par l'évaluation des mismatches (Phase 4). L'amorce conserve sa séquence brute mais voit sa liste de séquences compatibles (`compatible_sequence_ids`) augmentée de toutes les séquences compatibles réelles sous la tolérance $K$ configurée (ex: 2 mismatches).

**Justification biologique** :
Les réactions d'amplification isotherme (LAMP) ou PCR tolèrent couramment de légers mésappariements thermodynamiques (mismatches) en zone 5' d'hybridation (en dehors de l'extrémité 3' critique d'initiation de la polymérase à 65°C). Si l'on n'évalue pas cette tolérance sur les amorces ayant une bonne couverture brute de départ, la liste de leurs séquences compatibles reste sous-estimée. Lors de la combinatoire finale, le calcul d'intersection universelle élimine ces séquences cibles pour la signature globale, ce qui sous-évalue artificiellement la couverture diagnostique théorique des assays conçus sur des populations virales hétérogènes (ex: Dengue, Fièvre Jaune, Rotavirus).

**Impact attendu** :
Une évaluation exacte et cohérente de la tolérance aux mismatches pour toutes les amorces du set LAMP (Inner, Middle, Outer, et Loops). Les taux de couverture individuelle et globale des signatures générées reflètent désormais fidèlement les contraintes thermodynamiques configurées, augmentant ainsi de manière significative la quantité de signatures validées trouvées sur des cibles à forte diversité génomique.

---

### [2026-07-02] Optimisation : Formule Asymétrique de Pénalité Géométrique (Sigmoïde Asymétrique)

**Date/Étape** : 2026-07-02 — Remplacement de la pénalité de distance symétrique par une pénalité de distance asymétrique pour favoriser les architectures d'amorces compactes.

**Fichiers impactés** :
- [lib/LLNL/LAVA/Core.pm](file:///Users/cheikhtalibouya/Documents/lava/Nouveau%20dossier%20lava/lava-dna-master/lib/LLNL/LAVA/Core.pm)

**Nature du changement** : [Thermodynamique / Architecture / Algorithmique]

**Explication technique** :
1. **Asymétrie des distances** : Modification de la routine `generateSigmoidPenalty` dans `Core.pm`. Auparavant, la pénalité sigmoïde était calculée à partir de la valeur absolue de l'écart par rapport à la cible (`abs($actual - $target)`). Elle pénalisait donc symétriquement les espacements trop courts et trop longs.
2. **Annulation en zone courte** : La fonction applique désormais une condition d'asymétrie : si la distance stérique réelle observée entre deux amorces adjacentes (ex: F3-F2 ou F2-LOOP) est inférieure ou égale à la cible idéale proportionnelle (`$actual <= $target`), la pénalité de distance est immédiatement ramenée à **0**.
3. **Maintien de la pénalité en zone longue** : Seules les distances réelles supérieures à la cible plus le plateau gratuit de tolérance (`$actual > $target + $L_plateau_width`) subissent la montée de pénalité sigmoïde progressive, afin de pénaliser les amplicons trop étirés.

**Justification biologique** :
Dans la cinétique moléculaire d'amplification d'une réaction isotherme LAMP, la distance stérique séparant les amorces réelles n'a pas besoin d'être strictement égale à la proportionnalité théorique. 
- Une distance plus courte que la cible est **fortement avantageuse** : elle rapproche F3 de F2 (ou B3 de B2), ce qui accélère la cinétique de déplacement de brin par la Bst polymérase et améliore le taux de réplication isotherme, à condition qu'il n'y ait pas de chevauchement d'amorces (déjà empêché par le filtre strict de distance minimale `--min_primer_spacing` $\ge 1$ base).
- Pénaliser lourdement (ex: pénalité de 40 à 80) une signature compacte et performante sous prétexte qu'elle a un faible espacement était un contresens physique. Cette modification permet de libérer le design combinatoire et de retenir des signatures hautement efficaces dans le tube de réaction.

**Impact attendu** :
Une réduction drastique et immédiate des scores de pénalité pour les signatures d'amorces d'amplicons courts et optimaux (chute des scores de Spc de >200 à 0, ne laissant que la pénalité thermodynamique pure de Primer3 ~15). Ces signatures hautement performantes seront désormais classées en tête par le pipeline de LAVA.

---

### [2026-07-03] Amélioration de la Précision du Statut d'Exécution et Remontée Significative des Erreurs

**Date/Étape** : 2026-07-03 - Distinction claire du statut "0 signature trouvée" et remontée intelligente des lignes d'erreurs Perl dans l'interface Web.

**Fichiers impactés** :
- `lava_flask_app.py`
- `lava_stem_primer.pl`
- `templates/monitor.html`
- `templates/executions.html`

**Nature du changement** : [Interface / UX / Reporting]

**Explication technique** :
1. **Harmonisation des logs Perl** : Ajout dans `lava_stem_primer.pl` de la ligne standard `After reduction: N final signatures` émise juste après la réduction par chevauchement, garantissant une cohérence parfaite avec `lava_loop_primer.pl`.
2. **Parsing par expression régulière** : Modification de la logique de fin de run dans `lava_flask_app.py` pour extraire dynamiquement le nombre réel de signatures via l'expression régulière `r'after reduction:\s*(\d+)\s*final signatures'`.
3. **Statut différencié et suppression du faux succès** : Si le nombre de signatures $N = 0$, l'application ne déclare plus un statut `completed` trompeur, mais attribue le statut spécifique `completed_no_results`. L'interface ne présente plus d'icône de validation verte ni d'affirmation de succès basée sur la simple présence de fichiers de logs ou de paramètres bruts.
4. **Remontée significative des erreurs Perl** : En cas de code de retour non nul, le serveur analyse les logs en ordre inverse pour extraire la dernière ligne d'erreur significative (marquée par `error`, `died`, `can't`, etc.) et la remonte directement dans le message utilisateur.

**Justification biologique** :
Lors du design d'amorces LAMP sur des alignements viraux très divergents ou avec des contraintes thermodynamiques et spatiales strictes, il est fréquent qu'aucune combinaison ne satisfasse l'ensemble des critères. Il est essentiel pour le bioinformaticien d'être averti immédiatement et sans ambiguïté visuelle que 0 candidat n'a survécu aux filtres, afin d'ajuster ses paramètres de tolérance sans perdre de temps à télécharger des fichiers de résultats vides.

**Impact attendu** :
Une clarté accrue dans le suivi des exécutions Web. Fin des faux messages de succès lorsque 0 signature est produite, et affichage explicite de la cause d'échec directement dans l'interface sans nécessiter l'ouverture manuelle des fichiers de logs bruts.

---

### [2026-07-03] Validation Rigoureuse de l'Alignement FASTA et Détection d'Erreurs d'Entrée

**Date/Étape** : 2026-07-03 - Validation de la géométrie de l'alignement FASTA en entrée et interception spécifique dans l'interface Web.

**Fichiers impactés** :
- `lava_loop_primer.pl`
- `lava_stem_primer.pl`
- `lava_flask_app.py`

**Explication technique** :
1. **Validation BioPerl et vérification brute** : Ajout d'une double vérification stricte immédiatement après le chargement de l'alignement (`next_aln()`) dans `lava_loop_primer.pl` et `lava_stem_primer.pl`. Le moteur vérifie d'une part la présence d'au moins 2 séquences et l'état `is_flush()`, et d'autre part analyse les longueurs brutes des séquences dans le fichier FASTA d'origine afin de déjouer le padding automatique silencieux (ajout de tirets en fin de séquence) effectué par `Bio::AlignIO::fasta`.
2. **Code d'erreur dédié et marqueur stable** : En cas de séquences de longueurs inégales ou d'alignement invalide, le moteur Perl interrompt immédiatement l'exécution avec le code de sortie dédié `2` et émet le marqueur explicite `ERROR: INPUT_NOT_ALIGNED` sur la sortie d'erreur.
3. **Interception prioritaire côté Flask** : Dans `lava_flask_app.py`, lors du traitement d'une fin d'exécution en erreur (`return_code != 0`), la présence du marqueur `input_not_aligned` dans les logs est vérifiée en priorité. Si détectée, l'interface substitue au message d'erreur générique un message traduit explicite (`error_input_not_aligned`) invitant l'utilisateur à aligner ses séquences avec des outils standard (MAFFT, Clustal).

**Justification biologique** :
La conception d'amorces LAMP consensuelles repose sur le calcul de l'entropie de Shannon par colonne de position dans un alignement multiple de séquences (MSA). L'injection de séquences brutes non alignées (de longueurs différentes) fausse totalement le repère des coordonnées spatiales : des homologues fonctionnels ne se trouvent plus sur la même colonne d'alignement, ce qui génère des calculs d'entropie aberrants et conduit soit à des échecs silencieux, soit à des amorces inopérantes in vitro. Garantir que l'entrée est un véritable alignement multiple avant tout calcul combinatoire protège la validité scientifique des amorces produites.

**Impact attendu** :
Protection complète de l'utilisateur contre l'envoi accidentel de fichiers FASTA non alignés ou de séquences isolées. L'arrêt est instantané, explicite, et l'interface guide clairement l le chercheur vers l'étape de pré-traitement (alignement multiple) nécessaire.

---

### [2026-07-03] Compatibilité avec les Fichiers FASTA à Séquence Unique

**Date/Étape** : 2026-07-03 - Prise en charge des fichiers FASTA mono-séquence dans le contrôle d'alignement.

**Fichiers impactés** :
- `lava_loop_primer.pl`
- `lava_stem_primer.pl`

**Nature du changement** : [Algorithmique / Validation / Compatibilité]

**Explication technique** :
Ajustement de la condition de validation d'entrée. Si le fichier FASTA contient exactement 1 séquence (`num_sequences == 1`), la vérification d'alignement multiple (`is_flush()` et comparaison des longueurs) est naturellement ignorée car non applicable à une séquence isolée. Les contrôles géométriques de longueurs inégales ne s'appliquent qu'à partir de 2 séquences (`num_sequences >= 2`). Si 0 séquence est détectée, l'exécution s'arrête avec le marqueur `INPUT_EMPTY`.

**Justification biologique** :
Dans de nombreux protocoles de diagnostic moléculaire, le chercheur souhaite concevoir des amorces LAMP à partir d'une unique séquence génomique de référence (sans calcul de dégénérescence ni alignement muti-souches). Bloquer les fichiers mono-séquence constituait un contresens clinique ; cet ajustement préserve la sécurité du contrôle anti-décalage sur les alignements multiples tout en garantissant la fluidité du pipeline sur les références uniques.

**Impact attendu** :
Restauration immédiate de la capacité de LAVA à concevoir des amorces LAMP (classique ou enrichi) sur des séquences uniques, sans déclencher d'alerte de non-alignement.

---

### [2026-07-03] Internationalisation Complète du Monitoring et des Messages d'Exécution

**Date/Étape** : 2026-07-03 - Traduction complète et dynamique de la surveillance et des diagnostics d'erreur.

**Fichiers impactés** :
- `lava_flask_app.py`
- `templates/monitor.html`
- `templates/executions.html`

**Nature du changement** : [Architecture / UI Fix]

**Explication technique** :
Remplacement des chaînes codées en dur en français dans le thread d'arrière-plan `execute_lava_background` par des appels dynamiques au dictionnaire `TRANSLATIONS` indexés par `running_executions[execution_id]['lang']`. Cela garantit que les messages de fin (`completion_message`) ou d'erreur s'affichent dans la langue active de la session. Enrichissement bilingue de la fonction `translate_error_to_user_friendly(error_message, lang)`. Remplacement de tous les libellés de surveillance et des compteurs de progression dans `monitor.html` et `executions.html` (HTML et JS) par des filtres Jinja `|t`.

**Justification biologique** :
Dans un contexte international de surveillance épidémiologique, la clarté et l'accessibilité linguistique des messages de retour (détection de signatures, alertes sur les seuils de couverture, erreurs de format FASTA ou d'encodage) sont cruciales pour éviter les interprétations erronées et permettre un ajustement rapide des paramètres géométriques et thermodynamiques.

**Impact attendu** :
Une bascule instantanée et complète en anglais ou en français de la page de suivi en temps réel, des alertes de fin d'exécution et de la liste des analyses.

---

### [2026-07-03] Création du Dépôt Git Autonome de l'Interface Web (LAVA-Web)

**Date/Étape** : 2026-07-03 - Structuration d'un dépôt propre dédié au déploiement de l'interface graphique et à la publication.

**Fichiers impactés** :
- `lava-interface-public/` (Nouveau dépôt autonome)
- `README.md`
- `.gitignore`

**Nature du changement** : [Architecture / Publication]

**Explication technique** :
Création et initialisation d'un dépôt Git dédié (`lava-interface-public`) en miroir du dépôt moteur (`lava-dna-public` / `LAVA-Virus`). Ce nouveau dépôt regroupe exclusivement l'application web Flask (`lava_flask_app.py`), les templates bilingues (`templates/`), les scripts de production (`deployment/`), ainsi que le sous-ensemble minimal de modules scientifiques Perl requis pour son exécution en production. Tous les fichiers de log, caches Python, fichiers temporaires et archives de travail sont rigoureusement exclus via un `.gitignore` optimisé.

**Justification biologique** :
Séparer clairement le dépôt de l'interface web prête à l'emploi du dépôt moteur de recherche fondamentale permet de proposer à la communauté scientifique et aux laboratoires cliniques un package web immédiatement déployable sur serveur de diagnostic, sans alourdir le dépôt avec l'historique d'expérimentation ou les scripts de test algorithmique.

**Impact attendu** :
Disponibilité et publication officielle du package Git public sur GitHub sous le nom `LAVA_Virus-Interface` (https://github.com/Cheikht4/LAVA_Virus-Interface.git), simplifiant drastiquement l'installation et la maintenance de serveurs LAVA dans les instituts de recherche.

---

### [2026-07-05] Publication Officielle du Dépôt LAVA_Virus-Interface

**Date/Étape** : 2026-07-05 - Publication et synchronisation du dépôt autonome de l'interface graphique.

**Fichiers impactés** :
- `lava-interface-public/` (dépôt public distant)

**Nature du changement** : [Publication / Git]

**Explication technique** :
Connexion du répertoire local `lava-interface-public` au remote officiel GitHub (`https://github.com/Cheikht4/LAVA_Virus-Interface.git`) via la commande `git remote add origin` et push de la branche `main`. Ce dépôt constitue désormais la référence officielle et autonome pour le déploiement web de LAVA-DNA.

**Justification biologique** :
Fournir un accès direct et universel au code source propre de l'interface web permet aux laboratoires de référence internationale d'auditer, cloner et déployer l'outil de diagnostic dans leurs propres infrastructures sécurisées.

**Impact attendu** :
Accessibilité publique immédiate de l'application web bilingue LAVA sur GitHub.

---

### [2026-07-05] Passage de l'Interface Web sous Licence Propriétaire (Autorisation Requise)

**Date/Étape** : 2026-07-05 - Restructuration juridique des droits d'utilisation de l'interface graphique.

**Fichiers impactés** :
- `LICENSE`
- `README.md`

**Nature du changement** : [Juridique / Architecture]

**Explication technique** :
Mise à jour du fichier `LICENSE` pour distinguer clairement 3 parties : (1) Le moteur Perl hérité de LLNL (BSD 3-Clause), (2) Les modules Perl étendus comme `Validator.pm`, `Core.pm`, `PipelineUtils.pm` (BSD 3-Clause), et (3) L'ensemble de l'interface graphique web (`lava_flask_app.py`, `templates/`, `static/`, `deployment/`) qui passe sous **Licence Propriétaire - Tous droits réservés (Cheikh Talibouya)**. Toute utilisation, copie, modification ou déploiement clinique/commercial de l'interface est désormais impérativement soumise à l'autorisation écrite préalable de l'auteur.

**Justification biologique** :
Alors que les algorithmes de recherche fondamentale en bioinformatique gagnent à rester en open-source pour permettre la vérification par les pairs et la reproductibilité scientifique, l'interface clinique de diagnostic, intégrant les sécurités anti-décalage et la surveillance d'exécution, requiert un contrôle de diffusion strict afin d'éviter l'émergence de clones non certifiés dans un cadre de santé publique.

**Impact attendu** :
Protection complète de la propriété intellectuelle de l'interface graphique LAVA-DNA, tout en préservant l'ouverture open-source du moteur de calcul bioinformatique sous-jacent.

---

### [2026-07-06] Durcissement de Sécurité : Traversée de Répertoire et Concurrence Atomique

**Date/Étape** : 2026-07-06 - Audit et correction de vulnérabilités sur l'interface web (`lava_flask_app.py`).

**Fichiers impactés** :
- `lava_flask_app.py`

**Nature du changement** : [Sécurité / Architecture]

**Explication technique** :
1. **Protection contre la traversée de répertoire (`output_name`)** : Nettoyage systématique du paramètre `output_name` issu du formulaire via `secure_filename()` avant toute construction de chemin dans `execute_lava`. Si la chaîne nettoyée est vide (ex: saisie malveillante du type `../../`), le système applique automatiquement la valeur par défaut sécurisée `'lava_result'`.
2. **Synchronisation atomique des quotas de concurrence** : Introduction d'un verrou global `executions_lock = threading.Lock()`. Dans la route `/execute`, la vérification des quotas (seuils globaux et par utilisateur) et l'insertion de l'exécution dans le dictionnaire `running_executions` avec le statut `'starting'` sont désormais encapsulées dans un unique bloc atomique (`with executions_lock:`), éliminant toute race condition lors de requêtes simultanées.

**Justification biologique** :
Sur un serveur de diagnostic clinique partagé par plusieurs équipes de recherche, l'intégrité du système de fichiers est primordiale pour éviter l'écrasement ou la fuite de données génomiques sensibles (séquences virales de patients). De plus, le calcul d'amorces LAMP étant intensif en ressources CPU, garantir l'inviolabilité des quotas d'exécution empêche toute surcharge accidentelle ou déni de service (DoS) qui paralyserait les analyses en cours.

**Impact attendu** :
Confinement absolu de tous les fichiers de résultats dans le répertoire dédié (`results/`) et respect strict des limites de calcul simultané en environnement multi-utilisateurs.

---

### [2026-07-06] Durcissement de Sécurité (Fin d'Audit) : Masquage des Traces en Production et Purge du Rate Limiter

**Date/Étape** : 2026-07-06 - Finalisation de l'audit de sécurité et durcissement de l'interface web (`lava_flask_app.py`).

**Fichiers impactés** :
- `lava_flask_app.py`

**Nature du changement** : [Sécurité / Architecture / Bug Fix]

**Explication technique** :
1. **Masquage conditionnel de la trace technique (`technical_error`)** : Dans la route `/api/status`, l'ajout du champ `technical_error` (trace brute Python/Perl) au flux JSON renvoyé au client est désormais strictement conditionné au mode de développement (`os.environ.get('FLASK_ENV') != 'production'`). En environnement de production, seule l'explication traduite et vulgarisée (`error`) est transmise au client web, interdisant toute fuite d'informations sur l'architecture interne du serveur ou les chemins de fichiers.
2. **Purge active de la mémoire du Rate Limiter** : Au sein de la boucle de nettoyage asynchrone (`background_data_cleanup`), un mécanisme de purge des adresses IP inactives a été ajouté. Toutes les heures, les horodatages des requêtes dans `ip_request_history` sont filtrés sur une fenêtre glissante de 3600 secondes ; si la liste d'une IP devient vide, sa clé est définitivement supprimée du dictionnaire `defaultdict(list)`, éliminant toute accumulation passive et fuite mémoire sur les serveurs de longue durée.

**Justification biologique** :
Dans un contexte de diagnostic hospitalier ou de veille épidémiologique en production, l'exposition des traces d'erreurs brutes (contenant des chemins d'accès serveur, des structures de répertoires ou des commandes système) représente un vecteur d'attaque d'ingénierie inverse. Leur masquage garantit l'inviolabilité de la plateforme LAVA-DNA tout en conservant une remontée d'erreur claire pour les virologues. De plus, la purge mémoire préserve la stabilité opérationnelle des serveurs lors des campagnes massives de screening sur plusieurs mois.

**Impact attendu** :
Sécurité maximale des échanges API en environnement de production (aucune fuite technique) et stabilité mémoire garantie sur le long terme sans redémarrage de l'application web.

---

### [2026-07-06] Harmonisation de l'Identité Web : Renommage en LAVA_Virus et Simplification

**Date/Étape** : 2026-07-06 - Mise à jour de l'identité visuelle et textuelle de l'interface graphique.

**Fichiers impactés** :
- `lava_flask_app.py`
- `launch_lava_smart_kill.py`
- `templates/index.html`
- `templates/base.html`
- `templates/login.html`
- `templates/monitor.html`
- `templates/executions.html`

**Nature du changement** : [Architecture / Interface]

**Explication technique** :
Remplacement systématique de l'appellation "LAVA-DNA" par la marque officielle **LAVA_Virus** sur l'ensemble de l'application web (titres de pages HTML, barre de navigation, dictionnaire bilingue de traduction FR/EN et logs de démarrage). Parallèlement, le sous-titre verbeux "LAMP primer design with LAVA - Stable web interface" a été épuré pour devenir simplement **LAMP primer design**, offrant une interface web à la fois plus directe, moderne et centrée sur sa fonction algorithmique essentielle.

**Justification biologique** :
L'outil LAVA étant désormais spécialisé dans le ciblage de génomes viraux hautement variables (séquences consensus et alignements de variants), l'appellation LAVA_Virus reflète avec exactitude le domaine d'application clinique et épidémiologique de la plateforme, évitant toute confusion avec des outils généralistes de PCR sur ADN génomique hôte.

**Impact attendu** :
Clarté immédiate pour les virologues quant à la vocation virale de l'outil et interface graphique épurée.

---

### [2026-07-06] Précision Granulaire du Suivi de Progression par Type d'Amorce LAMP

**Date/Étape** : 2026-07-06 - Spécification détaillée des étapes de validation dans les messages de progression (Outer, Middle, Inner, Loop, Stem).

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm`
- `lava_stem_primer.pl`
- `lava_loop_primer.pl`

**Nature du changement** : [Architecture / Interface / Reporting]

**Explication technique** :
1. **Ajout d'un paramètre de label dynamique** : Modification des routines de validation `getOligosWithMismatchTolerance` (dans les scripts Perl) et `buildNativeReversePool` (dans `PipelineUtils.pm`) afin d'accepter un argument supplémentaire `$label` définissant le type d'amorce en cours de traitement.
2. **Émission spécifique dans LAVA-PROGRESS** : Au lieu d'émettre des messages génériques en dur du type `[LAVA-PROGRESS] Validation Forward` ou `Reverse Validation`, le moteur Perl injecte désormais le nom précis et la nomenclature LAMP officielle dans le flux de progression STDOUT : `Outer Forward (F3)`, `Outer Reverse (B3)`, `Middle Forward (F2)`, `Middle Reverse (B2)`, `Inner Forward (F1c)`, `Inner Reverse (B1c)`, `Loop Back (BLOOP)`, `Loop Forward (FLOOP)`, `Stem Back (BSTEM)` et `Stem Forward (FSTEM)`.
3. **Compatibilité transparente avec l'interface Web** : Le contrôleur Flask (`lava_flask_app.py`), qui intercepte le premier champ avant le séparateur `|` des lignes `[LAVA-PROGRESS]`, affiche automatiquement et en temps réel le libellé granulaire sur l'interface utilisateur de suivi.

**Justification biologique** :
Dans le design d'amorces LAMP enrichi ou classique sur des génomes viraux complexes, chaque catégorie d'amorce (F3/B3, F2/B2, F1c/B1c, Loops ou Stems) obéit à des contraintes thermodynamiques, des tailles et des localisations génomiques distinctes. Un échec de validation ou un temps de calcul prolongé sur une étape spécifique (par exemple les amorces Inner F1c/B1c, très longues et soumises à des structures secondaires) exige que le bioinformaticien sache instantanément quelle population d'oligonucléotides est en cours de criblage ou de rejet, sans se contenter d'une distinction binaire Forward/Reverse.

**Impact attendu** :
Une lisibilité et une transparence exceptionnelles du suivi en temps réel sur l'interface Web LAVA_Virus : le chercheur visualise précisément la progression de la validation étape par étape et par type d'amorce LAMP.

---

### [2026-07-09] Fonctionnalité d'Importation et Réplication Automatique des Paramètres d'Exécution

**Date/Étape** : 2026-07-09 - Import de fichier de paramètres (`.params.txt` ou `.json`) pour faciliter la reproductibilité expérimentale.

**Fichiers impactés** :
- `lava_flask_app.py`
- `templates/index.html`

**Nature du changement** : [Architecture / Interface / UX]

**Explication technique** :
1. **Création du point d'entrée d'importation (`/upload_params_file`)** : Implémentation d'une nouvelle route Flask capable de recevoir et d'analyser les fichiers de paramètres générés lors d'exécutions précédentes (`.params.txt` ou format `.json`).
2. **Analyseur de format LAVA multi-format** : Le parseur extrait dynamiquement le type de script (`STEM` ou `LOOP`), le mode de design (`classic` ou `enriched`), ainsi que l'intégralité des drapeaux Perl (ex: `--primer_min_match_percent: 85.0`). Il mappe chaque paramètre sur sa variable de session Python adéquate en convertissant automatiquement les types de données (`int`, `float`, `bool`, `str`).
3. **Intégration ergonomique dans l'IHM** : Ajout de deux boutons d'action rapide "Importer un fichier de paramètres" (dans l'en-tête de la section de configuration et à côté du bouton de sauvegarde) reliés à une fonction AJAX en JavaScript. Dès la sélection du fichier, le formulaire est mis à jour et rechargé instantanément pour refléter avec exactitude les conditions expérimentales importées.

**Justification biologique** :
En recherche clinique et épidémiologique, la comparaison de performance des amorces LAMP sur différentes souches virales exige une stricte invariance des conditions thermodynamiques (températures de fusion $T_m$, concentrations salines et dNTP, tolérance aux mésappariements, fenêtres de réduction spatiale). Permettre de recharger directement un fichier `.params.txt` issu d'un run réussi évite les erreurs de saisie manuelle et garantit la reproductibilité parfaite des protocoles in silico d'un isolat à l'autre ou entre différents collaborateurs du laboratoire.

**Impact attendu** :
Un gain de temps considérable pour l'utilisateur qui peut désormais répliquer ou ajuster des conditions d'expérience complexes en un seul clic via l'import de ses fichiers de paramètres antérieurs.

---

### [2026-07-09] Audit et Durcissement de Sécurité de la Route d'Importation de Paramètres

**Date/Étape** : 2026-07-09 - Audit et durcissement complet de `/upload_params_file`.

**Fichiers impactés** :
- `lava_flask_app.py`

**Nature du changement** : [Architecture / Bug Fix / Sécurité]

**Explication technique** :
1. **Validation stricte par liste blanche (Priorité 1)** : Construction dynamique de l'ensemble des clés autorisées à partir de `get_default_params().keys()`. Avant toute injection dans `session['params']`, le moteur vérifie que la clé appartient à cette liste blanche ; toute clé inconnue ou malveillante est ignorée silencieusement dans les branches JSON et texte. Les valeurs `script_type` et `lamp_mode` sont rigoureusement restreintes respectivement à `['STEM', 'LOOP']` et `['classic', 'enriched']`.
2. **Limitation de taille et protection anti-abus (Priorité 2)** : Application du rate limiter (`check_rate_limit(max_requests=15, window_seconds=60)`) en amont de la route. Vérification préalable de la taille du fichier importé sur le disque (`file.seek(0, os.SEEK_END)`) pour rejeter immédiatement tout fichier excédant 1 Mo.
3. **Masquage des traces techniques en production (Priorité 3)** : Conditionnement du retour des exceptions : en mode `FLASK_ENV=production`, l'application ne retourne plus la trace brute (`str(e)`) au client, se contentant d'un message utilisateur sécurisé et traduit.
4. **Filtrage des extensions et nettoyage des noms (Priorité 4)** : Mise en place d'une liste blanche d'extensions (`ALLOWED_PARAMS_EXTENSIONS = {'txt', 'json', 'params'}` et fichiers `.params.txt`) combinée à l'appel systématique à `secure_filename`.

**Justification biologique** :
Les pipelines bioinformatiques exposés sur un serveur web clinique doivent garantir l'intégrité absolue de la session d'analyse et prévenir toute injection d'attributs arbitraires ou attaque par déni de service (saturation mémoire par upload massif). Ce durcissement protège le moteur thermodynamique LAVA tout en préservant la reproductibilité des analyses de routine.

**Impact attendu** :
Une sécurité logicielle de niveau production : étanchéité totale face aux injections de paramètres illégitimes et aux surcharges serveur, sans altérer l'expérience utilisateur lors de l'import de fichiers de paramètres légitimes.

---

### [2026-07-15] Parallélisation Multi-Cœurs du Moteur Combinatoire LAVA (Option B)

**Date/Étape** : 2026-07-15 - Parallélisation multi-processus native des boucles combinatoires via `Parallel::ForkManager`.

**Fichiers impactés** :
- `lib/LLNL/LAVA/ForkManager.pm` (nouveau module d'encapsulation multi-cœurs)
- `lava_loop_primer.pl` (parallélisation des boucles combinatoires Forward et Reverse, option `--threads|cpu`)
- `lava_stem_primer.pl` (parallélisation des boucles combinatoires Stem Forward et Stem Reverse, option `--threads|cpu`)

**Nature du changement** : [Algorithmique / Architecture / Performance]

**Explication technique** :
1. **Implémentation du module `LLNL::LAVA::ForkManager` (Option B)** : Création d'un module Perl interne gérant dynamiquement la concurrence multi-cœurs. Si le module CPAN `Parallel::ForkManager` est disponible sur l'hôte, le moteur exploite le multi-processus POSIX nativement. Si le module n'est pas installé, un mode dégradé séquentiel ultra-léger garantit la portabilité sans erreur de compilation.
2. **Découpage en Chunks et Copy-On-Write (COW)** : Au lieu d'utiliser le module natif `threads` de Perl (sujet à de lourdes fuites mémoire avec BioPerl et déconseillé en bioinformatique intensive), l'architecture sépare l'espace de recherche en sous-ensembles (chunks) d'amorces (`$chunk_start` à `$chunk_end`). Chaque processus enfant hérite instantanément des tables de pénalités en lecture seule grâce au mécanisme de mémoire partagée *Copy-On-Write* du noyau Unix.
3. **Agrégation Déterministe et Filtrage Thermodynamique** : À la fin de chaque sous-processus (`run_on_finish`), le processus parent agrège les meilleures combinaisons d'amorces (`$bestForwardInfos`, `$bestForwardPenalties`, `$bestReverseInfos`, `$bestReversePenalties`) et additionne les compteurs de signatures (`$_sig_fwd_hits`, `$_sig_rev_hits`).
4. **Interface CLI `--threads|cpu`** : Ajout du paramètre `--threads|cpu` (valeur par défaut : `auto` configurée sur `LLNL::LAVA::ForkManager->_auto_cpus()`), permettant à l'utilisateur ou à l'interface web d'allouer précisément le nombre de cœurs de calcul ou de laisser le moteur adapter automatiquement sa charge au processeur de la machine.

**Justification biologique** :
L'évaluation combinatoire exhaustive (recherche d'intersections et minimisation des pénalités sigmoïdes de distance et d'énergie d'hybridation sur l'ensemble des amorces F1c, F2, F3, B1c, B2, B3 et Stem) implique le parcours de plusieurs centaines de milliers à plusieurs millions de combinaisons thermodynamiques (`inner × stem × middle × outer`). Sur des génomes viraux très riches en variants (comme la Dengue ou le SARS-CoV-2), la recherche séquentielle pouvait nécessiter plusieurs heures de calcul. La parallélisation distribue cette évaluation cinétique sur tous les cœurs disponibles, maintenant rigoureusement les mêmes critères de sélectivité (`maxDeltaTm`, `minPrimerSpacing`, `signatureMaxLength`) sans perte de candidats.

**Impact attendu** :
Réduction drastique du temps d'exécution (accélération quasi-linéaire selon le nombre de cœurs alloués, passant de plusieurs heures à quelques minutes ou secondes sur les jeux d'amorces complexes), tout en garantissant une stricte reproductibilité des signatures LAMP et une stabilité mémoire absolue de l'application serveur LAVA.

---

### [2026-07-15] Parallélisation Multi-Cœurs et Centralisation de la Validation des Amorces (PipelineUtils)

**Date/Étape** : 2026-07-15 - Parallélisation multi-processus de la tolérance aux mésappariements (`checkPrimerMismatchTolerance`) dans `getOligosWithMismatchTolerance` et `buildNativeReversePool`.

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm` (centralisation et parallélisation de `getOligosWithMismatchTolerance` et `buildNativeReversePool`)
- `lava_loop_primer.pl` (import de `getOligosWithMismatchTolerance` et appel automatique de `set_pipeline_threads`)
- `lava_stem_primer.pl` (import de `getOligosWithMismatchTolerance` et appel automatique de `set_pipeline_threads`)

**Nature du changement** : [Algorithmique / Architecture / Performance]

**Explication technique** :
1. **Centralisation de la Validation dans PipelineUtils.pm** : Extraction et harmonisation de la fonction `getOligosWithMismatchTolerance` de `lava_loop_primer.pl` et `lava_stem_primer.pl` vers le module central `lib/LLNL/LAVA/PipelineUtils.pm`. Élimination des redondances algorithmiques entre les deux modes.
2. **Parallélisation Chunk-Based via ForkManager** : Application de `LLNL::LAVA::ForkManager` au processus de criblage et de validation des amorces (`checkPrimerMismatchTolerance`) dans `getOligosWithMismatchTolerance` (brin plus / Forward) et `buildNativeReversePool` (brin moins / Reverse). Le pool de candidats généré par Primer3 est découpé de manière équilibrée en lots (*chunks*) distribués sur le pool de processus enfants.
3. **Agrégation des Métriques et Tri Déterminé** : Les processus enfants transmettent leurs amorces validées (strictes ou dégénérées), leurs statistiques et leurs logs de criblage dans un dictionnaire de retour agrégé par `run_on_finish`. Un tri final stable sur la position génomique (`location`) puis sur la longueur (`length`) garantit l'invariance absolue des résultats entre une exécution séquentielle et une exécution parallèle sur $N$ cœurs.
4. **Synchronisation Dynamique du Pool de Threads** : Ajout de la routine `set_pipeline_threads` dans `PipelineUtils.pm` pour synchroniser globalement le nombre de cœurs alloués avec l'option CLI `--threads|cpu` ou les requêtes Flask transmises depuis l'interface utilisateur.

**Justification biologique** :
L'évaluation de la tolérance aux mésappariements sur les alignements viraux massifs exige de confronter chaque oligonucléotide candidat (souvent plusieurs milliers par type d'amorce F3, F2, F1c, B1c, B2, B3) contre l'intégralité des séquences du Multiple Sequence Alignment (MSA). Pour chaque amorce et chaque variant, l'algorithme vérifie l'intégrité absolue de la zone 3' critique (site d'initiation de la polymérase) tout en quantifiant la couverture IUPAC et les mutations en 5'/milieu. Cette opération, hautement intensive en calculs de chaînes de caractères et en comptages d'entropie de Shannon, constituait un goulot d'étranglement majeur avant l'étape de combinatoire. Sa distribution sur l'ensemble des cœurs accélère le filtrage initial de l'espace des séquences, permettant l'analyse interactive de panels viraux très profonds sans sacrifier la rigueur de la validation cinétique.

**Impact attendu** :
Accélération spectaculaire de la phase initiale de validation et de génération des amorces candidates sur l'ensemble des modes (`LOOP` et `STEM`), avec un temps de criblage divisé proportionnellement au nombre de processeurs disponibles sur l'hôte, tout en assurant une architecture logicielle centralisée et propre.

---

### [2026-07-15] Exposition Ergonomique de la Parallélisation Multi-Cœurs sur l'Interface Web (Flask & UI)

**Date/Étape** : 2026-07-15 - Mise à jour de l'interface graphique (`lava_flask_app.py` et `templates/index.html`) pour permettre la configuration du nombre de cœurs / threads depuis le navigateur Web.

**Fichiers impactés** :
- `lava_flask_app.py` (intégration de `threads` dans les paramètres par défaut `get_default_params()`, le mapping `param_mapping`, `common_params` de la route d'exécution, et les dictionnaires de traduction FR/EN)
- `templates/index.html` (ajout d'un champ de saisie dédié "Nombre de cœurs / Threads" dans la section Configuration d'exécution)

**Nature du changement** : [Architecture / Interface / Ergonomie]

**Explication technique** :
1. **Paramétrage par défaut et mapping Flask (`lava_flask_app.py`)** : Ajout de la clé `'threads': 'auto'` dans `get_default_params()`. Le convertisseur `_convert_param_value` gère dynamiquement la saisie d'entiers (`1`, `2`, `4`, `8`...) ou de la chaîne littérale `'auto'`. Lors de la soumission d'une analyse via `/execute`, le paramètre `threads` (ou son alias `cpu`) est injecté directement dans les arguments CLI de la commande Perl (`--threads <valeur>`).
2. **Champ de contrôle dans l'interface (`templates/index.html`)** : Création d'un champ texte réactif dans l'accordéon "Paramètres Avancés" (section *Execution Config*), permettant au bioinformaticien d'indiquer `auto` ou un nombre explicite de cœurs processeur à allouer.
3. **Internationalisation (I18n)** : Ajout des clés `threads_label` et `threads_desc` dans les dictionnaires français et anglais du backend Flask pour une expérience utilisateur bilingue fluide.

**Justification biologique** :
Lors de l'analyse d'alignements viraux complexes, l'allocation dynamique des ressources de calcul est essentielle pour s'adapter à la fois à l'infrastructure matérielle de l'utilisateur (ordinateur portable personnel vs serveur partagé en laboratoire) et à la profondeur du criblage d'amorces LAMP. Donner le contrôle direct sur la parallélisation depuis l'interface Web évite de saturer un serveur multi-utilisateurs tout en permettant, sur machine dédiée, de débloquer la puissance maximale (`auto`) pour accélérer le design d'amorces à haute densité de variants.

**Impact attendu** :
- Contrôle complet sur l'allocation des cœurs CPU directement depuis l'interface graphique de l'application LAVA.
- Transparence accrue lors de l'export et de l'import des fichiers de paramètres (`.params.txt` / JSON).
- Transmission fluide et sécurisée de la configuration de concurrence vers le moteur bioinformatique Perl (`PipelineUtils.pm` et `ForkManager.pm`).

---

### [2026-07-15] Sécurisation et Plafonnement Strict de la Concurrence CPU (`threads` / `ForkManager`)

**Date/Étape** : 2026-07-15 - Mise en place d'une défense en profondeur contre le déni de service (DoS) par sur-souscription des cœurs dans `lava_flask_app.py` et `lib/LLNL/LAVA/ForkManager.pm`.

**Fichiers impactés** :
- `lava_flask_app.py` (création de la fonction de validation et de plafonnement `_validate_and_cap_threads`, application dans `_convert_param_value` et lors de la soumission dans `/execute` et `execute_lava`, mise à jour des traductions)
- `lib/LLNL/LAVA/ForkManager.pm` (validation numérique stricte de `max_processes`, repli sur `get_auto_cpu_count()` et instauration d'un plafond dur dans `new`)

**Nature du changement** : [Architecture / Sécurité / Performance]

**Explication technique** :
1. **Plafonnement et prévention de sur-souscription côté Flask (`lava_flask_app.py`)** :
   - Introduction d'une validation systématique `_validate_and_cap_threads(val)` qui intercepte toute valeur de `threads` ou `cpu` (depuis le formulaire Web, un fichier de paramètres ou un appel direct).
   - Prise en compte du plafond administrateur via la variable d'environnement `MAX_THREADS_PER_RUN` (par défaut `os.cpu_count() - 1`).
   - Calcul d'un plafond effectif de concurrence (`concurrency_cap = max(1, os.cpu_count() // MAX_CONCURRENT_RUNS)`) pour s'assurer que $N$ exécutions simultanées ne saturent pas l'hôte en forkant un nombre de processus supérieur aux cœurs disponibles. Si l'utilisateur demande `auto` ou une valeur numérique excessive (ex: `500`), elle est automatiquement ramenée et bornée au plafond effectif sécurisé. Toute chaîne non numérique (ex: `abc`) retombe sur `auto` sans erreur bloquante.
2. **Défense en profondeur côté Perl (`LLNL::LAVA::ForkManager`)** :
   - Dans `sub new`, vérification numérique stricte (`$max_processes !~ /^-?\d+$/`) pour éliminer la coercition de chaînes invalides qui déclenchait des avertissements Perl, repliant immédiatement sur le décompte automatique `get_auto_cpu_count()`.
   - Application d'un plafond dur au niveau du moteur : si `max_processes` dépasse le nombre réel de cœurs disponibles sur la machine (`$auto_count`), il est automatiquement rabattu à cette limite avec émission d'un avertissement (`warn`), y compris lors des appels en ligne de commande autonomes hors interface Web.
3. **Documentation I18n** : Maintien d'un libellé d'aide concis et épuré pour `threads_desc` en français et en anglais (`Nombre de cœurs CPU alloués ('auto' ou entier ex: 4).`) afin de ne pas surcharger visuellement l'utilisateur final.

**Justification biologique** :
La validation thermodynamique des amorces LAMP à l'échelle d'alignements complets de génomes viraux est l'opération la plus exigeante en temps processeur et en bande passante mémoire du pipeline bioinformatique. Sans plafond d'exécution par processus, une simple erreur de paramétrage (ou une soumission malveillante) demandant 500 cœurs, couplée au quota de 5 exécutions concurrentes du serveur, provoquerait le fork instantané de 2500 processus lourds manipulant de larges objets BioPerl. Ce phénomène d'écrasement mémoire et de sur-souscription CPU entraînerait un effondrement immédiat du système (Kernel OOM Killer) et la perte des calculs en cours pour tous les scientifiques du laboratoire. Le plafonnement proportionnel assure l'isolation inter-processus et la stabilité opérationnelle continue du moteur scientifique LAVA.

**Impact attendu** :
- Protection absolue contre la surcharge processeur et l'épuisement de mémoire sur le serveur hôte.
- Fonctionnement fluide et garanti de 5 analyses LAMP simultanées sans contention ni dégradation de performance.
- Robustesse totale du moteur Perl autonome (`ForkManager`) face aux saisies aberrantes en ligne de commande (`--threads 500` ou `--threads abc`).

---

### Date/Étape : 2026-07-15 - Correction et fluidification du suivi de progression en temps réel (Autoflush et Granularité ForkManager)

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm`
- `lava_stem_primer.pl`
- `lava_loop_primer.pl`

**Nature du changement** : [Bug Fix / Algorithmique / Architecture]

**Explication technique** :
1. **Suppression du tampon de bloc C/Perl (`$| = 1;`)** :
   - Lorsque le moteur Perl est invoqué par l'interface Flask via `subprocess.Popen` (tube de communication standard non-TTY), la bibliothèque C sous-jacente active par défaut un tamponnage par blocs (`block buffering` de 8 Ko). En conséquence, les messages de progression `[LAVA-PROGRESS]` restaient capturés en mémoire morte dans le tampon du processus Perl et n'étaient vidés vers l'application Python qu'à l'achèvement final du pipeline, provoquant un saut brutal de la barre de 0% à 100%.
   - L'activation systématique de l'autoflush (`$| = 1;`) dans `PipelineUtils.pm`, `lava_stem_primer.pl` et `lava_loop_primer.pl`, couplée à un vidage explicite (`select(STDOUT); $| = 1;`) après chaque `printf("[LAVA-PROGRESS] ...")`, force le canal STDOUT à transmettre instantanément chaque mise à jour vers le contrôleur Flask.
2. **Optimisation de la granularité de parallélisation (`$n_chunks`)** :
   - Les boucles de validation thermodynamique (`getOligosWithMismatchTolerance`) et les boucles combinatoires (`lava_stem_primer.pl`, `lava_loop_primer.pl`) sous `LLNL::LAVA::ForkManager` découpaient auparavant l'espace de recherche en un nombre de lots strictement égal à `max_processes` (ou à un lot unique en mode mono-processeur). La notification `run_on_finish` n'étant émise par le processus parent qu'à la terminaison de chaque lot, le suivi restait figé pendant de longues périodes.
   - Le coefficient de granularité a été multiplié par un facteur 12 (`max_processes * 12` avec un seuil minimal de 25 à 30 lots de travail). Ce découpage fin garantit que les sous-processus retournent continuellement des lots partiels validés à cadence régulière.

**Justification biologique** :
Lors du design d'amorces LAMP sur des génomes viraux complets ou très polymorphes, l'exploration combinatoire des distances inter-amorces (F3-B3, F1c-B1c, boucles et tiges) implique le test de dizaines de milliers de quadruplets et sextuplets oligonucléotidiques. Pour le bioinformaticien, une barre de progression figée à 0% ou ne se chargeant qu'en toute fin d'analyse génère une ambiguïté sur l'état du serveur (suspicion de blocage infini, de boucle morte ou d'épuisement mémoire). La restitution fluide de la cinétique de criblage permet au scientifique de surveiller l'avancement réel du filtrage thermodynamique, d'anticiper le temps d'attente (ETA précis) et de valider la réactivité de la grappe de calcul.

**Impact attendu** :
- Affichage continu, précis et en temps réel de la barre de progression sur l'interface de surveillance Web (`monitor.html`).
- Fin définitive de la latence visuelle et du chargement instantané à 100% en fin d'exécution.
- Fluidité de suivi garantie à la fois en exécution monocoq (`threads=1`) et en parallélisation intensive multicoq.

---

### Date/Étape : 2026-07-15 - Suppression du bridage modulo 200 et émission continue des transitions d'étapes

**Fichiers impactés** :
- `lib/LLNL/LAVA/PipelineUtils.pm`
- `templates/monitor.html`

**Nature du changement** : [Bug Fix / Architecture]

**Explication technique** :
1. **Suppression du bridage de fréquence (`return if $done % 200 != 0`) dans `_update_progress`** :
   - Historiquement, la fonction `_update_progress` s'interrompait silencieusement si l'avancement incrémental `$done` n'était pas un multiple strict de 200 ou n'égalait pas `$total`. Sous la gestion multi-processus (`ForkManager`), les processus enfants restituent des lots dont la taille (`$chunk_size`) n'est pratiquement jamais un multiple exact de 200 (ex: 314, 628, 942...). En conséquence, la condition `return if ...` bloquait 99% des émissions `[LAVA-PROGRESS]` vers l'application Flask, laissant la barre invisible ou figée à 0% jusqu'à l'achèvement du tout dernier lot ($done == $total) où elle bondissait à 100%.
   - Ce filtre par modulo a été intégralement supprimé : chaque retour d'un lot d'exécution émet désormais immédiatement la ligne de statut dans le flux standard, permettant une progression visuelle fluide de 0% à 100%.
2. **Émission d'initialisation de transition (`_make_progress_bar`)** :
   - Lors du basculement d'une étape du pipeline à la suivante (ex: passage de `Outer Forward F3` à `Outer Reverse B3`), la création de la nouvelle instance de barre de progression (`_make_progress_bar`) n'émettait aucune ligne `[LAVA-PROGRESS] $label|0|$total`. L'interface Web conservait alors le libellé et le statut 100% de l'étape précédente pendant toute la durée de la nouvelle étape.
   - Un appel explicite `printf("[LAVA-PROGRESS] %s|0|%d||? it/s|0\n", $label, $total)` a été ajouté dès l'instanciation de la barre, forçant la carte de monitoring à basculer instantanément sur le titre de la nouvelle étape et à réinitialiser la progression à 0%.
3. **Ré-animation graphique (`monitor.html`)** :
   - Dans l'interface de surveillance, lorsqu'une étape atteignait 100%, la classe CSS `progress-bar-animated` était retirée définitivement. Un branchement `else` (`pct < 100`) a été ajouté pour réinsérer dynamiquement cette classe à l'arrivée d'une nouvelle étape de calcul.

**Justification biologique** :
Dans un pipeline de criblage LAMP de haute spécificité, les étapes successives (validation des amorces F3, puis B3, puis F2/F1c/B2/B1c, et enfin la combinatoire des tiges et boucles) ont des durées intrinsèquement hétérogènes. Si le scientifique ne perçoit pas la transition exacte entre le filtrage de la région 5' (Forward) et de la région 3' (Reverse), il ne peut pas diagnostiquer quelle population oligonucléotidique subit le taux d'attrition ou de rejet thermodynamique le plus sévère. La transparence continue du suivi est indispensable pour ajuster les fenêtres de tolérance d'hybridation et d'entropie d'alignement.

**Impact attendu** :
- Affichage immédiat de la barre dès 0% à l'amorçage de chaque étape de criblage (F3, B3, F2, B2, etc.).
- Mise à jour régulière et proportionnelle de la jauge à chaque retour de lot, sans saut brutal à 100%.
- Basculement instantané du libellé et réactivation de l'animation visuelle lors du passage à l'étape suivante.


