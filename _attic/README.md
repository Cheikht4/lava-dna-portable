# Répertoire d'archive technique (_attic)

Date du déplacement : 2026-10-04.

Rien dans ce répertoire n'est autoritatif : aucun fichier ne doit être exécuté, importé ni édité.
Ces fichiers sont des reliques de refactorings, prototypes et scripts de transition utilisés lors des chantiers d'évolution du moteur LAVA. Leurs fonctionnalités ont toutes été intégrées, consolidées ou dépassées dans l'historique officiel Git.

## Inventaire des composants archivés

| Nom du fichier | Date modif. | Lignes | Rôle réel | Intégration dans l'historique |
| :--- | :--- | :--- | :--- | :--- |
| `old_block_loop.pl` | 2026-08-02 21:10:59 | 332 | Section 5 (Combinaison/Validation) originale de loop extraite pour remplacement | Issu de `63dd6b0`, remplacé au commit `40f94b8` |
| `old_block_stem.pl` | 2026-08-02 21:12:37 | 357 | Section 5 (Combinaison/Validation) originale de stem extraite pour remplacement | Issu de `63dd6b0`, remplacé au commit `40f94b8` |
| `dump_block.py` | 2026-08-02 21:10:33 | 14 | Extrait les lignes de la section 5 de `lava_loop_primer.pl` vers `old_block_loop.pl` | Utilisé pour préparer le commit `40f94b8` |
| `fix_refactor.py` | 2026-08-02 21:12:22 | 11 | Supprime le filtrage strict `$inner_gap` dans `lava_loop_primer.pl` | Consolidation du commit `40f94b8` |
| `fix_stem.py` | 2026-08-02 21:14:02 | 23 | Corrige l'accès aux variables `$f_outer` et `$b_outer` dans le batching stem | Consolidation du commit `40f94b8` |
| `patch.py` | 2026-08-03 13:35:35 | 130 | Tente d'injecter la vérification de longueur `$totalLen_v` et les compteurs de rejets | Précurseur de `5bd3fc1` |
| `patch_bit_vector.py` | 2026-08-02 23:21:37 | 103 | Remplace les hachages par des bitstrings (`vec`) dans `PipelineUtils.pm` | Intégré au commit `a9ccc06` |
| `patch_exact.py` | 2026-08-03 13:36:21 | 186 | Remplacements textuels pour borner la longueur de signature | Intégré au commit `5bd3fc1` |
| `patch_inner_spacing.py` | 2026-08-03 02:19:30 | 27 | Injecte le calcul de pénalité de distance interne F1c-B1c dans la boucle loop | Intégré au commit `f2feb0c` |
| `patch_intersection_list.py` | 2026-08-02 23:42:54 | 214 | Ajoute l'option `$return_list` pour éviter la reconstruction coûteuse des identifiants | Intégré au commit `b253309` |
| `patch_length_bound.py` | 2026-08-03 13:34:33 | 118 | Première version injectant l'option `--signature_min_length` et les bornes dures | Intégré au commit `0c60c8d` |
| `patch_length_bound_safe.py` | 2026-08-03 13:35:10 | 121 | Variante sécurisée injectant les bornes dures et compteurs dans ForkManager | Intégré au commit `0c60c8d` |
| `patch_phase5.py` | 2026-08-03 00:42:18 | 277 | Prototype de refonte Phase 5 vers un partitionnement statique par worker | Précurseur de `dc6f952` |
| `patch_phase5_clean.py` | 2026-08-03 00:42:50 | 314 | Version nettoyée du patch de partitionnement statique pour `lava_loop_primer.pl` | Précurseur de `dc6f952` |
| `patch_phase5_loop.py` | 2026-08-03 00:46:04 | 277 | Injecte le Map-Reduce par partitionnement statique dans `lava_loop_primer.pl` | Intégré au commit `dc6f952` |
| `patch_phase5_stem.py` | 2026-08-03 00:52:17 | 290 | Injecte le Map-Reduce par partitionnement statique dans `lava_stem_primer.pl` | Intégré au commit `dc6f952` |
| `patch_progress.py` | 2026-08-02 22:12:04 | 75 | Ajuste l'affichage de la barre `[LAVA-PROGRESS]` (évaluées et retenues) | Intégré au commit `87b2b5a` |
| `patch_progress_init.py` | 2026-08-02 19:50:53 | 79 | Ajoute l'état initial (0 itération) de la barre de progression | Intégré au commit `63dd6b0` |
| `refactor_loop.py` | 2026-08-02 21:11:33 | 329 | Remplace la section 5 de `lava_loop_primer.pl` par la version batchée | Intégré au commit `40f94b8` |
| `refactor_stem.py` | 2026-08-02 21:13:46 | 378 | Remplace la section 5 de `lava_stem_primer.pl` par la version batchée | Intégré au commit `40f94b8` |
| `run_dengue.sh` | 2026-08-02 21:16:30 | 3 | Profilage ponctuel de charge mémoire macOS (`/usr/bin/time -l`, `Devel::Size`) | Outil ponctuel de diagnostic Chantier 6 |
| `test_bit_vector.pl` | 2026-08-02 23:20:59 | 19 | Micro-test Perl validant l'opérateur `&` et `unpack("%32b*")` | Prototypage Chantier 6 Perf |
| `test_bit_vector2.pl` | 2026-08-02 23:21:11 | 11 | Micro-test Perl vérifiant l'opérateur `&` sur vecteurs de longueurs inégales | Prototypage Chantier 6 Perf |
| `test_distance_math.pl` | 2026-08-04 12:04:14 | 18 | Brouillon de test de la formule `$outerLocation - $outerLength + 1 - $middleLocation` | Remplacé par `t/33_testDistanceMath.t` au commit `0769fed` |
| `test_primer_obj.pl` | 2026-08-04 12:04:36 | 30 | Script de test validant les coordonnées et la distance entre objets Oligo | Prototypage de validation géométrique |
| `update_journal.py` | 2026-08-02 21:14:31 | 13 | Ajoute au journal l'entrée du Chantier 6 (Batching et contrôle mémoire) | Journalisation commit `40f94b8` |
| `update_journal_bit_vector.py` | 2026-08-02 23:23:05 | 13 | Ajoute au journal l'entrée de la Phase 1 vecteurs de bits | Journalisation commit `a9ccc06` |
| `update_journal_bit_vector2.py` | 2026-08-02 23:45:01 | 13 | Ajoute au journal l'entrée du contournement de reconstruction du tableau | Journalisation commit `b253309` |
| `update_journal_perf.py` | 2026-08-02 22:47:36 | 13 | Ajoute au journal l'entrée de suppression du tri dans `calculateSignatureIntersection` | Journalisation commit `ba8ee4e` |
