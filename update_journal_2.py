import sys

with open("LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("""
### Date/Étape : 2026-08-02 - Correctif Mémoire sur la Parallélisation de Validation
- **Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`
- **Nature du changement** : Architecture / Bug Fix
- **Explication technique** : 
  1. Correction d'un problème de croissance exponentielle de la mémoire dans les enfants `ForkManager`. Au lieu d'accumuler tous les résultats volumineux (`signature_intersection_ids` et `primer_coverage_details`) pour l'ensemble d'un thread, les processus enfants n'émettent plus que 4 valeurs scalaires : `[$idx, $coverage, $status, $target_count]`.
  2. Le `ForkManager` gère désormais les processus enfant **par chunk**, garantissant qu'un processus de validation est tué et libère sa mémoire toutes les X signatures.
  3. L'intersection détaillée complète (incluant la liste de toutes les séquences cibles pour le fasta généré) est recalculée dynamiquement, en mono-thread, UNIQUEMENT sur la liste restreinte de signatures retenues par la fonction de réduction de chevauchement (`reduceSignaturesByOverlap`).
  4. Ajout explicite de l'instruction `use POSIX;` pour sécuriser l'appel à `POSIX::ceil`.
- **Justification biologique** : 
  - La sérialisation d'objets extrêmement lourds tels que la cartographie complète de l'intersection de toutes les amorces avec toutes les séquences Dengue cibles pour chaque signature engorgeait les processus (1.2 million de signatures * 100 séquences = saturation de la RAM et surcharge CPU IPC). Ce correctif sépare le screening rapide de la couverture (multi-thread) de l'extraction des données métaboliques de l'intersection (mono-thread post-filtrage).
- **Impact attendu** :
  - Disparition totale des plantages par dépassement de mémoire (OOM) lors des runs sur les virus à large spectre (Dengue). Les 20 tests Canary et de Déterminisme confirment la préservation de l'exactitude stricte de la sortie.
""")
