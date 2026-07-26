import os

with open('LAVA_EVOLUTION_JOURNAL.md', 'a') as f:
    f.write("\n### Date/Étape : 2026-07-25 - Correction des amorces fixées sur brin moins et ajout de l'option d'optimisation\n")
    f.write("- **Fichiers impactés** : `lib/LLNL/LAVA/PipelineUtils.pm`, `lava_stem_primer.pl`, `lava_loop_primer.pl`\n")
    f.write("- **Nature du changement** : [Bug Fix / Algorithmique]\n")
    f.write("- **Explication technique** : \n")
    f.write("  1. Correction géométrique du champ `location` pour les amorces fixées sur le brin moins. La valeur `$position` (bord gauche) a été remplacée par `$position + length($primer_seq) - 1` (bord droit), pour respecter la convention absolue du moteur (cf. `buildNativeReversePool`).\n")
    f.write("  2. Ajout de l'option `--fixed_primer_optimize` (0 ou 1, par défaut 1). Lorsqu'elle est désactivée (0), le pipeline calcule tout de même la couverture via `checkPrimerMismatchTolerance` pour assurer l'exactitude de `compatible_seq_ids`, mais ignore la séquence dégénérée produite par le B&B, conservant l'amorce originale inchangée et forçant son inclusion.\n")
    f.write("- **Justification biologique** : Une amorce sur le brin antisens était décalée dans la fenêtre génomique, faussant l'évaluation de sa couverture sur l'alignement et provoquant une dégénérescence excessive inutile par le B&B. L'option d'optimisation offre aux utilisateurs la garantie de conserver des séquences d'amorces certifiées en paillasse sans altération algorithmique, tout en bénéficiant du rapport de couverture.\n")
    f.write("- **Impact attendu** : Rendu fidèle des coordonnées des amorces Reverse fixées, couverture précise de ces amorces, et possibilité d'injecter des amorces validées empiriquement sans risque de dégénérescence non désirée.\n")
