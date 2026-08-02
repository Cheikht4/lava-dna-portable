import sys

with open("LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("""
### Date/Étape : 2026-08-02 - Refonte de la Journalisation de Validation
- **Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`, `lib/LLNL/LAVA/PipelineUtils.pm`
- **Nature du changement** : Architecture / Bug Fix
- **Explication technique** : 
  1. Remplacement des `print` STDOUT exhaustifs pour chaque signature évaluée par une barre de progression en temps réel `[LAVA-PROGRESS]` dans le processus parent `ForkManager`.
  2. Implémentation d'un bloc de résumé statistique de la validation (Validées, Rejetées, Couverture MAX des rejetées, Distribution).
  3. Ajout de l'option `--verbose_validation` qui redirige les logs détaillés de `calculateSignatureIntersection` directement vers un fichier compressé à la volée (`<output_base>_validation_detail.log.gz`), évitant toute saturation de la mémoire du terminal.
- **Justification biologique** : 
  - La validation des génomes hautement variables (ex: Dengue) générait plus de 31 millions de lignes de logs vers STDOUT. Bien que le moteur Perl gère parfaitement la mémoire des données, le tampon du terminal client saturait la RAM système (plus de 94 Go consommés par Terminal.app). Cette refonte limite drastiquement le trafic I/O texte sans perte d'information utile. Les données détaillées restent disponibles à la demande pour l'analyse des cas marginaux via `--verbose_validation`.
- **Impact attendu** :
  - Disparition totale des plantages du terminal et des logs monstrueux impossibles à ouvrir. La console n'affiche plus que la progression propre et un résumé analytique actionnable.
""")
