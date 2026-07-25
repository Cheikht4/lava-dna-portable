import os

with open('LAVA_EVOLUTION_JOURNAL.md', 'a') as f:
    f.write("\n### Date/Étape : 2026-07-25 - Correction finale des compteurs de progression\n")
    f.write("- **Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`, `apply_optimizations.py`\n")
    f.write("- **Nature du changement** : [Bug Fix]\n")
    f.write("- **Explication technique** : Le système d'affichage de la progression LAVA-PROGRESS a été restructuré. Au lieu de cumuler le statut de chaque chunk via `flock` sur un fichier unique (ce qui provoquait des verrous lents et une lecture asymétrique aboutissant à des compteurs incohérents), le script crée désormais un répertoire de progression avec un fichier `.prog` distinct par chunk (`chunk_$chunk_id.prog`). La boucle de lecture utilise `glob` pour lire dynamiquement l'état agrégé.\n")
    f.write("- **Justification biologique** : Les affichages incohérents de l'avancement induisent en erreur l'utilisateur lors de grands runs combinatoires. Ce correctif redonne un indicateur fiable et proportionné sans affecter les performances.\n")
    f.write("- **Impact attendu** : Affichage correct et strictement borné à 100% dans l'interface de chargement.\n")
