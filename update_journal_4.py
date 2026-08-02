import sys

with open("LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("""
### Date/Étape : 2026-08-02 - Amélioration UX de la phase de Combinaison
- **Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`
- **Nature du changement** : Architecture / UX
- **Explication technique** : 
  1. Ajout d'une barre de progression interactive `[LAVA-PROGRESS]` dans la boucle `for` principale de la phase `Combining Best F/R Halves to create LAMP Signatures...` pour les deux scripts (Loop et Stem).
  2. Le suivi affiche désormais le nombre de combinaisons évaluées, le total, le nombre de signatures créées, la vitesse (itérations par seconde) et le temps restant estimé (ETA).
- **Justification biologique** : 
  - Lors de la combinaison des moitiés Forward et Reverse, l'espace des possibles peut être extrêmement vaste (plusieurs dizaines de milliers de paires potentielles), surtout sur des alignements permissifs. Auparavant, l'utilisateur restait face à un terminal figé sans savoir si le programme cherchait activement ou était bloqué. Cette barre de progression rend cette étape de calcul combinatoire transparente.
- **Impact attendu** :
  - Meilleure expérience utilisateur avec un retour visuel sur l'avancement de la création des signatures avant de passer à l'étape de validation.
""")
