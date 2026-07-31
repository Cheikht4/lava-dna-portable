# Journal d'Évolution du Projet LAVA (Version 2026)

Amélioration de la flexibilité, de la gestion de la diversité génomique et modernisation thermodynamique.

## Date/Étape : 2026-07-31 - Résolution de l'Explosion Combinatoire et Optimisation de la Recherche Top-K

**Fichiers impactés** : 
- `lava_loop_primer.pl`
- `lava_stem_primer.pl`
- `t/baseline/` (fichiers canary)

**Nature du changement** : Algorithmique / Architecture

**Explication technique** : 
- Refonte complète de la boucle d'assemblage (`Combining Best F/R Halves`). 
- Les contraintes géométriques (espacement inter-amorces) et thermodynamiques (différence de Tm) sont désormais évaluées *avant* la boucle imbriquée quadratique des combinaisons candidates, ce qui permet un élagage massif et très rapide de l'espace de recherche (hoisting de filtre).
- Par conséquent, la complexité O(K²) n'affecte plus le temps de calcul des séquences valides.
- Augmentation du paramètre `--half_signature_candidates` par défaut de 3 à 5, désormais sans impact délétère sur la performance globale (l'anomalie de temps de K=5 a été éliminée).

**Justification biologique** : 
L'exploration de combinaisons supplémentaires (K=5) augmente drastiquement la couverture des variants au sein de virus hautement variables, en donnant au moteur l'occasion d'assembler des demi-signatures qui, bien que classées un peu plus bas individuellement, forment une paire thermodynamiquement plus stable. Une meilleure tolérance aux mutations est ainsi garantie tout en respectant la cinétique d'hybridation et enzymatique de la réaction LAMP à 65°C.

**Impact attendu** : 
Le temps d'exécution redevient faible et plat quelle que soit la valeur de K (le tri et l'exploration des paires valides sont devenus marginaux dans le temps total). Les jeux viraux complexes bénéficient désormais systématiquement d'un top-K=5 par défaut pour maximiser la couverture, tout en évitant les surcoûts exponentiels injustifiés à K=10, limitant la saturation de la RAM (RSSI).
