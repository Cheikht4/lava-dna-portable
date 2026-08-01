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

### [2026-08-01] Éradication Globale du Non-Déterminisme
**Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`, `lib/LLNL/LAVA/Validator.pm`, `lib/LLNL/LAVA/PipelineUtils.pm`, `t/canary_regression.t`
**Nature du changement** : [Architecture / Bug Fix Critique]

**Explication technique** : 
1. **Tris Stables Exhaustifs** : Remplacement de 32 occurrences de tris instables (basés uniquement sur la pénalité `sort {$a->[1] <=> $b->[1]}`) par des tris complets incluant des critères géométriques et lexicographiques pour départager les ex-æquo :
   - Pour les amorces : Pénalité -> Position (`getLocation()`) -> Longueur (`getLength()`) -> Séquence (`getSequence()`).
   - Pour les signatures (LAMP) : Pénalité -> Position de départ (`getStartLocation()`) -> Longueur (`getLength()`) -> Résumé de position (`getLocationSummary()`).
2. **Itérations Ordonnées** : Remplacement de toutes les itérations non ordonnées sur les clés de hachages (`keys %hash`) introduisant un aléa de parcours (comportement par défaut depuis Perl 5.18). Les clés numériques sont désormais itérées avec `sort { $a <=> $b }` et les clés textuelles avec `sort keys`.
3. **Tests de non-régression** : Modification du sous-système de test `canary_regression.t` pour inclure la validation stricte de reproductibilité de bout en bout (y compris sur les fichiers générés comme `_amplified.fasta`), et ajout d'un double test consécutif garantissant le maintien du déterminisme pour les développements futurs.

**Justification biologique** : 
Un outil de conception d'amorces diagnostiques doit fournir un résultat strictement reproductible à partir d'un ensemble de séquences d'entrée constant. Les variations aléatoires observées précédemment dans les sorties (et le fichier `_amplified.fasta`) minaent la confiance de l'utilisateur final et entravaient la validation des processus qualités liés à l'optimisation des amorces, empêchant l'identification formelle des candidats identiques d'un lancement à l'autre.

**Impact attendu** : 
- **Reproductibilité totale** : Deux exécutions identiques produiront désormais toujours un hash SHA-256 de sortie identique, quel que soit le nombre de cœurs alloués (`--threads 1` ou `--threads 8`) ou l'environnement système.
- Régénération nécessaire des baselines en raison des variations de départage introduites.
