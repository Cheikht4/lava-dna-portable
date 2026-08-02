with open("LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("""
### [2026-07-30] Implémentation Top-K (Sélection par Couverture)
**Fichiers impactés** : `lava_loop_primer.pl`, `lava_stem_primer.pl`, `lib/LLNL/LAVA/PipelineUtils.pm`
**Nature du changement** : [Algorithmique / Optimisation de Couverture]

**Explication technique** : 
Le comportement historique consistait à ne retenir qu'un seul "meilleur" candidat thermodynamique (`$bestSetPenalty`) pour chaque paire d'amorces (ex: F1c/F2/F3) pendant la phase de Branch & Bound. 
Nous avons introduit une logique **Top-K** : 
1. `PipelineUtils.pm` extrait les `compatible_sequence_ids` et les compacte dans un bit-vector de 32 bits (via `pack`) pour chaque amorce lors du pré-traitement (`flattenInfoData`).
2. Dans les boucles Forward et Reverse, l'algorithme intercepte tous les candidats valides et calcule leur **couverture globale** par intersection binaire rapide (`&`) des bit-vectors.
3. Un tableau `@topCandidates` conserve les $K$ meilleures combinaisons (défaut `half_signature_candidates = 5`) en les triant en priorité par leur couverture (décroissante), puis par leur pénalité thermodynamique (croissante).
4. La phase de combinaison finale assemble ensuite de manière exhaustive ces $K$ meilleures demi-signatures.

**Justification biologique** : 
Privilégier aveuglément l'amorce thermodynamiquement "parfaite" pouvait conduire à l'élimination prématurée d'une excellente candidate couvrant 95% des variants viraux, au profit d'une amorce couvrant 50% des variants mais ayant un $T_m$ très légèrement meilleur. La logique Top-K garantit la survie des amorces à forte couverture virale jusqu'à l'étape d'assemblage final, sans faire exploser l'espace combinatoire (puisque limité à 5 par racine).

**Impact attendu** : 
- Capacité inédite du logiciel à "rattraper" des signatures très couvrantes sur des virus hautement variables.
- Aucune perte de vitesse notable grâce à l'implémentation par opérations binaires natives (`unpack("%32b*")`).
- Paramétrage fin accessible via `--half_signature_candidates`.
""")
