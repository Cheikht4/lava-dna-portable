import datetime

entry = f"""
Date/Étape : {datetime.date.today().strftime('%Y-%m-%d')} / Implémentation du diagnostic explicite d'échec
Fichiers impactés : lava_loop_primer.pl, lava_stem_primer.pl
Nature du changement : Algorithmique
Explication technique : Ajout de compteurs de rejets (`rej_geometry`, `rej_spacing`, `rej_loopgap`, `rej_tm_*`) et d'enregistreurs d'écarts minimaux (`min_delta_tm_*`) au sein des boucles combinatoires Forward et Reverse. Ces métriques, remontées sans surcoût par `ForkManager` via `run_on_finish`, sont analysées lors d'un échec (0 signature) pour générer un rapport (`_diagnostic.txt`) identifiant la contrainte limitante (Cas 1 : amont vide, Cas 2 : problème d'assemblage).
Justification biologique : Un échec silencieux ("No valid forward primer combinations found") n'aide pas le biologiste. En pointant avec précision quelle contrainte thermodynamique (ex: delta Tm trop restrictif entre Inner et Loop) ou géométrique (ex: empan minimal requis supérieur au max autorisé) cause le rejet de toutes les amorces, l'utilisateur peut ajuster intelligemment les paramètres de LAVA pour cibler le virus ou la souche en question.
Impact attendu : Les utilisateurs recevront un diagnostic clair lors d'un run échoué, leur suggérant des modifications paramétriques concrètes (ex: ajuster `--max_tm_diff` ou `--signature_max_length`) pour récupérer des signatures, sans impacter les performances de l'outil.
"""

with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write(entry)

print("Journal updated.")
