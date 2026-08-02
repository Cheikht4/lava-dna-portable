import datetime

entry = f"""
Date/Étape : {datetime.date.today().strftime('%Y-%m-%d')} / Bug Fix - Écrasement silencieux des bornes de Tm
Fichiers impactés : lava_loop_primer.pl, lava_stem_primer.pl
Nature du changement : Bug Fix
Explication technique : L'affectation des contraintes `min_tm` et `max_tm` écrasait silencieusement la plage demandée par l'utilisateur si la valeur cible (`target_tm`) s'y trouvait en-dehors. Une fonction `clamp_tm_target` a été factorisée pour corriger les 16 occurrences (4 types d'amorces × 2 bornes × 2 scripts). Dorénavant, c'est le Tm cible qui est recadré pour rentrer dans les bornes spécifiées par l'utilisateur, et ce recadrage déclenche l'émission d'un message d'avertissement explicite (`AVERTISSEMENT : --[param]_target_tm ... hors de la plage demandée ...`).
Justification biologique : Les températures de fusion (Tm) dictent la cinétique d'hybridation à 65°C. Les biologistes définissent souvent des plages précises (ex: 62°C - 65°C) pour sélectionner uniquement les candidats les plus stables thermodynamiquement, typiquement pour les amorces externes. L'écrasement silencieux ramenait la borne à 58°C, générant de mauvaises amorces sans le notifier. Le recadrage du target respecte la contrainte tout en satisfaisant les garde-fous de Primer3.
Impact attendu : Le pipeline avertira l'utilisateur s'il détecte une incohérence et ne sacrifiera plus jamais les bornes de Tm. Les tests "Cas 1" du diagnostic (Pools vides) peuvent désormais être déclenchés par de fortes contraintes Tm.
"""

with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write(entry)

print("Journal updated.")
