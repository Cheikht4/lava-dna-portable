with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Paramètre Threads\n")
    f.write("Fichiers impactés : lava_flask_app.py\n")
    f.write("Nature du changement : Bug Fix\n")
    f.write("Explication technique : Remise du paramètre `threads` sur `'auto'` par défaut au lieu de la valeur `9` extraite du run précédent.\n")
    f.write("Justification biologique : Laisser le système décider automatiquement du nombre de cœurs à utiliser en fonction du matériel serveur pour éviter de saturer inutilement la machine ou d'imposer une limite en dur inadaptée à tous les hôtes.\n")
    f.write("Impact attendu : Meilleure portabilité du moteur LAVA sans affecter la reproductibilité du comportement algorithmique.\n")
