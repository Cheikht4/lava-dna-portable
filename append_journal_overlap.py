with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Fix Importation Paramètre Chevauchement\n")
    f.write("Fichiers impactés : lava_flask_app.py\n")
    f.write("Nature du changement : Bug Fix\n")
    f.write("Explication technique : Ajout de la clé `resolve_overlap_by` dans le dictionnaire de `get_default_params`. Le processus d'importation des paramètres se basait sur les clés de ce dictionnaire pour filtrer les paramètres valides, ce qui entraînait l'ignorance silencieuse de cette option lors de l'upload.\n")
    f.write("Justification biologique : Permettre la reproductibilité complète des exécutions, notamment la stratégie de réduction spatiale des fenêtres d'amorces (pénalité ou couverture maximale).\n")
    f.write("Impact attendu : Le paramètre `resolve_overlap_by` est désormais correctement restauré lors du chargement d'un fichier `.params.txt` précédent.\n")
