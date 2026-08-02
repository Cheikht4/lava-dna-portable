with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Fix CSRF Import Paramètres\n")
    f.write("Fichiers impactés : templates/index.html\n")
    f.write("Nature du changement : Bug Fix\n")
    f.write("Explication technique : Ajout du token CSRF (`{{ csrf_token() }}`) à l'objet `FormData` dans la fonction JavaScript `uploadParamsFile` qui utilise `fetch`.\n")
    f.write("Justification biologique : Correction d'un défaut d'expérience utilisateur suite au renforcement de sécurité.\n")
    f.write("Impact attendu : Les utilisateurs peuvent à nouveau importer leurs fichiers de paramètres via le bouton d'import.\n")
