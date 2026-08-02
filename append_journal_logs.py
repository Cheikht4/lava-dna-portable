with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Purge des logs de débogage (Penalty Guard & B&B)\n")
    f.write("Fichiers impactés : lava_loop_primer.pl\n")
    f.write("Nature du changement : Bug Fix / Optimisation\n")
    f.write("Explication technique : Mise en commentaire des lignes d'avertissement `[PENALTY GUARD]` ainsi que des statistiques d'élagage (`[Forward B&B] Elagage` et `[Reverse B&B] Elagage`) qui venaient polluer la sortie standard.\n")
    f.write("Justification biologique : Ces messages, très fréquents lors de la conception d'amorces d'enrichissement LOOP ou STEM sur des régions denses (où la distance loopToMiddle peut être négative), étaient purement diagnostiques. L'utilisateur final n'a besoin de voir que la progression principale du pipeline, sans être submergé par les avertissements internes de l'arbre thermodynamique de Branch & Bound.\n")
    f.write("Impact attendu : Une sortie standard et une interface de logs beaucoup plus propres et lisibles, évitant la confusion avec de véritables erreurs.\n")
