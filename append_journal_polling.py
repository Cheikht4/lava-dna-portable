with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Optimisation Fréquence de Polling\n")
    f.write("Fichiers impactés : templates/monitor.html\n")
    f.write("Nature du changement : Architecture / Optimisation\n")
    f.write("Explication technique : Remplacement de l'intervalle statique (`setInterval` à 2000ms) par un délai dynamique (`setTimeout`) permettant de moduler la fréquence de rafraîchissement selon l'état du processus : 10 secondes si le statut est `queued`, et 2 secondes s'il est `starting` ou `running`.\n")
    f.write("Justification biologique : L'optimisation des requêtes asynchrones allège considérablement la charge sur le serveur web (Gunicorn/Flask). En effet, interroger l'état d'un processus en file d'attente (qui peut le rester de nombreuses minutes selon la saturation des ressources CPU pour les alignments génomiques complexes) toutes les 2 secondes générait un trafic inutile. Une récurrence de 10 secondes est amplement suffisante pour suivre la position dans la file.\n")
    f.write("Impact attendu : Réduction de la pression sur le serveur pour les requêtes de monitoring en période de forte affluence, tout en maintenant une excellente réactivité (2s) lors de l'exécution proprement dite.\n")
