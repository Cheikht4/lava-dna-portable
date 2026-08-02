with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write("\nDate/Étape : 2026-07-26 / Fix Auto-Refresh File d'Attente\n")
    f.write("Fichiers impactés : templates/monitor.html\n")
    f.write("Nature du changement : Bug Fix\n")
    f.write("Explication technique : Ajout du statut `queued` dans les conditions de déclenchement de l'intervalle JavaScript (`setInterval`) responsable de l'actualisation asynchrone de la page de monitoring.\n")
    f.write("Justification biologique : Suite à l'introduction récente du système de file d'attente (FIFO) pour la gestion de l'afflux des tâches et la limite CPU, les nouveaux jobs commençaient souvent par l'état `queued`. Or, le script de la page ne s'actualisait automatiquement que pour `starting` ou `running`, forçant l'utilisateur à rafraîchir manuellement pour voir quand son job démarrait réellement.\n")
    f.write("Impact attendu : Le monitoring s'actualise désormais tout seul dès la soumission du job, affichant en temps réel le passage de `En attente` à `En cours` sans aucune action de l'utilisateur.\n")
