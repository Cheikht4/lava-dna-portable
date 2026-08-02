import datetime

entry = f"""
Date/Étape : {datetime.date.today().strftime('%Y-%m-%d')} / Bug Fix - Portée de $options_r dans le bloc de diagnostic
Fichiers impactés : lava_loop_primer.pl, lava_stem_primer.pl
Nature du changement : Bug Fix
Explication technique : Le bloc de diagnostic accédait à `$options_r` via `$main::options_r` qui est inexistant puisque la variable est définie localement (lexicale). J'ai passé `$options_r` explicitement en argument à `print_zero_signature_diagnostic`. J'ai également ajusté les écarts minimaux de Tm pour les afficher avec 2 décimales via `sprintf("%.2f", ...)`.
Justification biologique : L'utilisation de `$main::options_r` entraînait le repli vers les valeurs par défaut dans le rapport pour les paramètres (100% couverture, 0 dégénérées), donnant de fausses informations à l'utilisateur lors du "Cas 1" (Pools vides). Un rapport contenant les vrais paramètres permet un diagnostic fiable pour que le biologiste relâche les contraintes ciblées. De plus, le format 2 décimales pour l'écart de Tm lève la contradiction visuelle entre une limite de 0.1°C atteinte à 0.14°C, favorisant l'analyse thermodynamique fine.
Impact attendu : Le fichier de rapport est maintenant créé exactement là où il le faut (`<output_file>_diagnostic.txt`), et contient les vrais paramètres du run. Les écarts Tm sont clairement distincts du seuil.
"""

with open("/Users/cheikhtalibouya/Documents/lava/LAVA_EVOLUTION_JOURNAL.md", "a") as f:
    f.write(entry)

print("Journal updated.")
