# LAVA-DNA Web Interface & Production Deployment Package

## Description (Français)

Ce dépôt constitue le package officiel autonome de l'interface graphique web (Flask / Gunicorn / Nginx) et du moteur scientifique LAVA pour le design avancé d'amorces LAMP (Loop-Mediated Isothermal Amplification). Il est optimisé pour un déploiement clé en main sur serveur de production ou une exécution locale simple par les chercheurs.

### Fonctionnalités principales
- **Interface bilingue complète (FR / EN)** : Bascule instantanée de la langue pour la configuration, la surveillance en temps réel (barre de progression LAVA-Progress) et les diagnostics d'échec.
- **Support des virus variables et références uniques** : Calcul combinatoire proportionnel à la taille de la région cible avec tolérance aux codes IUPAC et validation d'alignement FASTA.
- **Sécurité et Robustesse** : Cookies de session sécurisés (`HttpOnly`, `SameSite=Lax`), contrôle d'intégrité des fichiers téléchargés et gestion d'arrière-plan résiliente.
- **Package de Déploiement Production** : Fichiers de configuration Nginx (`deployment/nginx_lava.conf`), service systemd (`deployment/lava-dna.service`) et Gunicorn (`deployment/gunicorn_config.py`).

---

## Description (English)

This repository provides the standalone web interface package (Flask / Gunicorn / Nginx) and LAVA scientific engine for advanced LAMP primer design. It is streamlined for turnkey production server deployment or simple local execution.

### Key Features
- **Full Bilingual Interface (FR / EN)** : Instant language switching across parameter configuration, live monitoring (LAVA-Progress bar), and technical error diagnosis.
- **Support for Highly Variable Viruses & Single References** : Proportional combinatorial search with unrestricted IUPAC degeneracy support and strict FASTA alignment validation.
- **Security & Robustness** : Hardened session cookies (`HttpOnly`, `SameSite=Lax`), upload validation, and resilient background task management.
- **Production Deployment Suite** : Includes Nginx configuration (`deployment/nginx_lava.conf`), systemd service files (`deployment/lava-dna.service`), and Gunicorn configuration (`deployment/gunicorn_config.py`).

---

## Guide d'Installation Rapide / Quick Start Guide

### 1. Installation locale (Développement / Utilisation directe)

```bash
# 1. Créer et activer l'environnement virtuel Python
python3 -m venv lava_env
source lava_env/bin/activate

# 2. Installer les dépendances Python
pip install -r requirements_flask.txt

# 3. Installer les dépendances Perl / Conda (si non présentes sur le système)
# Assurez-vous d'avoir Perl, Primer3 et BioPerl installés via environment.yml

# 4. Lancer l'application Flask
python3 lava_flask_app.py
```
Ouvrez ensuite votre navigateur sur : `http://127.0.0.1:5000`

---

### 2. Déploiement en Production (Linux / Nginx / Gunicorn)

Consultez le dossier `deployment/` qui contient les scripts automatisés :
```bash
# Exécution du script de déploiement (nécessite les privilèges root/sudo)
cd deployment/
chmod +x deploy.sh
sudo ./deploy.sh
```

**⚠️ Étape de sécurité critique (Service Public)**
Pour sécuriser les cookies de session et activer le CSRF, vous DEVEZ générer une clé secrète statique.
Créez le fichier `/etc/lava-dna/lava.env` :
```bash
sudo mkdir -p /etc/lava-dna
echo "SECRET_KEY=$(python3 -c 'import secrets; print(secrets.token_hex(24))')" | sudo tee /etc/lava-dna/lava.env
sudo chown lavauser:lavauser /etc/lava-dna/lava.env
sudo chmod 600 /etc/lava-dna/lava.env
sudo systemctl restart lava-dna
```

> **Note :** Tout redémarrage de Gunicorn (`systemctl restart lava-dna`) purge intégralement la file d'attente en cours, car celle-ci est gérée de manière volatile en mémoire pour des raisons de performance.

---
## Architecture du Répertoire

- `lava_flask_app.py` : Contrôleur principal Flask et gestionnaire de tâches asynchrones.
- `launch_lava_smart_kill.py` : Gestionnaire d'arrêt propre des processus en arrière-plan.
- `templates/` : Pages HTML bilingues (`index.html`, `monitor.html`, `executions.html`).
- `static/` : Styles CSS et assets graphiques.
- `deployment/` : Scripts et fichiers de configuration pour la production (Nginx/systemd).
- `lava_loop_primer.pl` / `lava_stem_primer.pl` : Moteurs scientifiques Perl.
- `lib/` : Modules algorithmiques et thermodynamiques Perl LAVA.

---
## Licence & Droits d'Utilisation / License & Terms of Use

* **Moteur scientifique Perl (`lava_loop_primer.pl`, `lava_stem_primer.pl`, modules hérités)** : sous licence open-source BSD 3-Clause (LLNL / Clinton Torres / Cheikh Talibouya).
* **Interface Graphique Web & Suite de Déploiement (`lava_flask_app.py`, `templates/`, `static/`, `deployment/`)** : **Licence Propriétaire - Tous droits réservés (Cheikh Talibouya)**. L'utilisation, la reproduction, la modification, la distribution ou le déploiement clinique/commercial de cette interface web est strictement soumis à l'autorisation écrite préalable de l'auteur. Voir le fichier `LICENSE` pour les détails complets.

## Vérification de l'installation

Vous pouvez vérifier que le moteur scientifique fonctionne correctement en lançant un test de fumée avec les données de test fournies (dossier `t_data/`) :

```bash
# Vérifier la conception d'amorces STEM sur un gène cible (S. aureus)
perl lava_stem_primer.pl --params t_data/s_aureus_parameters.xml
```

Si la commande réussit, vous verrez des fichiers de sortie générés dans le dossier courant (ex: `results_stem_*.primers`). Assurez-vous d'avoir installé les modules Perl `BioPerl` et `XML::LibXML` et que `primer3_core` est bien dans le PATH.
