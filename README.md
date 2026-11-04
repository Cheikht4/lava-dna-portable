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

## Vérification de l'installation et Utilisation CLI / Command-Line Usage & Presets

LAVA fournit deux scripts de conception indépendants :
- `lava_loop_primer.pl` : conçoit des signatures LAMP avec amorces de boucle (`F3, B3, FIP [F1c+F2], BIP [B1c+B2], FLOOP, BLOOP`) ou en mode classique à 6 amorces (`--include_loop_primers 0`).
- `lava_stem_primer.pl` : conçoit des signatures LAMP avec amorces de tige (`F3, B3, FIP [F1c+F2], BIP [B1c+B2], FSTEM, BSTEM`) ou en mode classique à 6 amorces (`--include_stem_primers 0`).

### 1. Commandes de démarrage rapide (Quick Start)

Vous pouvez vérifier que le moteur scientifique fonctionne correctement sur l'alignement de référence fourni (`t/fixtures/rota_canary_40.fasta`, 40 séquences de Rotavirus A) :

```bash
# 1. Conception loop-LAMP (8 amorces : F3, B3, F2, B2, F1, B1, FLOOP, BLOOP)
perl -Ilib lava_loop_primer.pl --option_file t_data/normal_parameters.xml \
  --alignment_fasta t/fixtures/rota_canary_40.fasta --output_file /tmp/out_loop.fasta

# 2. Conception stem-LAMP (8 amorces : F3, B3, F2, B2, F1, B1, FSTEM, BSTEM)
perl -Ilib lava_stem_primer.pl --option_file t_data/normal_stem_parameters.xml \
  --alignment_fasta t/fixtures/rota_canary_40.fasta --output_file /tmp/out_stem.fasta

# 3. Conception LAMP classique à 6 amorces (sans amorces de boucle : F3, B3, F2, B2, F1, B1)
perl -Ilib lava_loop_primer.pl --option_file t_data/normal_parameters.xml \
  --alignment_fasta t/fixtures/rota_canary_40.fasta --output_file /tmp/out_classic.fasta \
  --include_loop_primers 0
```

### 2. Fichiers de préréglages (`t_data/`) et priorité des options

- **Priorité de la ligne de commande** : toute option passée explicitement sur la ligne de commande l'emporte sur la balise correspondante du fichier `--option_file` (par exemple `--alignment_fasta`, `--output_file` ou `--include_loop_primers 0`).
- **Gradation des préréglages loop-LAMP** sur `t/fixtures/rota_canary_40.fasta` (40 séquences) :
  - `t_data/strict_parameters.xml` : **1 121** combinaisons candidates évaluées -> **1** signature finale retenue (`37/40` séquences, couverture `92,50 %`, pénalité `7,90`).
  - `t_data/normal_parameters.xml` : **5 418** combinaisons candidates évaluées -> **1** signature finale retenue (`37/40` séquences, couverture `92,50 %`, pénalité `7,90`).
  - `t_data/loose_parameters.xml` : **17 333** combinaisons candidates évaluées -> **1** signature finale retenue (`37/40` séquences, couverture `92,50 %`, pénalité `7,90`).
  - En mode classique à 6 amorces (`--option_file t_data/normal_parameters.xml --include_loop_primers 0`) : **27 433** combinaisons candidates évaluées -> **1** signature finale retenue (`37/40` séquences, couverture `92,50 %`, pénalité `3,70`).
- **Note sur `t_data/normal_stem_parameters.xml`** : la fixture `t/fixtures/rota_canary_40.fasta` n'offre que **205 nucléotides exploitables**, ce qui est très court pour empiler 8 amorces sans chevauchement (`F3, F2, F1, FSTEM, BSTEM, B1, B2, B3`). Le fichier `t_data/normal_stem_parameters.xml` utilise donc des longueurs d'amorces stem réduites (`15/18/22 nt`) adaptées à cette fenêtre courte (**417** candidats évalués -> **1** signature finale, couverture `92,50 %`, pénalité `7,89`) ; il constitue un point de départ pour vérifier l'installation plutôt qu'un préréglage de production universel.

Assurez-vous d'avoir installé les modules Perl requis (`Bio::SeqIO`, `Bio::AlignIO`, `Bio::SimpleAlign`, `Bio::LocatableSeq`, `XML::LibXML`, ainsi que `Clone` via `Bio::Root::Root`) et que `primer3_core` (>= 2.5) est bien dans le `PATH`.
