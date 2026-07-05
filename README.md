# 🧬 LAVA-DNA Interface

> **Interface graphique moderne pour la conception de primers LAMP avec LAVA**

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://python.org)
[![Flask](https://img.shields.io/badge/Flask-3.0+-green.svg)](https://flask.palletsprojects.com)
[![Perl](https://img.shields.io/badge/Perl-5.26+-orange.svg)](https://perl.org)
[![Primer3](https://img.shields.io/badge/Primer3-2.6.1-purple.svg)](https://primer3.org)
[![License](https://img.shields.io/badge/License-BSD--3--Clause_%2F_Proprietary-blue.svg)](LICENSE)

## 🎯 **Vue d'ensemble**

LAVA-DNA Interface est une interface web moderne et intuitive pour le logiciel LAVA (LLNL LAMP Assay Validation and Analysis), permettant la conception automatisée de primers LAMP (Loop-Mediated Isothermal Amplification) pour la détection d'ADN.

### ✨ **Fonctionnalités principales**

- 🧬 **Conception de primers STEM et LOOP** LAMP
- 🎛️ **Interface web intuitive** avec Flask
- 📊 **Monitoring en temps réel** des exécutions
- ⚙️ **Configuration avancée** de tous les paramètres LAVA
- 📁 **Gestion des fichiers** FASTA et résultats
- 🚀 **Lancement automatique** du navigateur
- 🛡️ **Protection anti-faux-positifs** pour la fermeture

## 🚀 **Installation rapide**

### **Prérequis système**
- Ubuntu 20.04+ / Debian 11+ / WSL2 / macOS (via Homebrew)
- Python 3.9+
- Perl 5.26+
- Primer3 **2.6.1** (version testée et validée)
- BioPerl 1.6.924+
- Flask 3.0.3+ (interface web uniquement)

### **Installation manuelle**

#### **Pour Linux / WSL2**
```bash
# Dépendances système
sudo apt update
sudo apt install -y python3 python3-pip python3-venv primer3 bioperl libbioperl-perl cpanminus

# Module Perl
sudo cpanm Bio::Tools::Run::Primer3

# Environnement Python
rm -rf lava_env
python3 -m venv lava_env
source lava_env/bin/activate
pip install -r requirements.txt
pip install -r requirements_flask.txt

# Configuration Perl
export PERL5LIB=$(pwd)/lib:$PERL5LIB
echo "export PERL5LIB=$(pwd)/lib:\$PERL5LIB" >> ~/.bashrc
```

#### **Pour macOS**
```bash
# Dépendances système (nécessite Homebrew)
brew update
brew install python primer3 perl cpanminus

# Module Perl
cpanm Bio::Tools::Run::Primer3

# Environnement Python
rm -rf lava_env
python3 -m venv lava_env
source lava_env/bin/activate
pip install -r requirements.txt
pip install -r requirements_flask.txt

# Configuration Perl
export PERL5LIB=$(pwd)/lib:$PERL5LIB
echo "export PERL5LIB=$(pwd)/lib:\$PERL5LIB" >> ~/.zshrc
```

## 🎮 **Utilisation**

### **Lancement de l'interface**

**🎯 5 SCRIPTS DE LANCEMENT DISPONIBLES :**

#### **1️⃣ Script RECOMMANDÉ - Anti-faux-positifs**
```bash
# Protection complète contre les fermetures accidentelles
python launch_lava_no_false_positive.py
```
**✅ Idéal pour :** Utilisation intensive, tests multiples, navigation complexe

#### **2️⃣ Script intelligent - Fermeture contrôlée**
```bash
# Fermeture automatique après inactivité (plus robuste)
python launch_lava_smart_kill.py
```
**✅ Idéal pour :** Utilisation normale avec fermeture automatique

#### **3️⃣ Script basique - Fermeture simple**
```bash
# Fermeture automatique basique
python launch_lava_browser_kill.py
```
**✅ Idéal pour :** Utilisation basique

#### **4️⃣ Script WSL2 - Ouverture automatique**
```bash
# Spécialement conçu pour WSL2
python launch_lava_wsl.py
```
**✅ Idéal pour :** Environnement WSL2, ouverture auto navigateur Windows

#### **5️⃣ Script manuel - Contrôle total**
```bash
# Lancement manuel sans fermeture automatique
python launch_lava.py
# Puis ouvrir manuellement http://localhost:5000
```
**✅ Idéal pour :** Développement, tests, utilisation avancée

### **Workflow typique**
1. **📁 Upload** d'un fichier FASTA d'alignement
2. **⚙️ Configuration** des paramètres (STEM ou LOOP)
3. **🚀 Lancement** de l'exécution LAVA
4. **📊 Monitoring** en temps réel
5. **📥 Téléchargement** des résultats

## 📁 **Structure du projet**

```
lava-dna-interface/
├── 📁 lib/                          # Modules Perl LAVA
├── 📄 lava_stem_primer.pl           # Script STEM-LAMP
├── 📄 lava_loop_primer.pl           # Script LOOP-LAMP
├── 📄 lava_flask_app.py             # Application Flask principale
├── 📁 templates/                     # Templates HTML
├── 📄 requirements_flask.txt         # Dépendances Python
├── 📄 DOCUMENTATION_LAVA.txt         # Documentation complète
├── 📄 GUIDE_INSTALLATION_NOUVELLE_MACHINE.md  # Guide d'installation
└── 🚀 SCRIPTS DE LANCEMENT (5 versions) :
    ├── 📱 launch_lava_no_false_positive.py  # RECOMMANDÉ
    ├── 🧠 launch_lava_smart_kill.py         # Intelligent
    ├── 🔒 launch_lava_browser_kill.py       # Basique
    ├── 🌐 launch_lava_wsl.py                # WSL2
    └── 🎮 launch_lava.py                    # Manuel
```

## 🔧 **Configuration**

### **Paramètres LAVA configurables**
- **Généraux** : Longueur signature, nombre max primers, % match
- **Primers** : Longueurs, températures de fusion (Tm), % GC
- **Espacement** : Gaps minimums entre primers
- **Architecture** : Calcul automatique par LAVA

### **Types de primers**
- **🧬 STEM** : Primers d'ancrage avec espacement automatique
- **🔄 LOOP** : Primers de boucle avec paramètres configurables

## 📊 **Fichiers de sortie**

Chaque exécution génère :
- **`nom.primers`** : Primers principaux
- **`nom.primers.primers`** : Format étendu
- **`nom.primers.dash`** : Format tableau
- **`nom.primers.all_signatures`** : Toutes les signatures

## 🌍 **Compatibilité**

- ✅ **Linux** (Ubuntu, Debian, CentOS)
- ✅ **WSL2** (Windows 10/11)
- ✅ **macOS** (Intel & Apple Silicon via Homebrew)
- ❌ **Windows natif** (utiliser WSL2)

## 🚀 **Déploiement en Production (Serveur / Internet)**

Pour un déploiement public et sécurisé accessible via une URL internet, l'équipe informatique doit déployer l'application avec un serveur HTTP robuste (Nginx) et un gestionnaire WSGI (Gunicorn), car le serveur Flask intégré n'est pas conçu pour faire face à du trafic réel.

### **1. Pré-requis Serveur**
* Serveur Linux (Ubuntu/Debian recommandé)
* Python 3.8+ et dépendances (`pip install -r requirements.txt`)
* Nginx, Gunicorn et Perl (`sudo apt install nginx gunicorn perl`)
* Primer3 et BioPerl configurés globalement.

### **2. Variables d'Environnement de Sécurité**
Créez un fichier `.env` ou définissez ces variables pour l'utilisateur qui exécute l'application :
```bash
# Désactive la console de debug distante (Crucial contre le Remote Code Execution)
export FLASK_DEBUG=false
# Générez une clé secrète forte pour chiffrer les sessions (ex: python3 -c 'import secrets; print(secrets.token_hex(24))')
export SECRET_KEY="votre_cle_secrete_generee_aleatoirement"
```

### **3. Démarrage via Gunicorn (WSGI)**
Gunicorn remplacera la commande basique `python3 lava_flask_app.py`. Il gérera des processus concurrents.
Démarrez le service web en arrière-plan à l'intérieur du dossier du projet :
```bash
gunicorn --bind 127.0.0.1:5001 --workers 2 lava_flask_app:app
```
*(Optionnel : Configurer cela sous forme de Service `systemd` pour que LAVA redémarre tout seul si le serveur reboot).*

### **4. Configuration du Reverse Proxy Nginx (Exposition URL)**
LAVA tournant sur un port local (5001), il faut configurer le serveur web Nginx pour capter le vrai trafic internet (Port 80/443) et le router vers LAVA.
Créez un fichier `/etc/nginx/sites-available/lava_app` :
```nginx
server {
    listen 80;
    server_name lava.votre-institut.fr; # Remplacez par le nom de domaine prévu ou IP publique
    
    location / {
        proxy_pass http://127.0.0.1:5001;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # Configuration pour les longs calculs LAVA
        proxy_read_timeout 1800s;
        proxy_connect_timeout 1800s;
    }
}
```
Activez la configuration et redémarrez Nginx :
```bash
sudo ln -s /etc/nginx/sites-available/lava_app /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```
Vous pouvez ensuite installer un certificat HTTPS gratuit via `certbot`. La plateforme est désormais accessible globalement à l'URL `http://lava.votre-institut.fr` (ou l'IP du serveur).

## 🐛 **Dépannage**

### **Problèmes courants**
1. **Module Perl manquant** : `sudo cpanm Bio::Tools::Run::Primer3` (Linux) ou `cpanm Bio::Tools::Run::Primer3` (macOS)
2. **Primer3 non trouvé** : `sudo apt install primer3` (Linux) ou `brew install primer3` (macOS)
3. **Port occupé (ex: AirPlay Receiver sur macOS)** : Le port a été configuré sur 5001 par défaut pour contourner ce problème, mais vous pouvez modifier le port dans `lava_flask_app.py` et `launch_lava_*.py`
4. **Interface se ferme** : Utiliser `launch_lava_no_false_positive.py`

### **Logs et débogage**
- Vérifier les logs dans le terminal
- Consulter `DOCUMENTATION_LAVA.txt`
- Vérifier les variables d'environnement

## 📚 **Documentation**

- **📖 [Documentation complète](DOCUMENTATION_LAVA.txt)** : Guide utilisateur détaillé
- **🔧 [Guide d'installation](GUIDE_INSTALLATION_NOUVELLE_MACHINE.md)** : Installation sur nouvelle machine
- **💻 [Interface web](http://localhost:5000)** : Documentation interactive

## 🤝 **Contribution**

1. Fork le projet
2. Créer une branche feature (`git checkout -b feature/AmazingFeature`)
3. Commit les changements (`git commit -m 'Add AmazingFeature'`)
4. Push vers la branche (`git push origin feature/AmazingFeature`)
5. Ouvrir une Pull Request

## 📄 **Licence**

Ce projet utilise une structure de licence **mixte** :

- Le **moteur LAVA hérité** (`lava_loop_primer.pl`, `lava_stem_primer.pl`, modules `lib/LLNL/LAVA/` hérités) est distribué sous licence open-source **BSD 3-Clause** - Copyright (c) 2010, Lawrence Livermore National Security, LLC. (Clinton Torres).
- Les **modules Perl étendus** (`lib/LLNL/LAVA/Validator.pm`, `Core.pm`, `PipelineUtils.pm`) sont distribués sous licence open-source **BSD 3-Clause** - Copyright (c) 2026, Cheikh Talibouya.
- L'**interface graphique web et la suite de déploiement** (`lava_flask_app.py`, `templates/`, `static/`, `deployment/`) sont sous **Licence Propriétaire (Tous droits réservés)** - Copyright (c) 2026, Cheikh Talibouya. L'utilisation, modification ou distribution de l'interface graphique est strictement soumise à l'autorisation écrite préalable de l'auteur.

Voir le fichier [`LICENSE`](LICENSE) pour les termes complets.

## 🙏 **Remerciements**

- **LLNL** pour le développement de LAVA
- **BioPerl** pour les outils bioinformatiques
- **Primer3** pour la conception de primers
- **Flask** pour le framework web

## 📞 **Support**

- **Issues** : [GitHub Issues](https://github.com/Cheikht4/lava-dna-portable/issues)
- **Documentation** : `DOCUMENTATION_LAVA.txt`
- **Email** : cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn

---

**⭐ Si ce projet vous est utile, n'oubliez pas de le star sur GitHub !**
