# 🧬 LAVA-DNA Interface Graphique

Interface web simple et portable pour les scripts LAVA STEM et LOOP.

## 🚀 Démarrage Rapide

### Option 1 : Lancement automatique (Recommandé)
```bash
python launch_gui.py
```

### Option 2 : Lancement manuel
```bash
# 1. Installer les dépendances
pip install -r requirements.txt

# 2. Lancer l'interface
streamlit run lava_gui.py
```

L'interface s'ouvrira automatiquement dans votre navigateur à l'adresse : **http://localhost:8501**

## 📋 Utilisation

### 1. Fichier FASTA
- Uploadez votre fichier FASTA contenant les séquences
- L'aperçu s'affiche automatiquement

### 2. Configuration
- **Type de script** : Choisissez STEM ou LOOP dans la barre latérale
- **Nombre max de primers** : Limite les primers générés (défaut: 5000)
- **% Match minimum** : Correspondance stricte requise (défaut: 70%)
- **% IUPAC minimum** : Couverture IUPAC requise (défaut: 80%)

### 3. Exécution
- Cliquez sur "🚀 Lancer LAVA"
- Suivez le monitoring en temps réel
- Les résultats apparaissent automatiquement

### 4. Résultats
- Téléchargez les fichiers générés (.primers, .all_signatures, .dash)
- Consultez l'aperçu directement dans l'interface

## ⚙️ Fonctionnalités

- ✅ **Interface simple** : Design épuré et intuitif
- ✅ **Support STEM & LOOP** : Commutation automatique entre scripts
- ✅ **Paramètres essentiels** : Configuration des options principales
- ✅ **Monitoring temps réel** : Suivi de l'exécution avec statut
- ✅ **Gestion résultats** : Téléchargement et aperçu intégrés
- ✅ **Cross-platform** : Fonctionne sur Windows, Mac, Linux
- ✅ **Portable** : Un seul dossier à copier

## 🛠️ Prérequis

- Python 3.7+
- Scripts LAVA (lava_stem_primer.pl, lava_loop_primer.pl)
- Bibliothèques LAVA (dossier lib/)
- Primer3 installé et accessible

## 📁 Structure des Fichiers

```
lava-dna-master/
├── 🎨 Interface
│   ├── lava_gui.py              # Interface principale
│   ├── launch_gui.py            # Script de lancement
│   ├── requirements.txt         # Dépendances Python
│   └── README_Interface.md      # Ce guide
├── 🧬 Scripts LAVA
│   ├── lava_stem_primer.pl      # Script STEM
│   ├── lava_loop_primer.pl      # Script LOOP
│   └── lib/                     # Bibliothèques Perl
└── 📊 Résultats (générés)
    ├── results_stem_*.primers
    └── results_loop_*.primers
```

## 🐛 Dépannage

### Port 8501 occupé
```bash
# Arrêter les processus Streamlit existants
pkill -f streamlit

# Ou utiliser un autre port
streamlit run lava_gui.py --server.port 8502
```

### Erreur "Scripts LAVA non trouvés"
- Vérifiez que vous êtes dans le bon dossier (lava-dna-master)
- Assurez-vous que les fichiers .pl et lib/ sont présents

### Erreur Perl/LAVA
- Vérifiez que PERL5LIB pointe vers ./lib
- Assurez-vous que Primer3 est installé : `/usr/bin/primer3_core`

## 💡 Avantages de cette Interface

1. **Simplicité** : Interface épurée focalisée sur l'essentiel
2. **Rapidité** : Paramètres pré-configurés pour un usage immédiat  
3. **Portabilité** : Zéro compilation, fonctionne partout
4. **Maintenabilité** : Code Python simple, facilement modifiable
5. **Extension** : Base solide pour ajouter des fonctionnalités

---

**🧬 Interface créée pour optimiser l'utilisation de LAVA-DNA**
