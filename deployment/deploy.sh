#!/bin/bash
# Script de déploiement automatique pour LAVA-DNA
set -e

echo "🚀 Déploiement LAVA-DNA sur serveur de production"
echo "================================================"

# Variables de configuration
DOMAIN="your-domain.com"  # Remplacez par votre domaine
EMAIL="admin@your-domain.com"  # Pour Let's Encrypt
LAVA_PATH="/opt/lava-dna"
LAVA_USER="lavauser"

# Couleurs pour l'affichage
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Vérifier que le script est exécuté en tant que root
if [[ $EUID -ne 0 ]]; then
   print_error "Ce script doit être exécuté en tant que root"
   exit 1
fi

print_status "Mise à jour du système..."
apt update && apt upgrade -y

print_status "Installation des dépendances système..."
apt install -y \
    nginx \
    python3 \
    python3-pip \
    python3-venv \
    perl \
    build-essential \
    cpanminus \
    bioperl \
    libmodule-build-perl \
    libdata-stag-perl \
    libgetopt-long-descriptive-perl \
    liblist-moreutils-perl \
    certbot \
    python3-certbot-nginx \
    redis-server \
    git \
    curl \
    htop \
    ufw

print_status "Installation des modules Perl..."
cpanm Bio::Perl Bio::SeqIO Bio::AlignIO

print_status "Création de l'utilisateur lavauser..."
if ! id "$LAVA_USER" &>/dev/null; then
    useradd -m -s /bin/bash $LAVA_USER
    print_status "Utilisateur $LAVA_USER créé"
else
    print_warning "Utilisateur $LAVA_USER existe déjà"
fi

print_status "Création des répertoires..."
mkdir -p $LAVA_PATH
mkdir -p /var/log/lava
mkdir -p /var/run/lava
chown -R $LAVA_USER:$LAVA_USER $LAVA_PATH /var/log/lava /var/run/lava

print_status "Copie des fichiers de l'application..."
# Supposons que le script est exécuté depuis le répertoire du projet
cp -r . $LAVA_PATH/
chown -R $LAVA_USER:$LAVA_USER $LAVA_PATH

print_status "Création de l'environnement virtuel Python..."
sudo -u $LAVA_USER python3 -m venv $LAVA_PATH/venv
sudo -u $LAVA_USER $LAVA_PATH/venv/bin/pip install --upgrade pip
sudo -u $LAVA_USER $LAVA_PATH/venv/bin/pip install -r $LAVA_PATH/requirements_secure.txt

print_status "Configuration Nginx..."
cp $LAVA_PATH/deployment/nginx_lava.conf /etc/nginx/sites-available/lava-dna
sed -i "s/your-domain.com/$DOMAIN/g" /etc/nginx/sites-available/lava-dna
ln -sf /etc/nginx/sites-available/lava-dna /etc/nginx/sites-enabled/
rm -f /etc/nginx/sites-enabled/default

print_status "Test de la configuration Nginx..."
nginx -t

print_status "Configuration du service systemd..."
cp $LAVA_PATH/deployment/lava-dna.service /etc/systemd/system/
systemctl daemon-reload

print_status "Configuration du firewall..."
ufw --force reset
ufw default deny incoming
ufw default allow outgoing
ufw allow ssh
ufw allow 80/tcp
ufw allow 443/tcp
ufw --force enable

print_status "Démarrage des services..."
systemctl enable redis-server
systemctl start redis-server
systemctl enable nginx
systemctl start nginx
systemctl enable lava-dna
systemctl start lava-dna

print_status "Obtention du certificat SSL avec Let's Encrypt..."
certbot --nginx -d $DOMAIN --non-interactive --agree-tos --email $EMAIL

print_status "Configuration du renouvellement automatique des certificats..."
(crontab -l 2>/dev/null; echo "0 12 * * * /usr/bin/certbot renew --quiet") | crontab -

print_status "Nettoyage automatique des fichiers temporaires..."
(crontab -u $LAVA_USER -l 2>/dev/null; echo "0 2 * * * find $LAVA_PATH/uploads -type f -mtime +7 -delete") | crontab -u $LAVA_USER -
(crontab -u $LAVA_USER -l 2>/dev/null; echo "0 3 * * * find $LAVA_PATH/results -type f -mtime +30 -delete") | crontab -u $LAVA_USER -

print_status "Configuration des logs avec logrotate..."
cat > /etc/logrotate.d/lava-dna << EOF
/var/log/lava/*.log {
    daily
    missingok
    rotate 30
    compress
    notifempty
    create 0644 $LAVA_USER $LAVA_USER
    postrotate
        systemctl reload lava-dna
    endscript
}
EOF

print_status "Génération d'un mot de passe admin sécurisé..."
ADMIN_PASSWORD=$(openssl rand -base64 32)
echo "FLASK_SECRET_KEY=$(openssl rand -hex 32)" > $LAVA_PATH/.env
echo "LAVA_ADMIN_PASSWORD=$ADMIN_PASSWORD" >> $LAVA_PATH/.env
chown $LAVA_USER:$LAVA_USER $LAVA_PATH/.env
chmod 600 $LAVA_PATH/.env

print_status "Redémarrage des services..."
systemctl restart lava-dna
systemctl restart nginx

# Vérification de l'état des services
sleep 5
if systemctl is-active --quiet lava-dna; then
    print_status "Service LAVA-DNA : ✅ Actif"
else
    print_error "Service LAVA-DNA : ❌ Erreur"
    systemctl status lava-dna
    exit 1
fi

if systemctl is-active --quiet nginx; then
    print_status "Service Nginx : ✅ Actif"
else
    print_error "Service Nginx : ❌ Erreur"
    systemctl status nginx
    exit 1
fi

echo ""
echo "🎉 Déploiement terminé avec succès !"
echo "======================================="
echo ""
echo "🌐 URL d'accès : https://$DOMAIN"
echo "👤 Compte admin : admin"
echo "🔑 Mot de passe admin : $ADMIN_PASSWORD"
echo ""
echo "📋 Commandes utiles :"
echo "  - Logs de l'application : journalctl -u lava-dna -f"
echo "  - Logs Nginx : tail -f /var/log/nginx/lava_error.log"
echo "  - Redémarrer l'app : systemctl restart lava-dna"
echo "  - Status des services : systemctl status lava-dna nginx"
echo ""
echo "⚠️  N'oubliez pas de :"
echo "  1. Sauvegarder le mot de passe admin affiché ci-dessus"
echo "  2. Configurer votre DNS pour pointer vers ce serveur"
echo "  3. Surveiller les logs lors des premiers tests"
echo ""
print_warning "IMPORTANT : Le mot de passe admin est stocké dans $LAVA_PATH/.env"
