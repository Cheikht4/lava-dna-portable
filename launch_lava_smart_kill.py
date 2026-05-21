#!/usr/bin/env python3
"""
Script de lancement LAVA avec arrêt intelligent (uniquement à la fermeture complète du navigateur)
"""
import subprocess
import sys
import os
import time
import threading
import signal
import atexit
from threading import Timer

def open_browser_wsl():
    """Ouvre le navigateur Windows depuis WSL"""
    try:
        print("🌐 Ouverture du navigateur Windows...")
        subprocess.run(['cmd.exe', '/c', 'start', 'http://localhost:5001'], 
                      stdout=subprocess.DEVNULL, 
                      stderr=subprocess.DEVNULL)
        print("✅ Navigateur ouvert dans Windows")
    except Exception as e:
        print(f"⚠️  Impossible d'ouvrir automatiquement le navigateur: {e}")
        print("🔗 Ouvrez manuellement : http://localhost:5001")

def add_smart_detection():
    """Ajoute la détection intelligente dans le template"""
    try:
        base_template_path = 'templates/base.html'
        if os.path.exists(base_template_path):
            with open(base_template_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Vérifier si la détection existe déjà
            if 'smart-browser-detection' not in content:
                detection_script = '''
    <!-- Détection intelligente fermeture navigateur -->
    <script id="smart-browser-detection">
    let heartbeatInterval;
    let isFormSubmission = false;
    let lastInteraction = Date.now();
    
    // Heartbeat moins fréquent (toutes les 60 secondes)
    function sendHeartbeat() {
        fetch('/stay-alive', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({timestamp: Date.now()})
        }).catch(() => {});
    }
    
    // Démarrer le heartbeat
    document.addEventListener('DOMContentLoaded', function() {
        heartbeatInterval = setInterval(sendHeartbeat, 60000); // 1 minute
        sendHeartbeat();
    });
    
    // Détecter les interactions utilisateur
    ['click', 'keydown', 'mousemove', 'scroll'].forEach(event => {
        document.addEventListener(event, () => {
            lastInteraction = Date.now();
        });
    });
    
    // Détecter spécifiquement les soumissions de formulaires
    document.addEventListener('submit', function() {
        isFormSubmission = true;
        // Reset après la soumission
        setTimeout(() => { isFormSubmission = false; }, 10000);
    });
    
    // Détecter les clics sur les liens et boutons
    document.addEventListener('click', function(e) {
        if (e.target.tagName === 'A' || e.target.type === 'submit' || 
            e.target.closest('button') || e.target.closest('form')) {
            isFormSubmission = true;
            setTimeout(() => { isFormSubmission = false; }, 10000);
        }
    });
    
    // Détection de fermeture plus intelligente
    window.addEventListener('beforeunload', function(e) {
        // Ne pas déclencher si c'est une soumission de formulaire récente
        if (isFormSubmission) {
            return;
        }
        
        // Ne pas déclencher si interaction récente (< 2 secondes)
        if (Date.now() - lastInteraction < 2000) {
            return;
        }
        
        // Envoyer signal de fermeture avec délai
        navigator.sendBeacon('/maybe-closing', JSON.stringify({
            timestamp: Date.now(),
            lastInteraction: lastInteraction
        }));
    });
    
    // Détecter la visibilité de la page
    document.addEventListener('visibilitychange', function() {
        if (document.hidden) {
            // Page cachée - peut être une fermeture
            setTimeout(() => {
                if (document.hidden && !isFormSubmission) {
                    // Si toujours cachée après 5 secondes et pas de formulaire
                    navigator.sendBeacon('/page-hidden', JSON.stringify({
                        timestamp: Date.now()
                    }));
                }
            }, 5000);
        }
    });
    </script>'''
                
                # Insérer avant </body>
                content = content.replace('</body>', detection_script + '\n</body>')
                
                with open(base_template_path, 'w', encoding='utf-8') as f:
                    f.write(content)
                
                print("✅ Détection intelligente ajoutée")
            else:
                print("ℹ️  Détection déjà présente")
                
    except Exception as e:
        print(f"⚠️  Erreur détection: {e}")

def add_smart_routes():
    """Ajoute les routes Flask intelligentes"""
    try:
        with open('lava_flask_app.py', 'r') as f:
            flask_content = f.read()
        
        if '/stay-alive' not in flask_content:
            routes_code = '''
# Variables globales pour le tracking
browser_sessions = {}
shutdown_timer = None

@app.route('/stay-alive', methods=['POST'])
def stay_alive():
    """Signal de vie du navigateur"""
    session_id = session.get('session_id', str(time.time()))
    session['session_id'] = session_id
    session['last_alive'] = datetime.now()
    
    browser_sessions[session_id] = datetime.now()
    
    # Annuler l'arrêt programmé s'il y en a un
    global shutdown_timer
    if shutdown_timer:
        shutdown_timer.cancel()
        shutdown_timer = None
    
    return {'status': 'alive', 'session': session_id}

@app.route('/maybe-closing', methods=['POST'])
def maybe_closing():
    """Signal possible de fermeture"""
    print("🤔 Signal possible de fermeture reçu...")
    
    def check_if_really_closed():
        time.sleep(15)  # Attendre 15 secondes
        
        # Vérifier s'il y a eu de l'activité récente
        session_id = session.get('session_id')
        if session_id and session_id in browser_sessions:
            last_alive = browser_sessions[session_id]
            time_since = (datetime.now() - last_alive).total_seconds()
            
            if time_since > 30:  # Pas d'activité depuis 30 secondes
                print("🔒 Fermeture confirmée - Arrêt de l'interface")
                os._exit(0)
            else:
                print("✅ Activité détectée - Pas de fermeture")
    
    global shutdown_timer
    if shutdown_timer:
        shutdown_timer.cancel()
    
    shutdown_timer = threading.Timer(0, check_if_really_closed)
    shutdown_timer.start()
    
    return {'status': 'checking'}

@app.route('/page-hidden', methods=['POST'])
def page_hidden():
    """Page cachée depuis longtemps"""
    print("👁️  Page cachée détectée...")
    
    def delayed_check():
        time.sleep(60)  # Attendre 1 minute
        
        session_id = session.get('session_id')
        if session_id and session_id in browser_sessions:
            last_alive = browser_sessions[session_id]
            time_since = (datetime.now() - last_alive).total_seconds()
            
            if time_since > 120:  # Pas d'activité depuis 2 minutes
                print("💤 Inactivité prolongée - Arrêt de l'interface")
                os._exit(0)
    
    threading.Thread(target=delayed_check, daemon=True).start()
    return {'status': 'monitoring'}

'''
            
            # Ajouter les imports nécessaires
            if 'import threading' not in flask_content:
                flask_content = 'import threading\n' + flask_content
            if 'import time' not in flask_content:
                flask_content = 'import time\n' + flask_content
            
            # Insérer les routes
            lines = flask_content.split('\n')
            new_lines = []
            
            for line in lines:
                if line.strip().startswith('if __name__ == \'__main__\''):
                    new_lines.append(routes_code)
                new_lines.append(line)
            
            with open('lava_flask_app.py', 'w') as f:
                f.write('\n'.join(new_lines))
            
            print("✅ Routes intelligentes ajoutées")
        else:
            print("ℹ️  Routes déjà présentes")
            
    except Exception as e:
        print(f"⚠️  Erreur routes: {e}")

def cleanup_and_exit():
    """Nettoyage et arrêt propre"""
    print("\n🧹 Nettoyage...")
    
    try:
        subprocess.run(['pkill', '-f', 'lava_flask_app.py'], 
                     stderr=subprocess.DEVNULL, 
                     stdout=subprocess.DEVNULL)
    except:
        pass
    
    print("✅ Interface fermée")
    os._exit(0)

def main():
    print("🧬 LAVA-DNA Interface - Détection Intelligente")
    print("=============================================")
    print("🔒 Arrêt automatique SEULEMENT si:")
    print("   - Vraie fermeture du navigateur")
    print("   - Inactivité > 2 minutes")
    print("   - Ctrl+C dans le terminal")
    print("=============================================")
    
    if not os.path.exists('lava_env'):
        print("❌ Environnement virtuel non trouvé")
        return 1
    
    try:
        print("✅ Configuration de la détection intelligente...")
        
        add_smart_detection()
        add_smart_routes()
        
        print("🚀 Lancement de Flask...")
        
        env = os.environ.copy()
        env['PATH'] = os.path.abspath('lava_env/bin') + ':' + env['PATH']
        env['VIRTUAL_ENV'] = os.path.abspath('lava_env')
        
        flask_process = subprocess.Popen([
            'lava_env/bin/python', 'lava_flask_app.py'
        ], env=env)
        
        time.sleep(3)
        
        # Ouvrir le navigateur
        timer = Timer(2.0, open_browser_wsl)
        timer.daemon = True
        timer.start()
        
        print("📱 Interface : http://localhost:5001")
        print("🧠 Détection intelligente activée")
        print("✅ Les actualisations/soumissions ne ferment PAS l'interface")
        print("⚠️  Ctrl+C pour arrêt manuel")
        print("=" * 50)
        
        flask_process.wait()
        return 0
        
    except KeyboardInterrupt:
        print("\n🛑 Arrêt manuel")
        return 0
    except Exception as e:
        print(f"❌ Erreur : {e}")
        return 1
    finally:
        cleanup_and_exit()

# Gestion des signaux
signal.signal(signal.SIGINT, lambda s, f: cleanup_and_exit())
signal.signal(signal.SIGTERM, lambda s, f: cleanup_and_exit())
atexit.register(cleanup_and_exit)

if __name__ == '__main__':
    sys.exit(main())
