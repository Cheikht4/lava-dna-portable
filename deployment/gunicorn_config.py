# Configuration Gunicorn pour LAVA-DNA
import multiprocessing
import os

# Serveur
bind = "127.0.0.1:8000"
# CRITIQUE: workers DOIT rester à 1. La file d'attente (FIFO) est gérée en mémoire 
# dans le processus Flask. Augmenter les workers créerait des files distinctes,
# cassant l'ordonnancement et démultipliant silencieusement les plafonds de cœurs.
workers = 1
threads = 4
worker_class = "gthread"
worker_connections = 1000
max_requests = 1000
max_requests_jitter = 100
timeout = 120  # 120 secondes pour les lancements LAVA
keepalive = 2

# Sécurité
user = "lavauser"
group = "lavauser"
umask = 0o077

# Logging
accesslog = "/var/log/lava/gunicorn_access.log"
errorlog = "/var/log/lava/gunicorn_error.log"
loglevel = "info"
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(D)s'

# Process naming
proc_name = "lava-dna-app"

# Paths
pythonpath = "/opt/lava-dna"
chdir = "/opt/lava-dna"

# PID file
pidfile = "/var/run/lava/gunicorn.pid"

# Daemon
daemon = False  # Géré par systemd

# Preload pour les performances
preload_app = True

# Hooks
def on_starting(server):
    server.log.info("Starting LAVA-DNA application")

def on_reload(server):
    server.log.info("Reloading LAVA-DNA application")

def worker_int(worker):
    worker.log.info("worker received INT or QUIT signal")

def pre_fork(server, worker):
    server.log.info("Worker spawned (pid: %s)", worker.pid)

def post_fork(server, worker):
    server.log.info("Worker spawned (pid: %s)", worker.pid)

def worker_abort(worker):
    worker.log.info("worker received SIGABRT signal")
