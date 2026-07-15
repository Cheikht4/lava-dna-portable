# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of the copyright holder nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
# Part of the LAVA-DNA project. See LICENSE for full terms.
# Fait partie du projet LAVA-DNA. Voir LICENSE pour les termes complets.

package LLNL::LAVA::ForkManager;

use strict;
use warnings;
use POSIX qw(:sys_wait_h WNOHANG);
use File::Path qw(mkpath rmtree);
use File::Spec;
use Storable qw(nstore retrieve);

# Constructeur du gestionnaire de multi-processus native par fork
sub new {
    my ($class, $max_processes, $temp_dir) = @_;
    
    my $auto_count = get_auto_cpu_count();
    $auto_count = 1 if $auto_count < 1;
    
    # Validation et gestion du mode automatique ou cas non numerique ('abc')
    if (!defined $max_processes || $max_processes =~ /^auto$/i || $max_processes !~ /^-?\d+$/) {
        $max_processes = $auto_count;
    } else {
        $max_processes = int($max_processes);
    }
    
    if ($max_processes <= 0) {
        $max_processes = $auto_count;
    }
    
    # Plafond dur (defense en profondeur) : ne jamais forker plus que le nombre de coeurs disponibles (auto_count)
    if ($max_processes > $auto_count) {
        warn "[ForkManager] Plafond de processus depasse ($max_processes > $auto_count). Ajustement au maximum raisonnable : $auto_count.\n";
        $max_processes = $auto_count;
    }
    $max_processes = 1 if $max_processes < 1;
    
    if (!defined $temp_dir || $temp_dir eq '') {
        $temp_dir = File::Spec->catdir("results", "tmp_threads_" . $$ . "_" . time());
    }
    
    # Creer le repertoire temporaire si necessaire pour l echange de donnees
    if ($max_processes > 1) {
        mkpath($temp_dir) unless -d $temp_dir;
    }
    
    my $self = {
        max_processes  => $max_processes,
        temp_dir       => $temp_dir,
        children       => {},
        on_finish      => undef,
        current_id     => undef,
        in_child       => 0,
        parent_pid     => $$,
    };
    
    bless $self, $class;
    return $self;
}

# Fonction utilitaire pour detecter automatiquement le nombre de coeurs CPU
sub get_auto_cpu_count {
    my $cpus = 1;
    if ($^O eq 'darwin') {
        my $out = `sysctl -n hw.ncpu 2>/dev/null`;
        chomp($out) if defined $out;
        $cpus = $out if (defined $out && $out =~ /^\d+$/ && $out > 0);
    } else {
        my $out = `nproc 2>/dev/null`;
        chomp($out) if defined $out;
        if (!defined $out || $out !~ /^\d+$/ || $out <= 0) {
            $out = `grep -c ^processor /proc/cpuinfo 2>/dev/null`;
            chomp($out) if defined $out;
        }
        $cpus = $out if (defined $out && $out =~ /^\d+$/ && $out > 0);
    }
    # En mode auto, on laisse 1 coeur libre par securite pour le systeme et le serveur web
    $cpus = $cpus - 1;
    $cpus = 1 if $cpus < 1;
    return $cpus;
}

# Enregistre un callback execute par le parent lors de la reception des resultats d un enfant
sub run_on_finish {
    my ($self, $code) = @_;
    $self->{on_finish} = $code;
}

# Demarre un nouveau processus ou execute localement si max_processes == 1
sub start {
    my ($self, $id) = @_;
    
    $id = time() unless defined $id;
    $self->{current_id} = $id;
    
    # Mode sequentiel mono-processus si 1 seul thread requis
    if ($self->{max_processes} <= 1) {
        $self->{in_child} = 0;
        return 0; # 0 indique a la boucle d executer le code localement
    }
    
    # Si le nombre maximal de processus enfants est atteint, attendre qu un enfant se termine
    while (scalar(keys %{$self->{children}}) >= $self->{max_processes}) {
        $self->wait_one_child();
    }
    
    # Creation du processus enfant
    my $pid = fork();
    if (!defined $pid) {
        warn "[ForkManager] Erreur de creation du processus enfant par fork ($!). Execution locale.\n";
        $self->{in_child} = 0;
        return 0;
    }
    
    if ($pid == 0) {
        # Nous sommes dans le processus enfant
        $self->{in_child} = 1;
        return 0; # 0 indique a la boucle du processus enfant d executer sa tranche
    } else {
        # Nous sommes dans le processus parent
        $self->{children}->{$pid} = $id;
        return $pid; # Le parent recoit le PID enfant (> 0) et passe au next de la boucle
    }
}

# Termine le processus enfant et sauvegarde les donnees produites
sub finish {
    my ($self, $exit_code, $data_ref) = @_;
    $exit_code = 0 unless defined $exit_code;
    
    # Mode sequentiel (pas d enfant) : on execute directement le callback du parent
    if ($self->{max_processes} <= 1 || !$self->{in_child}) {
        if (defined $self->{on_finish}) {
            $self->{on_finish}->($$, $exit_code, $self->{current_id}, 0, 0, $data_ref);
        }
        return;
    }
    
    # Mode multi-processus : le processus enfant serialise ses resultats et quitte
    if (defined $data_ref && defined $self->{temp_dir} && -d $self->{temp_dir}) {
        my $tmp_file = File::Spec->catfile($self->{temp_dir}, "pfm_" . $$ . "_" . $self->{current_id} . ".dat");
        eval {
            nstore($data_ref, $tmp_file);
        };
        if ($@) {
            warn "[ForkManager] Erreur de serialisation dans l enfant $$ : $@\n";
        }
    }
    
    CORE::exit($exit_code);
}

# Attend et traite la fin d au moins un processus enfant
sub wait_one_child {
    my ($self) = @_;
    return 0 if scalar(keys %{$self->{children}}) == 0;
    
    my $pid = waitpid(-1, 0);
    return 0 if $pid <= 0;
    
    if (exists $self->{children}->{$pid}) {
        my $id = delete $self->{children}->{$pid};
        my $exit_code = $? >> 8;
        my $exit_signal = $? & 127;
        my $core_dump = $? & 128 ? 1 : 0;
        
        my $data_ref = undef;
        my $tmp_file = File::Spec->catfile($self->{temp_dir}, "pfm_" . $pid . "_" . $id . ".dat");
        if (-f $tmp_file) {
            eval {
                $data_ref = retrieve($tmp_file);
            };
            if ($@) {
                warn "[ForkManager] Erreur de lecture des donnees de l enfant $pid : $@\n";
            }
            unlink($tmp_file);
        }
        
        if (defined $self->{on_finish}) {
            $self->{on_finish}->($pid, $exit_code, $id, $exit_signal, $core_dump, $data_ref);
        }
        return $pid;
    }
    return 0;
}

# Attend la fin de tous les processus enfants et nettoie le repertoire temporaire
sub wait_all_children {
    my ($self) = @_;
    while (scalar(keys %{$self->{children}}) > 0) {
        $self->wait_one_child();
    }
    
    # Nettoyer le repertoire temporaire
    if (defined $self->{temp_dir} && -d $self->{temp_dir} && $self->{parent_pid} == $$) {
        rmtree($self->{temp_dir});
    }
}

# Destructeur : securite pour supprimer le dossier temporaire
sub DESTROY {
    my ($self) = @_;
    if (defined $self->{temp_dir} && -d $self->{temp_dir} && $self->{parent_pid} == $$) {
        rmtree($self->{temp_dir});
    }
}

1;
