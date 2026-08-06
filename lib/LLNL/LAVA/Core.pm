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

package LLNL::LAVA::Core;

use strict;
use warnings;
use vars qw(@ISA @EXPORT_OK);

require Exporter;
@ISA = qw(Exporter);

@EXPORT_OK = qw(
    calculate_proportional_geometry
    generateSigmoidPenalty
    generateDistancePenalties
    countDegenerateBases
);

#=============================================================================
# CALCULATE PROPORTIONAL GEOMETRY
#=============================================================================
# Calcule les distances cibles basées sur la longueur totale de la signature. / Calculates target distances based on total signature length.
# Ratios: F3-F2 (12%), F2-F1 (18%), F1-B1 (40%)
# L'utilisateur fournit L (longueur totale estimée de la signature LAMP) / User provides L (estimated total length of LAMP signature)
sub calculate_proportional_geometry {
    my ($L, $max_L, $sum_primer_lengths) = @_;
    
    # Valeurs par défaut si non fournies
    $L = 250 unless (defined $L && $L > 50);
    $max_L = $L unless (defined $max_L && $max_L >= $L);
    $sum_primer_lengths = 144 unless defined $sum_primer_lengths; # 8 amorces de 18nt par défaut en LOOP
    
    # Calcul de l'espace libre réel / Calculate real available spacing
    my $available_total = $L - $sum_primer_lengths;
    $available_total = 10 if $available_total < 10; # Sécurité / Safety
    
    my $available_max = $max_L - $sum_primer_lengths;
    $available_max = $available_total if $available_max < $available_total;

    # Conversion en entiers pour éviter les problèmes d'arrondi
    my $geometry = {
        'f3_f2_target'       => int($available_total * 0.12),
        'f3_f2_borne_haute'  => int($available_max * 0.12),
        'f2_f1_target'       => int($available_total * 0.18),
        'f2_f1_borne_haute'  => int($available_max * 0.18),
        'inner_target'       => int($available_total * 0.40), # Distance F1c-B1c
        'inner_borne_haute'  => int($available_max * 0.40),
        'b1_b2_target'       => int($available_total * 0.18),
        'b1_b2_borne_haute'  => int($available_max * 0.18),
        'b2_b3_target'       => int($available_total * 0.12),
        'b2_b3_borne_haute'  => int($available_max * 0.12)
    };

    return $geometry;
}

#=============================================================================
# GENERATE SIGMOID PENALTY (ASYMMETRIC)
#=============================================================================
# Calcule une pénalité basée sur une courbe sigmoïde asymétrique.
# Calculates a penalty based on an asymmetric sigmoid curve.
#
# Logique biologique / Biological logic:
# - Les distances plus courtes que la cible sont favorables stériquement et cinétiquement -> Pénalité = 0.
#   Distances shorter than the target are sterically and kinetically favorable -> Penalty = 0.
# - Les distances plus longues (amplicon étiré) nuisent à la réaction -> Pénalité progressive.
#   Longer distances (stretched amplicon) hinder the reaction -> Progressive penalty.
#
# Paramètres / Parameters:
#   actual         - Distance réelle observée / Actual observed distance
#   threshold      - Distance limite gratuite / Free tolerance threshold
#   k_slope        - Facteur de pente pour la montée de pénalité / Slope factor for penalty increase
#
sub generateSigmoidPenalty {
    my ($actual, $threshold, $k_slope) = @_;
    
    # Valeurs par défaut si non fournies / Default values if not provided
    $k_slope = 0.15 unless defined $k_slope;
    
    return 100 if $actual < 0; 
    
    # Les distances plus courtes ou égales à la limite sont idéales cinétiquement (pas de pénalité)
    return 0 if $actual <= $threshold;
    
    my $max_penalty = 100;
    my $excess = $actual - $threshold;
    
    # Formule Sigmoïde corrigée (commence à 0 après la limite)
    # P(x) = max_penalty * [ (2 / (1 + exp(-k * x))) - 1 ]
    my $penalty = $max_penalty * ( (2 / (1 + exp(-$k_slope * $excess))) - 1 );
    
    return $penalty;
}

#=============================================================================
# GENERATE DISTANCE PENALTIES (MODERNIZED)
#=============================================================================
# Remplace l'ancienne fonction basée sur les paraboles. / Replaces the old parabola-based function.
# Génère un tableau de pénalités pour toutes les distances possibles jusqu'à maxDistance. / Generates an array of penalties for all possible distances up to maxDistance.
sub generateDistancePenalties {
    my ($maxDistance, $threshold, $k_slope) = @_;
    
    my @penalties = ();
    
    for (my $i = 0; $i < $maxDistance; $i++) {
        $penalties[$i] = generateSigmoidPenalty($i, $threshold, $k_slope);
    }
    
    return \@penalties;
}

#=============================================================================
# COUNT DEGENERATE BASES
#=============================================================================
# Compte le nombre de bases non-standard (non A, C, G, T) dans une chaine.
# Utilise pour trier les signatures par "proprete".
sub countDegenerateBases {
    my ($sequence) = @_;
    return 0 unless defined $sequence;
    
    # Compter tout ce qui n'est pas A, C, G, T (insensible a la casse)
    my $count = ($sequence =~ tr/BDHVKMNRSWYbdhvkmnrswy//);
    return $count;
}

1;
