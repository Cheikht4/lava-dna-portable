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
    my ($L) = @_;
    
    # Valeurs par défaut si L non fourni ou trop petit / Default values if L not provided or too small
    $L = 250 unless (defined $L && $L > 50);

    # Conversion en entiers pour éviter les problèmes d'arrondi plus tard / Convert to integers to avoid rounding issues later
    my $geometry = {
        'f3_f2_target' => int($L * 0.12),
        'f2_f1_target' => int($L * 0.18),
        'inner_target' => int($L * 0.40), # Distance F1c-B1c
        'b1_b2_target' => int($L * 0.18), # Symétrique
        'b2_b3_target' => int($L * 0.12)  # Symétrique
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
#   target         - Distance cible idéale proportionnelle / Ideal proportional target distance
#   plateau_ratio  - Pourcentage au-dessus de la cible toléré sans pénalité / Plateau ratio above target
#   k_slope        - Facteur de pente pour la montée de pénalité / Slope factor for penalty increase
#
sub generateSigmoidPenalty {
    my ($actual, $target, $plateau_ratio, $k_slope) = @_;
    
    # Valeurs par défaut si non fournies / Default values if not provided
    $plateau_ratio = 0.25 unless defined $plateau_ratio;
    $k_slope = 0.15 unless defined $k_slope;
    
    return 100 if $actual < 0; 
    
    # Les distances plus courtes ou égales à la cible sont idéales cinétiquement (pas de pénalité)
    # Distances shorter than or equal to the target are kinetically ideal (no penalty)
    return 0 if $actual <= $target;
    
    # Paramètres de la montée progressive au-dessus de la cible
    # Parameters for progressive rise above target
    my $L_plateau_width = $target * $plateau_ratio; # Plateau "gratuit" de + X% au-dessus de la cible
    my $max_penalty = 100;
    
    my $diff = $actual - $target;
    
    # Si l'excès reste dans le plateau de tolérance gratuit, pas de pénalité
    # If the excess remains within the free tolerance plateau, no penalty
    return 0 if $diff <= $L_plateau_width;

    # Calcul de la pénalité progressive au-delà du plateau
    # Calculation of progressive penalty beyond the plateau
    my $excess = $diff - $L_plateau_width;
    
    # Formule Sigmoïde corrigée (commence à 0 après la limite du plateau)
    # Corrected Sigmoid Formula (starts at 0 after the plateau limit)
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
    my ($maxDistance, $targetLength, $plateau_ratio, $k_slope) = @_;
    
    my @penalties = ();
    
    for (my $i = 0; $i < $maxDistance; $i++) {
        $penalties[$i] = generateSigmoidPenalty($i, $targetLength, $plateau_ratio, $k_slope);
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
