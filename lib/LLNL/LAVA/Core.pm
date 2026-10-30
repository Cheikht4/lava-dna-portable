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
    generateSigmoidPenalty
    generateDistancePenalties
    countDegenerateBases
);

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
