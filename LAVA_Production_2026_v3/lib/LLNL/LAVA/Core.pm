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
# Calcule les distances cibles basées sur la longueur totale de la signature.
# Ratios: F3-F2 (12%), F2-F1 (18%), F1-B1 (40%)
# L'utilisateur fournit L (longueur totale estimée de la signature LAMP)
sub calculate_proportional_geometry {
    my ($L) = @_;
    
    # Valeurs par défaut si L non fourni ou trop petit
    $L = 250 unless (defined $L && $L > 50);

    # Conversion en entiers pour éviter les problèmes d'arrondi plus tard
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
# GENERATE SIGMOID PENALTY
#=============================================================================
# Calcule une pénalité basée sur une courbe sigmoïde.
# Zone de confort (plateau 0) : +/- 15% de la cible.
# Au-delà : augmentation fluide logistique.
# SIGMOÏDE GÉNÉRALISÉE PERMISSIVE (LAVA 2026 UPDATE)
sub generateSigmoidPenalty {
    my ($actual, $target, $plateau_ratio, $k_slope) = @_;
    
    # Valeurs par défaut si non fournies
    $plateau_ratio = 0.25 unless defined $plateau_ratio;
    $k_slope = 0.15 unless defined $k_slope;
    
    return 100 if $actual < 0; 
    
    # Paramètres du modèle "Colline Douce"
    my $L_plateau_width = $target * $plateau_ratio; # Plateau "gratuit" de +/- X%
    my $max_penalty = 100;
    
    my $diff = abs($actual - $target);
    
    # Formule Sigmoïde Généralisée : P(x) = 100 / (1 + exp(-k * (|x - T| - L)))
    # Si diff < L, l'exposant est positif et grand -> exp() grand -> P ~ 0
    # Si diff > L, l'exposant devient négatif -> exp() petit -> P augmente vers 100
    
    # Note : Pour assurer un VRAI zéro dans le plateau, on garde une condition explicite
    return 0 if $diff <= $L_plateau_width;

    # Calcul de la pénalité progressive au-delà du plateau
    # On décale 'x' de la largeur du plateau pour que la montée commence à 0 après la limite
    my $excess = $diff - $L_plateau_width;
    
    # Nouvelle Formule Sigmoïde corrigée (commence à 0 mathématiquement)
    # P(x) = max_penalty * [ (2 / (1 + exp(-k * x))) - 1 ]
    # Si x = 0 => P(0) = max_penalty * [ 2/2 - 1 ] = 0
    # Si x = inf => P(inf) = max_penalty * [ 2/1 - 1 ] = max_penalty
    my $penalty = $max_penalty * ( (2 / (1 + exp(-$k_slope * $excess))) - 1 );
    
    return $penalty;
}

#=============================================================================
# GENERATE DISTANCE PENALTIES (MODERNIZED)
#=============================================================================
# Remplace l'ancienne fonction basée sur les paraboles.
# Génère un tableau de pénalités pour toutes les distances possibles jusqu'à maxDistance.
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
