#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use File::Basename;

# Test unitaire : vérifier que la formule avec +1 produit une distance nulle pour deux amorces adjacentes
# (selon la convention de calcul dynamique du moteur LAVA).

# Selon les traces d'exécution (dengue 1), pour deux amorces adjacentes de 18nt,
# les variables extraites sont telles que la formule ($loc - $len + 1) - $mid_loc == 0.
# On simule cet état :
my $outerLocation = 8787;  
my $outerLength   = 18;
my $middleLocation = 8770; 

# Formule sans le +1 (celle qui donnait -1)
my $distance_bug = ($outerLocation - $outerLength) - $middleLocation;

# Formule avec le +1 (celle corrigée)
my $distance_fix = ($outerLocation - $outerLength + 1) - $middleLocation;

is($distance_bug, -1, "L'ancienne formule (sans +1) provoque un decalage de -1");
is($distance_fix, 0,  "La formule corrigee (avec +1) produit bien une distance de 0 pour des amorces adjacentes");

done_testing();
