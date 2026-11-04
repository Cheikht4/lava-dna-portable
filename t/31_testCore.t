#!/usr/bin/perl
# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>.
# BSD 3-Clause License - See LICENSE at root of project.
#
# Tests unitaires pour LLNL::LAVA::Core
# Unit tests for LLNL::LAVA::Core
#
# Fonctions testées / Tested functions :
#   - generateSigmoidPenalty : asymétrie (distance <= cible -> pénalité 0), croissance au-delà

use strict;
use warnings;
use lib '../lib';
use LLNL::LAVA::Core qw(generateSigmoidPenalty);

my $testCount = 0;
my $passCount = 0;

sub ok {
    my ($result, $label) = @_;
    $label //= '';
    $testCount++;
    if ($result) {
        $passCount++;
        print "ok $testCount - $label\n";
    } else {
        print "not ok $testCount - $label\n";
    }
}

BEGIN { print "1..10\n"; }

# ─────────────────────────────────────────────────────────
# SECTION 1 : generateSigmoidPenalty - Asymétrie
# Vérification que les distances <= cible donnent pénalité 0
# Verify that distances <= target yield penalty 0
# ─────────────────────────────────────────────────────────

# Distance exactement égale à la cible -> pénalité = 0
ok(generateSigmoidPenalty(24, 24) == 0, 'generateSigmoidPenalty: actual == target -> penalty 0');

# Distance inférieure à la cible (amplicon compact) -> pénalité = 0
ok(generateSigmoidPenalty(3, 24)  == 0, 'generateSigmoidPenalty: actual < target (3 vs 24) -> penalty 0');
ok(generateSigmoidPenalty(1, 18)  == 0, 'generateSigmoidPenalty: actual < target (1 vs 18) -> penalty 0');
ok(generateSigmoidPenalty(0, 10)  == 0, 'generateSigmoidPenalty: actual = 0 < target -> penalty 0');

# Distance négative -> pénalité = 100 (maximum)
ok(generateSigmoidPenalty(-1, 24) == 100, 'generateSigmoidPenalty: negative actual -> penalty 100');

# ─────────────────────────────────────────────────────────
# SECTION 2 : generateSigmoidPenalty - Croissance au-delà du seuil
# Penalty growth beyond threshold
# ─────────────────────────────────────────────────────────

# Juste au-delà du seuil (threshold + 1) -> pénalité > 0 mais < 100
# Just beyond threshold (threshold + 1) -> 0 < penalty < 100
my $just_beyond = generateSigmoidPenalty(21, 20);
ok($just_beyond > 0 && $just_beyond < 100, 'generateSigmoidPenalty: just beyond threshold (threshold+1)');

# Effet du facteur de pente k_slope : une pente plus forte donne une pénalité plus élevée
# Slope factor k_slope effect: higher slope yields higher penalty
ok(generateSigmoidPenalty(25, 20, 0.10) < generateSigmoidPenalty(25, 20, 0.30),
   'generateSigmoidPenalty: higher k_slope yields steeper penalty (25,20,0.10 < 25,20,0.30)');

# Très au-delà -> pénalité élevée (proche de 100)
my $far_beyond = generateSigmoidPenalty(200, 20);
ok($far_beyond > 90, 'generateSigmoidPenalty: far beyond target -> penalty > 90');

# ─────────────────────────────────────────────────────────
# SECTION 3 : generateSigmoidPenalty - Monotonie et asymptote de saturation
# La pénalité doit être croissante au-delà du seuil et atteindre ~98,01 au seuil de saturation
# Penalty must be monotonically increasing beyond threshold and reach ~98.01 at saturation
# ─────────────────────────────────────────────────────────

my $p30 = generateSigmoidPenalty(30, 20);
my $p40 = generateSigmoidPenalty(40, 20);
my $p60 = generateSigmoidPenalty(60, 20);
ok($p30 < $p40 && $p40 < $p60, 'generateSigmoidPenalty: monotonically increasing beyond plateau');

my $k = 4.6 / (86 - 49);
my $p_at_sat = generateSigmoidPenalty(86, 49, $k);
ok(abs($p_at_sat - 98.0096) < 0.01,
   'generateSigmoidPenalty: au seuil de saturation, penalite = 98,01 (asymptote)');

print "\n# Résultat : $passCount/$testCount tests passés\n";
print "# Result  : $passCount/$testCount tests passed\n";
