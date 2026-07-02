#!/usr/bin/perl
# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>.
# BSD 3-Clause License — See LICENSE at root of project.
#
# Tests unitaires pour LLNL::LAVA::Core
# Unit tests for LLNL::LAVA::Core
#
# Fonctions testées / Tested functions :
#   - generateSigmoidPenalty : asymétrie (distance <= cible -> pénalité 0), croissance au-delà
#   - calculate_proportional_geometry : ratios 12%/18%/40%

use strict;
use warnings;
use lib '../lib';
use LLNL::LAVA::Core qw(generateSigmoidPenalty calculate_proportional_geometry);

my $testCount = 0;
my $passCount = 0;

sub ok {
    my ($result, $label) = @_;
    $testCount++;
    if ($result) {
        $passCount++;
        print "ok $testCount - $label\n";
    } else {
        print "not ok $testCount - $label\n";
    }
}

BEGIN { print "1..14\n"; }

# ─────────────────────────────────────────────────────────
# SECTION 1 : generateSigmoidPenalty — Asymétrie
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
# SECTION 2 : generateSigmoidPenalty — Croissance au-delà du plateau
# Penalty growth beyond plateau
# ─────────────────────────────────────────────────────────

# Distance dans le plateau de tolérance (target + 25% max) -> pénalité = 0
# Cible = 20, plateau_ratio = 0.25 => plateau jusqu'à 20 + 5 = 25
ok(generateSigmoidPenalty(25, 20) == 0, 'generateSigmoidPenalty: within plateau (target+25%) -> penalty 0');

# Juste au-delà du plateau -> pénalité > 0 mais < 100
my $just_beyond = generateSigmoidPenalty(26, 20);
ok($just_beyond > 0 && $just_beyond < 100, 'generateSigmoidPenalty: just beyond plateau -> 0 < penalty < 100');

# Très au-delà -> pénalité élevée (proche de 100)
my $far_beyond = generateSigmoidPenalty(200, 20);
ok($far_beyond > 90, 'generateSigmoidPenalty: far beyond target -> penalty > 90');

# ─────────────────────────────────────────────────────────
# SECTION 3 : generateSigmoidPenalty — Monotonie
# La pénalité doit être croissante au-delà du plateau
# Penalty must be monotonically increasing beyond plateau
# ─────────────────────────────────────────────────────────

my $p30 = generateSigmoidPenalty(30, 20);
my $p40 = generateSigmoidPenalty(40, 20);
my $p60 = generateSigmoidPenalty(60, 20);
ok($p30 < $p40 && $p40 < $p60, 'generateSigmoidPenalty: monotonically increasing beyond plateau');

# ─────────────────────────────────────────────────────────
# SECTION 4 : calculate_proportional_geometry
# Vérification des ratios LAMP / LAMP ratio verification
# F3-F2: 12%, F2-F1: 18%, F1c-B1c: 40%, symétrique côté B
# ─────────────────────────────────────────────────────────

my $g200 = calculate_proportional_geometry(200);

# F3-F2 = 12% de 200 = 24
ok($g200->{'f3_f2_target'} == 24, 'calc_proportional_geometry L=200: f3_f2_target = 24 (12%)');

# F2-F1 = 18% de 200 = 36
ok($g200->{'f2_f1_target'} == 36, 'calc_proportional_geometry L=200: f2_f1_target = 36 (18%)');

# F1c-B1c = 40% de 200 = 80
ok($g200->{'inner_target'} == 80, 'calc_proportional_geometry L=200: inner_target = 80 (40%)');

# Valeur par défaut si L trop petit / Default if L too small
my $g_default = calculate_proportional_geometry(10);  # < 50 -> défaut L=250
ok($g_default->{'f3_f2_target'} == int(250 * 0.12), 'calc_proportional_geometry: L<50 uses default L=250');

print "\n# Résultat : $passCount/$testCount tests passés\n";
print "# Result  : $passCount/$testCount tests passed\n";
