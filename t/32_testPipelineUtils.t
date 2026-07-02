#!/usr/bin/perl
# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>.
# BSD 3-Clause License — See LICENSE at root of project.
#
# Tests unitaires pour LLNL::LAVA::PipelineUtils
# Unit tests for LLNL::LAVA::PipelineUtils
#
# Fonctions testées / Tested functions :
#   - reducePrimersByWindow : logique NMS (zone d'exclusion, quota max_per_window)

use strict;
use warnings;
use lib '../lib';
use LLNL::LAVA::PipelineUtils qw(reducePrimersByWindow);

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

BEGIN { print "1..8\n"; }

# ─────────────────────────────────────────────────────────
# Mock d'un objet primer minimal pour les tests
# Minimal mock primer object for tests
# ─────────────────────────────────────────────────────────
package MockPrimer;
sub new {
    my ($class, %args) = @_;
    return bless { location => $args{location}, penalty => $args{penalty} }, $class;
}
sub getLocation { return $_[0]->{location}; }
sub getPenalty  { return $_[0]->{penalty};  }

package main;

# ─────────────────────────────────────────────────────────
# SECTION 1 : Liste vide -> retour immédiat
# Empty list -> immediate return
# ─────────────────────────────────────────────────────────
my $empty_result = reducePrimersByWindow([], 20, 2);
ok(scalar(@$empty_result) == 0, 'reducePrimersByWindow: empty list returns empty list');

# ─────────────────────────────────────────────────────────
# SECTION 2 : Un seul candidat -> toujours sélectionné
# Single candidate -> always selected
# ─────────────────────────────────────────────────────────
my @single = (MockPrimer->new(location => 50, penalty => 1.0));
my $single_result = reducePrimersByWindow(\@single, 20, 2);
ok(scalar(@$single_result) == 1, 'reducePrimersByWindow: single primer is always selected');
ok($single_result->[0]->getLocation() == 50, 'reducePrimersByWindow: single primer keeps its location');

# ─────────────────────────────────────────────────────────
# SECTION 3 : Quota max_per_window respecté
# Zone d'exclusion = window_size = 20 centré sur le premier candidat
# max_per_window = 1 -> un seul candidat par zone
# ─────────────────────────────────────────────────────────
# 3 candidats dans la même fenêtre (positions 50, 55, 58) — pénalités croissantes
# 3 candidates in the same window (positions 50, 55, 58) — increasing penalties
my @same_window = (
    MockPrimer->new(location => 50, penalty => 1.0), # meilleur / best
    MockPrimer->new(location => 55, penalty => 2.0),
    MockPrimer->new(location => 58, penalty => 3.0),
);
my $quota1_result = reducePrimersByWindow(\@same_window, 20, 1);
ok(scalar(@$quota1_result) == 1, 'reducePrimersByWindow: max_per_window=1 keeps only best in zone');
ok($quota1_result->[0]->getLocation() == 50, 'reducePrimersByWindow: selected is best candidate (lowest penalty)');

# ─────────────────────────────────────────────────────────
# SECTION 4 : max_per_window=2 -> 2 candidats acceptés dans la zone
# max_per_window=2 -> 2 candidates accepted in the zone
# ─────────────────────────────────────────────────────────
my $quota2_result = reducePrimersByWindow(\@same_window, 20, 2);
ok(scalar(@$quota2_result) == 2, 'reducePrimersByWindow: max_per_window=2 keeps 2 best in zone');

# ─────────────────────────────────────────────────────────
# SECTION 5 : Deux zones distinctes -> les meilleurs de chaque zone sélectionnés
# Two distinct zones -> best from each zone is selected
# ─────────────────────────────────────────────────────────
# Zone A : positions 10, 15 (dans fenêtre [0..20])
# Zone B : positions 100, 105 (dans fenêtre [90..110])
my @two_zones = (
    MockPrimer->new(location => 10,  penalty => 1.0), # best zone A
    MockPrimer->new(location => 15,  penalty => 2.0), # second zone A
    MockPrimer->new(location => 100, penalty => 3.0), # best zone B
    MockPrimer->new(location => 105, penalty => 4.0), # second zone B
);
my $two_zones_result = reducePrimersByWindow(\@two_zones, 20, 1);
ok(scalar(@$two_zones_result) == 2, 'reducePrimersByWindow: two distinct zones -> one candidate each');

# Les deux sélectionnés doivent être les meilleurs (penalty 1.0 et 3.0)
my @selected_penalties = sort { $a <=> $b } map { $_->getPenalty() } @$two_zones_result;
ok($selected_penalties[0] == 1.0 && $selected_penalties[1] == 3.0,
    'reducePrimersByWindow: selected are best candidates from each zone (penalty 1.0 and 3.0)');

print "\n# Résultat : $passCount/$testCount tests passés\n";
print "# Result  : $passCount/$testCount tests passed\n";
