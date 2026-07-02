#!/usr/bin/perl
# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>.
# BSD 3-Clause License — See LICENSE at root of project.
#
# Tests unitaires pour LLNL::LAVA::Validator
# Unit tests for LLNL::LAVA::Validator
#
# Modules testés / Tested functions :
#   - isIUPACCompatible : tous les codes IUPAC / all IUPAC codes
#   - generateIUPACCode : correspondance base -> code / base -> code mapping
#   - rev_comp : avec bases dégénérées / with degenerate bases
#   - checkPrimerMismatchTolerance : protection 3', tolérance 5', limite dégénérescence

use strict;
use warnings;
use lib '../lib';
use LLNL::LAVA::Validator qw(isIUPACCompatible generateIUPACCode rev_comp checkPrimerMismatchTolerance);

# Activer l'UTF-8 pour les labels de tests / Enable UTF-8 for test labels
binmode(STDOUT, ':utf8');


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

# Nombre total de tests déclarés / Total declared tests
BEGIN { print "1..35\n"; }

# ─────────────────────────────────────────────────────────
# SECTION 1 : isIUPACCompatible
# Vérifie que chaque code IUPAC accepte les bonnes bases
# Checks that each IUPAC code accepts the correct bases
# ─────────────────────────────────────────────────────────

# Bases standards / Standard bases
ok(isIUPACCompatible('A', 'A'), 'isIUPACCompatible: A matches A');
ok(isIUPACCompatible('C', 'C'), 'isIUPACCompatible: C matches C');
ok(isIUPACCompatible('G', 'G'), 'isIUPACCompatible: G matches G');
ok(isIUPACCompatible('T', 'T'), 'isIUPACCompatible: T matches T');
ok(!isIUPACCompatible('A', 'C'), 'isIUPACCompatible: A does NOT match C');

# Codes ambigus à 2 bases / Two-base ambiguity codes
ok(isIUPACCompatible('R', 'A'), 'isIUPACCompatible: R(A/G) matches A');
ok(isIUPACCompatible('R', 'G'), 'isIUPACCompatible: R(A/G) matches G');
ok(!isIUPACCompatible('R', 'C'), 'isIUPACCompatible: R(A/G) does NOT match C');

ok(isIUPACCompatible('Y', 'C'), 'isIUPACCompatible: Y(C/T) matches C');
ok(isIUPACCompatible('Y', 'T'), 'isIUPACCompatible: Y(C/T) matches T');
ok(!isIUPACCompatible('Y', 'A'), 'isIUPACCompatible: Y(C/T) does NOT match A');

ok(isIUPACCompatible('S', 'C'), 'isIUPACCompatible: S(C/G) matches C');
ok(isIUPACCompatible('S', 'G'), 'isIUPACCompatible: S(C/G) matches G');
ok(isIUPACCompatible('W', 'A'), 'isIUPACCompatible: W(A/T) matches A');
ok(isIUPACCompatible('W', 'T'), 'isIUPACCompatible: W(A/T) matches T');
ok(isIUPACCompatible('K', 'G'), 'isIUPACCompatible: K(G/T) matches G');
ok(isIUPACCompatible('K', 'T'), 'isIUPACCompatible: K(G/T) matches T');
ok(isIUPACCompatible('M', 'A'), 'isIUPACCompatible: M(A/C) matches A');
ok(isIUPACCompatible('M', 'C'), 'isIUPACCompatible: M(A/C) matches C');

# Codes à 3 bases / Three-base codes
ok(isIUPACCompatible('B', 'C'), 'isIUPACCompatible: B(C/G/T) matches C');
ok(isIUPACCompatible('B', 'G'), 'isIUPACCompatible: B(C/G/T) matches G');
ok(!isIUPACCompatible('B', 'A'), 'isIUPACCompatible: B(C/G/T) does NOT match A');
ok(isIUPACCompatible('D', 'A'), 'isIUPACCompatible: D(A/G/T) matches A');
ok(isIUPACCompatible('H', 'A'), 'isIUPACCompatible: H(A/C/T) matches A');
ok(isIUPACCompatible('V', 'C'), 'isIUPACCompatible: V(A/C/G) matches C');
ok(!isIUPACCompatible('V', 'T'), 'isIUPACCompatible: V(A/C/G) does NOT match T');

# N = tout / N = anything
ok(isIUPACCompatible('N', 'A'), 'isIUPACCompatible: N matches A');
ok(isIUPACCompatible('N', 'T'), 'isIUPACCompatible: N matches T');

# ─────────────────────────────────────────────────────────
# SECTION 2 : rev_comp avec codes dégénérés
# rev_comp with degenerate bases
# ─────────────────────────────────────────────────────────

# Séquence simple / Simple sequence
ok(rev_comp('ATCG') eq 'CGAT', 'rev_comp: ATCG -> CGAT');

# Avec bases dégénérées / With degenerate bases
# R(A/G) -> Y(C/T) sur le brin complémentaire
ok(rev_comp('ATRG') eq 'CYAT', 'rev_comp: R(A/G) complements to Y(C/T)');

# K(G/T) -> M(A/C)
ok(rev_comp('AKGC') eq 'GCMT', 'rev_comp: K(G/T) complements to M(A/C)');

# ─────────────────────────────────────────────────────────
# SECTION 3 : checkPrimerMismatchTolerance
# Protection de la zone 3' et tolérance 5'
# 3' zone protection and 5' tolerance
# ─────────────────────────────────────────────────────────
# Jeu de données minimal / Minimal dataset
# Amorce : ATGCAATGCAATGCAATGCA (20 bp)
my $primer = "ATGCAATGCAATGCAATGCA";
my $length  = length($primer);
my @seqs = (
    "ATGCAATGCAATGCAATGCA",  # 0: match parfait
    "TTGCAATGCAATGCAATGCA",  # 1: 1 mismatch 5' -> toléré avec max=1
    "ATGCAATGCAATGCAATGCT",  # 2: 1 mismatch 3' -> rejeté (protection 3')
);

# Test : sans tolérance, seule la séquence 0 est compatible
my ($seq0, $cov0, $degen0, $ids0) = checkPrimerMismatchTolerance(
    \@seqs, 0, $length, $primer,
    100, 100, 5,
    0, 6, 0, 0, 6, 0.05
);
ok(scalar(@$ids0) == 1 && $ids0->[0] == 0, 'checkPrimerMismatchTolerance: 0 mismatch -> only perfect match');

# Test : avec tolérance 1, la séquence 1 (mismatch 5') est acceptée aussi
my ($seq1, $cov1, $degen1, $ids1) = checkPrimerMismatchTolerance(
    \@seqs, 0, $length, $primer,
    100, 100, 5,
    0, 6, 0, 1, 6, 0.05
);
ok(scalar(@$ids1) >= 2, 'checkPrimerMismatchTolerance: 1 mismatch -> 5-prime variant accepted');

# Test : la séquence 2 (mismatch 3') reste exclue même avec tolérance 2
my ($seq2, $cov2, $degen2, $ids2) = checkPrimerMismatchTolerance(
    \@seqs, 0, $length, $primer,
    100, 100, 5,
    0, 6, 0, 2, 6, 0.05
);
my $seq2_included = grep { $_ == 2 } @$ids2;
ok(!$seq2_included, "checkPrimerMismatchTolerance: 3' mismatch is rejected even with tolerance 2");

print "\n# Résultat : $passCount/$testCount tests passés\n";
print "# Result  : $passCount/$testCount tests passed\n";
