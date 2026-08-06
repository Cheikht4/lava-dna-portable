#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Compare qw(compare);
use File::Copy    qw(copy);

# =============================================================================
# canary_regression.t -- LAVA-DNA v2026
# Tests de regression + equivalence baseline pour lava_loop_primer.pl
# et lava_stem_primer.pl.
#
# Variables d'environnement :
#   LAVA_UPDATE_BASELINE=1  -> regenere les baselines au lieu de comparer.
#                              A utiliser UNIQUEMENT apres validation humaine
#                              d'une nouvelle reference voulue.
#   PERL5LIB                -> ajoute lib/ automatiquement si non defini.
# =============================================================================

# --- Subroutines utilitaires ---

sub longest_common_substring {
    my ($s1, $s2) = @_;
    my $max_len = 0;
    for my $i (0 .. length($s1) - 1) {
        for my $j (1 .. length($s1) - $i) {
            my $sub = substr($s1, $i, $j);
            if (index($s2, $sub) != -1) {
                $max_len = length($sub) if length($sub) > $max_len;
            }
        }
    }
    return $max_len;
}

sub reverse_complement {
    my ($seq) = @_;
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

# Compte le nombre de sequences dans un FASTA
sub count_fasta_sequences {
    my ($file) = @_;
    return 0 unless -f $file;
    open my $fh, '<', $file or return 0;
    my $count = 0;
    $count++ while <$fh> =~ /^>/;
    close $fh;
    return $count;
}

# Compare deux fichiers texte ligne a ligne, retourne undef si identiques
# ou un message d'erreur avec les premieres lignes divergentes.
sub strict_diff {
    my ($got, $expected) = @_;
    return "Fichier genere absent : $got"       unless -f $got;
    return "Fichier de reference absent : $expected" unless -f $expected;

    open my $fh_g, '<', $got      or return "Impossible de lire $got : $!";
    open my $fh_e, '<', $expected or return "Impossible de lire $expected : $!";

    my $line_no = 0;
    while (1) {
        my $lg = <$fh_g>;
        my $le = <$fh_e>;
        last unless defined $lg || defined $le;
        $line_no++;
        $lg //= '';
        $le //= '';
        if ($lg ne $le) {
            chomp $lg; chomp $le;
            return "Divergence ligne $line_no :\n  OBTENU   : $lg\n  ATTENDU  : $le";
        }
    }
    close $fh_g; close $fh_e;
    return undef;  # identiques
}

# Verifie ou regenere un fichier de baseline
sub check_or_update_baseline {
    my ($generated, $baseline, $label) = @_;

    if ($ENV{LAVA_UPDATE_BASELINE}) {
        copy($generated, $baseline)
            or die "Impossible de copier $generated -> $baseline : $!";
        diag("BASELINE MISE A JOUR : $baseline");
        pass("$label (baseline regeneree)");
    } else {
        my $err = strict_diff($generated, $baseline);
        if (!defined $err) {
            pass($label);
        } else {
            fail($label);
            diag($err);
        }
    }
}

# =============================================================================
# Initialisation
# =============================================================================

$ENV{PERL5LIB} = "./lib:" . ($ENV{PERL5LIB} || "");

my $UPDATE_BASELINE = $ENV{LAVA_UPDATE_BASELINE} // 0;

if ($UPDATE_BASELINE) {
    diag("MODE REGENERATION BASELINE actif (LAVA_UPDATE_BASELINE=1)");
}

# Parametres communs aux deux scripts
my $common_args =
    "--alignment_fasta t/fixtures/rota_canary_40.fasta --threads 1 " .
    "--signature_max_length 206 " .
    "--dna_conc 400.0 --dntp_conc 1.4 --entropy_threshold 1.5 " .
    "--inner_primer_min_length 15 --inner_primer_target_length 18 --inner_primer_max_length 22 " .
    "--inner_primer_min_tm 59.0 --inner_primer_target_tm 60.0 --inner_primer_max_tm 65.0 " .
    "--middle_primer_min_length 15 --middle_primer_target_length 18 --middle_primer_max_length 22 " .
    "--middle_primer_min_tm 59.0 --middle_primer_target_tm 60.0 --middle_primer_max_tm 65.0 " .
    "--outer_primer_min_length 15 --outer_primer_target_length 18 --outer_primer_max_length 22 " .
    "--outer_primer_min_tm 57.0 --outer_primer_target_tm 58.0 --outer_primer_max_tm 59.0 " .
    "--min_primer_coverage 70.0 " .
    "--resolve_overlap_by coverage " .
    "--max_3prime_degenerate_bases 0 --max_consecutive_degenerate_bases 2 " .
    "--max_total_degenerate_bases 3 --max_poly_bases 5 --max_tolerated_mismatches 2 " .
    "--max_tm_diff 5.0 --primer_iupac_min_percent 80.0 --primer_min_match_percent 70.0 " .
    "--min_base_frequency 0.2 --three_prime_zone_size 3 --salt_monovalent 50.0 " .
    "--salt_divalent 8.0 " .
    "--max_overlap_percent 0.0 --max_per_window 0 --max_primer_gen 10000.0 --window_size 0";

# =============================================================================
# PARTIE 1 : LOOP
# =============================================================================

my $loop_out  = "t/canary_loop";
my $loop_args = $common_args .
    " --output_file $loop_out --loop_min_gap 20 " .
    "--loop_primer_min_length 15 --loop_primer_target_length 18 --loop_primer_max_length 22 " .
    "--loop_primer_min_tm 59.0 --loop_primer_target_tm 60.0 --loop_primer_max_tm 61.0";

diag("=== LOOP canary : execution (--threads 1) ===");
my $loop_exit = system("perl lava_loop_primer.pl $loop_args > /dev/null 2>&1");
is($loop_exit, 0, "lava_loop_primer.pl executes successfully");

my @loop_sigs = glob("t/canary_loop_signatures_individuelles/signature_*_VALID_*.txt");
ok(scalar(@loop_sigs) > 0, "REGRESSION LOOP : 0 signature sur la reference rota");

# --- Tests d'equivalence baseline LOOP ---
diag("=== LOOP : comparaison avec baseline ===");

for my $ext (qw(.primers .all_signatures .dash _amplified.fasta)) {
    my $generated = "${loop_out}${ext}";
    my $baseline  = "t/baseline/canary_loop${ext}";
    check_or_update_baseline($generated, $baseline,
        "BASELINE LOOP${ext} : sortie identique a la reference");
}

# =============================================================================
# PARTIE 2 : STEM
# =============================================================================

my $stem_out  = "t/canary_stem";
my $stem_args = $common_args .
    " --output_file $stem_out --include_stem_primers 1 " .
    "--stem_primer_min_length 15 --stem_primer_target_length 18 --stem_primer_max_length 22 " .
    "--stem_primer_min_tm 59.0 --stem_primer_target_tm 60.0 --stem_primer_max_tm 61.0";

diag("=== STEM canary : execution (--threads 1) ===");
my $stem_exit = system("perl lava_stem_primer.pl $stem_args > /dev/null 2>&1");
is($stem_exit, 0, "lava_stem_primer.pl executes successfully");

my @stem_sigs = glob("t/canary_stem_signatures_individuelles/signature_*_VALID_*.txt");
ok(scalar(@stem_sigs) > 0, "REGRESSION STEM : 0 signature sur la reference rota");

# --- Tests d'equivalence baseline STEM ---
diag("=== STEM : comparaison avec baseline ===");

for my $ext (qw(.primers .all_signatures .dash _amplified.fasta)) {
    my $generated = "${stem_out}${ext}";
    my $baseline  = "t/baseline/canary_stem${ext}";
    check_or_update_baseline($generated, $baseline,
        "BASELINE STEM${ext} : sortie identique a la reference");
}

# =============================================================================
# PARTIE 3 : Anti-dimer STEM
# =============================================================================

diag("=== STEM : verification anti-dimere ===");

if (scalar(@stem_sigs) > 0) {
    my $sig_file = $stem_sigs[0];
    open my $fh, '<', $sig_file or die "Cannot open $sig_file: $!";
    my ($fstem, $bstem) = ("", "");
    while (<$fh>) {
        if (/^# FSTEM:\s+([A-Za-z]+)/) { $fstem = $1; }
        elsif (/^# BSTEM:\s+([A-Za-z]+)/) { $bstem = $1; }
    }
    close $fh;

    ok($fstem ne "", "FSTEM sequence found in signature");
    ok($bstem ne "", "BSTEM sequence found in signature");

    if ($fstem ne "" && $bstem ne "") {
        my $rc_bstem  = reverse_complement($bstem);
        my $dimer_len = longest_common_substring($fstem, $rc_bstem);
        cmp_ok($dimer_len, '<', 8,
            "REGRESSION STEM : dimere FSTEM/BSTEM detecte (complementarite = $dimer_len nt)");
    }
}

# =============================================================================
# PARTIE 4 : Reproductibilite stricte
# =============================================================================

diag("=== DETERMINISME : verification de reproductibilite (RUN 2) ===");
my $loop_out_run2 = "t/canary_loop_run2";
my $loop_args_run2 = $common_args .
    " --output_file $loop_out_run2 --loop_min_gap 20 " .
    "--loop_primer_min_length 15 --loop_primer_target_length 18 --loop_primer_max_length 22 " .
    "--loop_primer_min_tm 59.0 --loop_primer_target_tm 60.0 --loop_primer_max_tm 61.0";

my $exit2 = system("perl lava_loop_primer.pl $loop_args_run2 > /dev/null 2>&1");
is($exit2, 0, "lava_loop_primer.pl (RUN 2) executes successfully");

for my $ext (qw(.primers .all_signatures .dash _amplified.fasta)) {
    my $run1 = "${loop_out}${ext}";
    my $run2 = "${loop_out_run2}${ext}";
    my $err = strict_diff($run2, $run1);
    if (!defined $err) {
        pass("DETERMINISME LOOP${ext} : le deuxieme run est strictement identique au premier");
    } else {
        fail("DETERMINISME LOOP${ext} : le deuxieme run differe du premier");
        diag($err);
    }
}

done_testing();
