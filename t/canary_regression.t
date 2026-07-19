#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;

# Subroutine to get the longest common substring
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

# Add lib to PERL5LIB
$ENV{PERL5LIB} = "./lib:" . ($ENV{PERL5LIB} || "");

# Parameters common to both
my $common_args = "--alignment_fasta t/fixtures/rota_canary_40.fasta --threads 1 --signature_max_length 206 " .
    "--dna_conc 400.0 --dntp_conc 1.4 --entropy_threshold 1.5 " .
    "--inner_primer_min_length 15 --inner_primer_target_length 18 --inner_primer_max_length 22 " .
    "--inner_primer_min_tm 59.0 --inner_primer_target_tm 60.0 --inner_primer_max_tm 65.0 " .
    "--middle_primer_min_length 15 --middle_primer_target_length 18 --middle_primer_max_length 22 " .
    "--middle_primer_min_tm 59.0 --middle_primer_target_tm 60.0 --middle_primer_max_tm 65.0 " .
    "--outer_primer_min_length 15 --outer_primer_target_length 18 --outer_primer_max_length 22 " .
    "--outer_primer_min_tm 57.0 --outer_primer_target_tm 58.0 --outer_primer_max_tm 59.0 " .
    "--max_dist_middle_inner 30 --max_dist_outer_middle 30 --min_primer_coverage 70.0 " .
    "--min_signatures_for_success 1 --resolve_overlap_by coverage " .
    "--max_3prime_degenerate_bases 0 --max_consecutive_degenerate_bases 2 " .
    "--max_total_degenerate_bases 3 --max_poly_bases 5 --max_tolerated_mismatches 2 " .
    "--max_tm_diff 5.0 --primer_iupac_min_percent 80.0 --primer_min_match_percent 70.0 " .
    "--min_base_frequency 0.2 --three_prime_zone_size 3 --salt_monovalent 50.0 " .
    "--salt_divalent 8.0 --penalty_plateau 0.25 --penalty_slope 0.15 " .
    "--max_overlap_percent 0.0 --max_per_window 0 --max_primer_gen 10000.0 --window_size 0";

# Run LOOP test
my $loop_args = $common_args . " --output_file t/canary_loop --loop_min_gap 20 " .
    "--loop_primer_min_length 15 --loop_primer_target_length 18 --loop_primer_max_length 22 " .
    "--loop_primer_min_tm 59.0 --loop_primer_target_tm 60.0 --loop_primer_max_tm 61.0";
    
diag("Running LOOP canary test...");
my $loop_exit = system("perl lava_loop_primer.pl $loop_args > /dev/null 2>&1");
is($loop_exit, 0, "lava_loop_primer.pl executes successfully");

my @loop_sigs = glob("t/canary_loop_signatures_individuelles/signature_*_VALID_*.txt");
ok(scalar(@loop_sigs) > 0, "REGRESSION LOOP : 0 signature sur la reference rota");

# Run STEM test
my $stem_args = $common_args . " --output_file t/canary_stem --include_stem_primers 1 " .
    "--stem_primer_min_length 15 --stem_primer_target_length 18 --stem_primer_max_length 22 " .
    "--stem_primer_min_tm 59.0 --stem_primer_target_tm 60.0 --stem_primer_max_tm 61.0";

diag("Running STEM canary test...");
my $stem_exit = system("perl lava_stem_primer.pl $stem_args > /dev/null 2>&1");
is($stem_exit, 0, "lava_stem_primer.pl executes successfully");

my @stem_sigs = glob("t/canary_stem_signatures_individuelles/signature_*_VALID_*.txt");
ok(scalar(@stem_sigs) > 0, "REGRESSION STEM : 0 signature sur la reference rota");

# Anti-dimer check
if (scalar(@stem_sigs) > 0) {
    my $sig_file = $stem_sigs[0];
    open my $fh, '<', $sig_file or die "Cannot open $sig_file: $!";
    my ($fstem, $bstem) = ("", "");
    while (<$fh>) {
        if (/^# FSTEM:\s+([A-Za-z]+)/) {
            $fstem = $1;
        } elsif (/^# BSTEM:\s+([A-Za-z]+)/) {
            $bstem = $1;
        }
    }
    close $fh;
    
    ok($fstem ne "", "FSTEM sequence found in signature");
    ok($bstem ne "", "BSTEM sequence found in signature");
    
    if ($fstem ne "" && $bstem ne "") {
        my $rc_bstem = reverse_complement($bstem);
        my $dimer_len = longest_common_substring($fstem, $rc_bstem);
        
        cmp_ok($dimer_len, '<', 8, "REGRESSION STEM : dimere FSTEM/BSTEM detecte (complementarite = $dimer_len nt)");
    }
}

done_testing();
