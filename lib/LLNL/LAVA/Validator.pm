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

package LLNL::LAVA::Validator;

use strict;
use warnings;
use vars qw(@ISA @EXPORT_OK);

require Exporter;
@ISA = qw(Exporter);

@EXPORT_OK = qw(
    checkPrimerMismatchTolerance
    isIUPACCompatible
    rev_comp
    generateIUPACCode
    getPrimerTargetedSequences
    validateCompleteSignatureSpacing
);


#=============================================================================
# IUPAC & SEQUENCE UTILITIES
#=============================================================================

sub rev_comp {
  my $seq = shift;
  my $rc = reverse($seq);
  $rc =~ tr/ACGTUacgtu/TGCAAtgcaa/;
  $rc =~ tr/RYSWKMBDHVNryswkmbdhvn/YRSWMKVHDBNyrswmkvhdbn/;
  return $rc;
}

sub isIUPACCompatible {
  my ($base, $iupac_code) = @_;
  
  my %iupac_map = (
    'A' => ['A'], 'T' => ['T'], 'G' => ['G'], 'C' => ['C'],
    'R' => ['A', 'G'], 'Y' => ['C', 'T'], 'S' => ['G', 'C'], 'W' => ['A', 'T'],
    'K' => ['G', 'T'], 'M' => ['A', 'C'],
    'B' => ['C', 'G', 'T'], 'D' => ['A', 'G', 'T'], 'H' => ['A', 'C', 'T'], 'V' => ['A', 'C', 'G'],
    'N' => ['A', 'C', 'G', 'T']
  );
  
  my $allowed_bases = $iupac_map{uc($iupac_code)} || [$iupac_code];
  my $target_base = uc($base);
  return grep { $_ eq $target_base } @{$allowed_bases};
}

sub generateIUPACCode {
    my ($bases_ref) = @_;
    return 'N' if !@$bases_ref;
    
    my $key = join('', sort @$bases_ref);
    my %table = (
        'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T',
        'AG'=>'R', 'CT'=>'Y', 'GT'=>'K', 'AC'=>'M', 'CG'=>'S', 'AT'=>'W',
        'ACG'=>'V', 'ACT'=>'H', 'AGT'=>'D', 'CGT'=>'B',
        'ACGT'=>'N'
    );
    return $table{$key} || 'N';
}

#=============================================================================
# HELPER: Get Primer Targeted Sequences
#=============================================================================
# Helper to check matching sequences based on IUPAC logic
sub getPrimerTargetedSequences {
  my ($sequences_r, $position, $length, $final_sequence) = @_;
  
  my $num_sequences = scalar(@{$sequences_r});
  my @targeted_sequences = ();
  
  for my $seq_idx (0 .. $num_sequences - 1) {
    my $sequence = $sequences_r->[$seq_idx];
    my $matches = 1;
    
    # Vérifier si ce primer matche cette séquence / Check if this primer matches this sequence
    for my $pos_offset (0 .. $length - 1) {
      my $abs_position = $position + $pos_offset;
      if ($abs_position >= length($sequence)) {
        $matches = 0;
        last;
      }
      
      my $seq_base = uc(substr($sequence, $abs_position, 1));
      my $primer_base = substr($final_sequence, $pos_offset, 1);
      
      # Vérifier compatibilité IUPAC / Check IUPAC compatibility
      if (!isIUPACCompatible($seq_base, $primer_base)) {
        $matches = 0;
        last;
      }
    }
    
    if ($matches) {
      push @targeted_sequences, $seq_idx;
    }
  }
  
  return @targeted_sequences;
}


#=============================================================================
# CORE VALIDATION LOGIC
#=============================================================================

sub checkPrimerMismatchTolerance {
  my ($sequences_r, $position, $length, $candidate_primer, $min_match_percent, $min_iupac_percent, $min_primer_coverage, $max_total_degen, $max_consec_degen, $max_3p_degen, $max_tolerated_mismatch, $zone_size, $min_base_freq) = @_;
  
  # Paramètres par défaut si non fournis / Default parameters if not provided
  $min_primer_coverage = 80 unless defined $min_primer_coverage;
  $min_base_freq = 0.05 unless defined $min_base_freq;
  
  my $num_sequences = scalar(@{$sequences_r});
  return ("", 0, 0, []) if $num_sequences == 0;
  
  my $candidate_primer_uc = uc($candidate_primer);
  
  # DEBUG LOGS (Condensed for shared usage)
  print "\n[LAVA::Validator] Checking Primer: $candidate_primer_uc @ Pos $position\n";
  
  # ========================================================================
  # PHASE 1: EXTRACTION DES RÉGIONS CIBLES (GAP-AWARE)
  # ========================================================================
  my @target_regions = ();
    for my $seq_idx (0 .. $num_sequences - 1) {
      my $sequence = $sequences_r->[$seq_idx];
    if ($position < length($sequence)) {
      my $region = substr($sequence, $position, $length);
      $region = uc($region);
      $region =~ s/-//g; # Séquence physique sans gaps
      
      if (length($region) == $length) {
         push @target_regions, { 
            seq_idx => $seq_idx, 
            region => $region 
         };
      }
    }
  }
  
  return ("", 0, 0, []) if @target_regions == 0;
  my $total_regions = scalar(@target_regions);
  
  # ========================================================================
  # PHASE 2: TEST D'ORIENTATION & COMPATIBILITÉ (SENSE vs ANTISENSE)
  # ========================================================================
  
  # Fonction locale
  my $test_orientation = sub {
      my ($prim, $targets_ref) = @_;
      my @perf = ();
      my @non = ();
      foreach my $t (@$targets_ref) {
          my $matches = 1;
          for(my $i=0; $i<$length; $i++) {
              my $pb = substr($prim, $i, 1);
              my $tb = substr($t->{region}, $i, 1);
              if (!isIUPACCompatible($tb, $pb)) {
                  $matches = 0;
                  last;
              }
          }
          if ($matches) { push @perf, $t; }
          else { push @non, $t; }
      }
      return (\@perf, \@non);
  };
  
  # 1. Test Sense
  my ($sense_perfect, $sense_non) = $test_orientation->($candidate_primer_uc, \@target_regions);
  my $sense_score = scalar(@$sense_perfect);
  
  # 2. Test Antisense (RC des cibles)
  my @target_regions_rc = ();
  foreach my $t (@target_regions) {
      push @target_regions_rc, { seq_idx => $t->{seq_idx}, region => rev_comp($t->{region}) };
  }
  my ($anti_perfect, $anti_non) = $test_orientation->($candidate_primer_uc, \@target_regions_rc);
  my $anti_score = scalar(@$anti_perfect);
  
  my @perfect_matches = ();
  my @non_matches = ();
  my $orientation = "SENSE";
  
  if ($anti_score > $sense_score) {
      @target_regions = @target_regions_rc;
      @perfect_matches = @$anti_perfect;
      @non_matches = @$anti_non;
      $orientation = "ANTISENSE";
  } else {
      @perfect_matches = @$sense_perfect;
      @non_matches = @$sense_non;
  }

  my $perfect_match_count = scalar(@perfect_matches);
  my $perfect_match_percent = ($perfect_match_count / $total_regions) * 100;
  
  
  # ========================================================================
  # PHASE 2b & PHASE 3: OPTIMISATION COMBINATOIRE OPTIMALE (Branch & Bound)
  # PHASE 2b & PHASE 3: OPTIMAL COMBINATORIAL OPTIMIZATION (Branch & Bound)
  # ========================================================================
  my $optimized_primer = $candidate_primer_uc;
  my $has_modifications = 0;
  my $three_prime_start_idx = $length - $zone_size;
  $three_prime_start_idx = 0 if $three_prime_start_idx < 0;
  my $three_prime_end_idx = $length - 1;
  
  # Fonction helper d'évaluation d'un candidat sur les cibles (Phase 4 intégrée)
  # Helper function to evaluate a candidate primer against targets (Integrated Phase 4)
  my $evaluate_candidate = sub {
    my ($test_primer_seq, $tolerance) = @_;
    $tolerance = 0 unless defined $tolerance;  # defaut = correspondance EXACTE (0 mismatch)
    my @compatible_seqs = ();
    
    for my $target (@target_regions) {
      my $is_fully_compatible = 1;
      my $mismatch_count = 0;
      
      for my $pos_offset (0 .. $length - 1) {
        my $target_base = substr($target->{region}, $pos_offset, 1);
        my $primer_base = substr($test_primer_seq, $pos_offset, 1);
        
        if (!isIUPACCompatible($target_base, $primer_base)) {
          # Vérifier si on est dans la zone 3' critique / Check if in the critical 3' zone
          if ($pos_offset >= $three_prime_start_idx && $pos_offset <= $three_prime_end_idx) {
            # Mismatch en 3' : REJET IMMÉDIAT / 3' mismatch: IMMEDIATE REJECT
            $is_fully_compatible = 0;
            last;
          } else {
            # Mismatch hors zone 3' : Toléré jusqu'à max_tolerated_mismatch / Mismatch outside 3' zone: Tolerated up to max_tolerated_mismatch
            $mismatch_count++;
            if ($mismatch_count > $tolerance) {
              $is_fully_compatible = 0;
              last;
            }
          }
        }
      }
      
      if ($is_fully_compatible) {
        push @compatible_seqs, $target->{seq_idx};
      }
    }
    return \@compatible_seqs;
  };

  my $best_primer = $candidate_primer_uc;
  my $best_compatible_seqs = $evaluate_candidate->($candidate_primer_uc);
  my $best_coverage_percent = (scalar(@$best_compatible_seqs) / $total_regions) * 100;
  my $best_mod_count = 0;
  my $best_3p_count = 0;

  # Si la couverture brute est inférieure à min_primer_coverage, on lance l'optimisation combinatoire
  # If raw coverage is below min_primer_coverage, run combinatorial optimization
  if ($best_coverage_percent < $min_primer_coverage && $max_total_degen > 0) {
      # 1. Identification de toutes les positions candidates modifiables
      # 1. Identify all candidate modifiable positions
      my @candidate_positions = ();
      
      for my $pos_offset (0 .. $length - 1) {
          my $primer_base = substr($candidate_primer_uc, $pos_offset, 1);
          my %base_counts = ();
          my $position_matches_count = 0;
          
          for my $target (@target_regions) {
              my $target_base = substr($target->{region}, $pos_offset, 1);
              $base_counts{$target_base}++;
              if ($target_base eq $primer_base) {
                  $position_matches_count++;
              }
          }
          
          my $primer_base_percent = ($position_matches_count / $total_regions) * 100;
          
          if ($primer_base_percent < $min_match_percent) {
              # FILTRE BRUIT (min_base_freq) / NOISE FILTER
              my $min_count_noise = $total_regions * $min_base_freq;
              my %significant_bases = ();
              $significant_bases{$primer_base} = 1;
              
              foreach my $b (keys %base_counts) {
                  if ($base_counts{$b} >= $min_count_noise) {
                      $significant_bases{$b} = 1;
                  }
              }
              my @all_bases = keys %significant_bases;
              my $iupac_code = generateIUPACCode(\@all_bases);
              
              if ($iupac_code ne 'N' && $iupac_code ne $primer_base) {
                  # Vérifier le pourcentage de compatibilité avec ce code IUPAC
                  my $iupac_matches_count = 0;
                  for my $target (@target_regions) {
                      my $target_base = substr($target->{region}, $pos_offset, 1);
                      if (isIUPACCompatible($target_base, $iupac_code)) {
                          $iupac_matches_count++;
                      }
                  }
                  my $iupac_percent = ($iupac_matches_count / $total_regions) * 100;
                  
                  if ($iupac_percent >= $min_iupac_percent) {
                      my $is_3p = ($pos_offset >= $three_prime_start_idx && $pos_offset <= $three_prime_end_idx) ? 1 : 0;
                      push @candidate_positions, {
                          pos   => $pos_offset,
                          code  => $iupac_code,
                          is_3p => $is_3p,
                          gain  => ($iupac_percent - $primer_base_percent)
                      };
                  }
              }
          }
      }

      # 2. Énumération combinatoire des sous-ensembles avec protection anti-explosion (Étape 1)
      # 2. Combinatorial subset enumeration with anti-explosion protection (Step 1)
      my $n_cand = scalar(@candidate_positions);
      if ($n_cand > 0) {
          # ÉTAPE 1 : Tri par gain de couverture décroissant et plafonnement aux Top-12 positions les plus impactantes
          # STEP 1: Sort by decreasing coverage gain and cap to Top-12 most impactful candidate positions
          @candidate_positions = sort { $b->{gain} <=> $a->{gain} || $a->{is_3p} <=> $b->{is_3p} } @candidate_positions;
          if ($n_cand > 12) {
              splice(@candidate_positions, 12);
              $n_cand = 12;
          }

          my $eval_count = 0;
          my $max_evaluations = 2000;
          my $recurse;
          $recurse = sub {
              my ($idx, $current_combo_ref, $count_3p) = @_;
              
              return if $eval_count >= $max_evaluations;
              
              # Évaluer si le sous-ensemble courant contient au moins une modification
              if (@$current_combo_ref > 0) {
                  $eval_count++;
                  my $test_primer = $candidate_primer_uc;
                  for my $item (@$current_combo_ref) {
                      substr($test_primer, $item->{pos}, 1) = $item->{code};
                  }
                  
                  my $seqs_ref = $evaluate_candidate->($test_primer);
                  my $cov_pct = (scalar(@$seqs_ref) / $total_regions) * 100;
                  my $mod_cnt = scalar(@$current_combo_ref);
                  
                  # Sélectionner si meilleure couverture, ou même couverture avec moins de modifications ou moins de dégénérescence en 3'
                  if ($cov_pct > $best_coverage_percent + 1e-6 ||
                      (abs($cov_pct - $best_coverage_percent) <= 1e-6 && $mod_cnt < $best_mod_count) ||
                      (abs($cov_pct - $best_coverage_percent) <= 1e-6 && $mod_cnt == $best_mod_count && $count_3p < $best_3p_count)) {
                      
                      $best_coverage_percent = $cov_pct;
                      $best_mod_count = $mod_cnt;
                      $best_3p_count = $count_3p;
                      $best_primer = $test_primer;
                      $best_compatible_seqs = $seqs_ref;
                  }
                  
                  # Si on atteint 100% de couverture, inutile d'ajouter d'autres dégénérescences sur cette branche
                  return if $best_coverage_percent >= 100.0 && $cov_pct >= 100.0;
              }
              
              # Arrêt si on atteint le budget max / Stop if max budget reached
              return if scalar(@$current_combo_ref) >= $max_total_degen;
              
              for my $i ($idx .. $n_cand - 1) {
                  return if $eval_count >= $max_evaluations;
                  my $item = $candidate_positions[$i];
                  my $new_count_3p = $count_3p + $item->{is_3p};
                  
                  # Vérification limite 3' / Check 3' limit
                  next if $new_count_3p > $max_3p_degen;
                  
                  # Vérification limite positions consécutives (indépendante de l'ordre de tri) / Check consecutive limit (order-independent)
                  if (@$current_combo_ref > 0) {
                      my @all_pos = map { $_->{pos} } @$current_combo_ref;
                      push @all_pos, $item->{pos};
                      @all_pos = sort { $a <=> $b } @all_pos;
                      
                      my $max_consec = 1;
                      my $curr_consec = 1;
                      for (my $j = 1; $j < scalar(@all_pos); $j++) {
                          if ($all_pos[$j] == $all_pos[$j-1] + 1) {
                              $curr_consec++;
                              $max_consec = $curr_consec if $curr_consec > $max_consec;
                          } else {
                              $curr_consec = 1;
                          }
                      }
                      next if $max_consec > $max_consec_degen;
                  }
                  
                  push @$current_combo_ref, $item;
                  $recurse->($i + 1, $current_combo_ref, $new_count_3p);
                  pop @$current_combo_ref;
              }
          };
          
          my @empty_combo = ();
          $recurse->(0, \@empty_combo, 0);
      }
  }

  $optimized_primer = $best_primer;
  $has_modifications = ($best_mod_count > 0) ? 1 : 0;

  # Validation finale : la tolerance aux mismatchs n'intervient qu'ICI, contre min_primer_coverage.
  # Final validation: mismatch tolerance applies ONLY here, against min_primer_coverage.
  my $tolerant_seqs_r = $evaluate_candidate->($optimized_primer, $max_tolerated_mismatch);
  my $final_coverage_percent = (scalar(@$tolerant_seqs_r) / $total_regions) * 100;

  if ($final_coverage_percent < $min_primer_coverage) {
    return ("", $final_coverage_percent, 0, []);
  }

  return ($optimized_primer, $final_coverage_percent, $has_modifications, $tolerant_seqs_r);
}


#=============================================================================
# SIGNATURE VALIDATION LOGIC
#=============================================================================

# Fonction pour valider que tous les primers d'une signature respectent l'espacement minimum
# et ne se chevauchent pas
sub validateCompleteSignatureSpacing
{
  my ($forwardPrimers_r, $reversePrimers_r, $minSpacing) = @_;
  
  # ─────────────────────────────────────────────────────────────────────────────
  # Helper interne : récupère le strand depuis un PrimerInfo ou un Oligo
  # Le strand est stocké sur l'Oligo (via setTag), PAS sur PrimerInfo directement.
  # Internal helper: get strand from a PrimerInfo or Oligo object.
  # The strand tag lives on the Oligo, NOT on PrimerInfo.
  # ─────────────────────────────────────────────────────────────────────────────
  my $get_strand = sub {
    my ($primer, $default_strand) = @_;
    # 1) Essayer via l'Oligo sous-jacent / Try via the underlying Oligo
    if ($primer->can('getAnalyzedPrimer')) {
      my $oligo = $primer->getAnalyzedPrimer();
      if (defined $oligo && $oligo->can('getTagExists') && $oligo->getTagExists('strand')) {
        return $oligo->getTag('strand');
      }
    }
    # 2) Essayer directement si c'est déjà un Oligo / Try directly if it's an Oligo
    if ($primer->can('getTagExists') && $primer->getTagExists('strand')) {
      return $primer->getTag('strand');
    }
    # 3) Fallback sur le nom du primer (F→plus, B→minus) / Fallback on primer name
    return $default_strand // 'plus';
  };

  # Créer une liste ordonnée de tous les primers avec leurs positions / Create ordered primer list
  my @allPrimers = ();
  
  # Ajouter les primers forward (F3, F2, F1, FSTEM) - strand par défaut : plus
  # Add forward primers (F3, F2, F1, FSTEM) - default strand: plus
  foreach my $primer (@{$forwardPrimers_r}) {
    next if (!defined $primer);
    my $location = $primer->getLocation();
    my $length   = $primer->getLength();
    my $strand   = $get_strand->($primer, 'plus');
    
    my ($start, $end);
    if ($strand eq 'minus') {
      # Primer sur strand minus : location = position 3' (extrémité droite)
      # Primer on minus strand: location = 3' end (rightmost position)
      $start = $location - $length + 1;
      $end   = $location;
    } else {
      # Primer sur strand plus : location = position 5' (extrémité gauche)
      # Primer on plus strand: location = 5' end (leftmost position)
      $start = $location;
      $end   = $location + $length - 1;
    }
    
    push @allPrimers, {
      'name'   => $primer->{name} // 'Forward',
      'start'  => $start,
      'end'    => $end,
      'length' => $length,
      'strand' => $strand,
    };
  }
  
  # Ajouter les primers reverse (BSTEM, B1, B2, B3) - strand par défaut : minus
  # Add reverse primers (BSTEM, B1, B2, B3) - default strand: minus
  foreach my $primer (@{$reversePrimers_r}) {
    next if (!defined $primer);
    my $location = $primer->getLocation();
    my $length   = $primer->getLength();
    my $strand   = $get_strand->($primer, 'minus');
    
    my ($start, $end);
    if ($strand eq 'minus') {
      # Reverse primer sur strand minus : location = position 3'
      # Reverse primer on minus strand: location = 3' end
      $start = $location - $length + 1;
      $end   = $location;
    } else {
      # Reverse primer sur strand plus : location = position 5'
      # Reverse primer on plus strand: location = 5' end
      $start = $location;
      $end   = $location + $length - 1;
    }
    
    push @allPrimers, {
      'name'   => $primer->{name} // 'Reverse',
      'start'  => $start,
      'end'    => $end,
      'length' => $length,
      'strand' => $strand,
    };
  }
  
  # Trier par position de début / Sort by start position
  @allPrimers = sort { $a->{'start'} <=> $b->{'start'} } @allPrimers;
  
  # Vérifier l'espacement entre primers consécutifs / Check spacing between consecutive primers
  for (my $i = 0; $i < @allPrimers - 1; $i++) {
    my $current = $allPrimers[$i];
    my $next    = $allPrimers[$i + 1];
    
    # Espacement = nb de bases entre la fin du primer courant et le début du suivant
    # Spacing = number of bases between end of current and start of next
    my $spacing = $next->{'start'} - $current->{'end'} - 1;
    
    # Chevauchement = espacement négatif → signature invalide
    # Overlap = negative spacing → invalid signature
    if ($spacing < 0) {
      return 0; # Échec de validation / Validation failed
    }
  }
  
  return 1; # Validation réussie / Validation passed
}
1;
