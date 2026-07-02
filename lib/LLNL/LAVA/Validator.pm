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
  # PHASE 2b: DÉCISION RAPIDE MODIFIÉE (SANS EARLY EXIT COMPLET)
  # ========================================================================
  my $optimized_primer = $candidate_primer_uc;
  my $has_modifications = 0;
  my ($three_prime_start_idx, $three_prime_end_idx);
  
  # Si la couverture brute est déjà suffisante, on évite l'optimisation dégénérée (Phase 3)
  # mais on passe quand même par la Phase 4 pour évaluer la tolérance aux mismatches réels
  # If raw coverage is sufficient, skip degeneracy optimization (Phase 3)
  # but proceed to Phase 4 to compute actual mismatch tolerance.
  if ($perfect_match_percent >= $min_primer_coverage) {
      # Skip Phase 3
  } else {
      # ========================================================================
      # PHASE 3: ANALYSE POSITION PAR POSITION (Avec contraintes de dégénérescence) / PHASE 3: POSITION BY POSITION ANALYSIS (With degeneracy constraints)
      # ========================================================================
      
      my @position_compatible_seqs = ();
      my $degen_total = 0;
      my $degen_consec = 0;
      my $degen_3p = 0;
      
      $three_prime_start_idx = $length - $zone_size;
      $three_prime_start_idx = 0 if $three_prime_start_idx < 0;
      $three_prime_end_idx = $length - 1;
      
      # Initialiser tracking
      for my $pos_offset (0 .. $length - 1) {
        $position_compatible_seqs[$pos_offset] = [ map { $_->{seq_idx} } @perfect_matches ];
      }
      
      for my $pos_offset (0 .. $length - 1) {
        my $primer_base = substr($candidate_primer_uc, $pos_offset, 1);
        my %base_counts = ();
        my @position_matches = @{$position_compatible_seqs[$pos_offset]};
        
        for my $target (@non_matches) {
          my $target_base = substr($target->{region}, $pos_offset, 1);
          $base_counts{$target_base}++;
          if ($target_base eq $primer_base) {
            push @position_matches, $target->{seq_idx};
          }
        }
        
        my $primer_base_count = $base_counts{$primer_base} || 0;
        my $total_primer_matches = scalar(@position_matches);
        my $primer_base_percent = ($total_primer_matches / $total_regions) * 100;
        
        if ($primer_base_percent < $min_match_percent) {
          # Besoin d'une base dégénérée / Need a degenerate base
          
          # Vérifier les limites de dégénérescence avant de générer le code / Check degeneracy limits before generating code
          $degen_total++;
          $degen_consec++;
          
          if ($pos_offset >= $three_prime_start_idx && $pos_offset <= $three_prime_end_idx) {
            $degen_3p++;
          }
          
          if ($degen_total > $max_total_degen || 
              $degen_consec > $max_consec_degen || 
              $degen_3p > $max_3p_degen) {
            # Limite de dégénérescence dépassée, on réinitialise l'amorce à sa forme brute et on arrête l'optimisation
            $optimized_primer = $candidate_primer_uc;
            $has_modifications = 0;
            last;
          }
          
          # FILTRE BRUIT (min_base_freq)
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
          
          if ($iupac_code eq 'N') {
            # Trop de variation, pas de modif possible fiable
            $optimized_primer = $candidate_primer_uc;
            $has_modifications = 0;
            last;
          }
          
          my @iupac_matches = @{$position_compatible_seqs[$pos_offset]};
          for my $target (@non_matches) {
            my $target_base = substr($target->{region}, $pos_offset, 1);
            if (isIUPACCompatible($target_base, $iupac_code)) {
              push @iupac_matches, $target->{seq_idx} unless grep { $_ == $target->{seq_idx} } @iupac_matches;
            }
          }
          
          my $iupac_percent = (scalar(@iupac_matches) / $total_regions) * 100;
          
          if ($iupac_percent < $min_iupac_percent) {
            $optimized_primer = $candidate_primer_uc;
            $has_modifications = 0;
            last;
          }
          
          substr($optimized_primer, $pos_offset, 1) = $iupac_code;
          $position_compatible_seqs[$pos_offset] = \@iupac_matches;
          $has_modifications = 1;
        } else {
          # Pas de dégénérescence à cette position / No degeneracy at this position
          $degen_consec = 0;
          $position_compatible_seqs[$pos_offset] = \@position_matches;
        }
      }
  }
  
  # ========================================================================
  # PHASE 4: VALIDATION FINALE (Tolérance Mismatches & Protection 3') / PHASE 4: FINAL VALIDATION (Mismatch Tolerance & 3' Protection)
  # ========================================================================
  my @final_compatible_sequences = ();
  
  # L'extrémité 3' est TOUJOURS à la fin de la chaîne (5' -> 3') / The 3' end is ALWAYS at the end of the chain (5' -> 3')
  $three_prime_start_idx = $length - $zone_size;
  $three_prime_start_idx = 0 if $three_prime_start_idx < 0;
  $three_prime_end_idx = $length - 1;
  
  for my $target (@target_regions) {
    my $is_fully_compatible = 1;
    my $mismatch_count = 0;
    
    for my $pos_offset (0 .. $length - 1) {
      my $target_base = substr($target->{region}, $pos_offset, 1);
      my $primer_base = substr($optimized_primer, $pos_offset, 1);
      
      if (!isIUPACCompatible($target_base, $primer_base)) {
        # Vérifier si on est dans la zone 3' critique / Check if in the critical 3' zone
        if ($pos_offset >= $three_prime_start_idx && $pos_offset <= $three_prime_end_idx) {
          # Mismatch en 3' : REJET IMMÉDIAT
          $is_fully_compatible = 0;
          last;
        } else {
          # Mismatch hors zone 3' : Toléré jusqu'à max_tolerated_mismatch / Mismatch outside 3' zone: Tolerated up to max_tolerated_mismatch
          $mismatch_count++;
          if ($mismatch_count > $max_tolerated_mismatch) {
            $is_fully_compatible = 0;
            last;
          }
        }
      }
    }
    
    if ($is_fully_compatible) {
      push @final_compatible_sequences, $target->{seq_idx};
    }
  }
  
  my $final_coverage_percent = (scalar(@final_compatible_sequences) / $total_regions) * 100;

  if ($final_coverage_percent < $min_primer_coverage) {
    return ("", $final_coverage_percent, 0, []);
  }
  
  return ($optimized_primer, $final_coverage_percent, $has_modifications, \@final_compatible_sequences);
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
  
  # Ajouter les primers forward (F3, F2, F1, FSTEM) — strand par défaut : plus
  # Add forward primers (F3, F2, F1, FSTEM) — default strand: plus
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
  
  # Ajouter les primers reverse (BSTEM, B1, B2, B3) — strand par défaut : minus
  # Add reverse primers (BSTEM, B1, B2, B3) — default strand: minus
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
