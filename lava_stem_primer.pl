#!/usr/bin/perl -w

################################################################################
#
#
# Version 0.1.2 (2016)
# Updated by Michaël Bekaert <michael.bekaert@stir.ac.uk>.
# Produced at the Institute of Aquacuture, University of Stirling, UK
#
# MODIFIED VERSION - STEM PRIMER ARCHITECTURE:
# 5' ←── Signature LAMP (Longueur Dynamique 'L') ──→ 3'
# F3 ──(≥1)── F2 ──(≥1)── F1──(≥1)──  FL════(≥1)════BL ──(≥1)── B1 ──(≥1)── B2 ──(≥1)── B3
# │                       └─     Inner+STEM Pair (~40% de L)  ─┘                       │
# │           └──────              Middle Pair (~76% de L)          ──────┘             │
# └──────────────────────────────── Outer Pair (100% de L) ──────────────────────────────┘
# NOTE: STEM primers (FL/BL) are placed BETWEEN inner primers, not between inner/middle
#
# Copyright (c) 2010, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Clinton Torres <clinton.torres@llnl.gov>.
# CODE-42036.
# All rights reserved.
#
# This file is part of LAVA (LAMP Assay Versatile Analysis). For details, 
# see http://code.google.com/p/lava-dna/ . 
# Please also read the Additional BSD Notice.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# . Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the disclaimer below.
# . Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the disclaimer (as noted below) in the 
#   documentation and/or other materials provided with the distribution.
# . Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#   used to endorse or promote products derived from this software without 
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, 
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Additional BSD Notice
# 1. This notice is required to be provided under our contract with the 
#    U.S. Department of Energy (DOE). This work was produced at Lawrence 
#    Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with 
#    the DOE.
# 2. Neither the United States Government nor Lawrence Livermore National 
#    Security, LLC nor any of their employees, makes any warranty, express or 
#    implied, or assumes any liability or responsibility for the accuracy, 
#    completeness, or usefulness of any information, apparatus, product, or 
#    process disclosed, or represents that its use would not infringe 
#    privately-owned rights.
# 3. Also, reference herein to any specific commercial products, process, or 
#    services by trade name, trademark, manufacturer or otherwise does not 
#    necessarily constitute or imply its endorsement, recommendation, or 
#    favoring by the United States Government or Lawrence Livermore National 
#    Security, LLC. The views and opinions of authors expressed herein do not 
#    necessarily state or reflect those of the United States Government or 
#    Lawrence Livermore National Security, LLC, and shall not be used for 
#    advertising or product endorsement purposes.
################################################################################

use strict;
use Time::HiRes qw(time);
use warnings;
use Carp;
use lib 'lib';

use Getopt::Long;


use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::LocatableSeq;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::OligoEnumerator::Primer3Conserved;

use LLNL::LAVA::PrimerAnalyzer::PCRPrimer;
use LLNL::LAVA::PrimerInfo;
use LLNL::LAVA::PrimerSet::PCRPair;

use LLNL::LAVA::PrimerSetAnalyzer::PCRPair;
use LLNL::LAVA::PrimerSetInfo::PCRPair;

use LLNL::LAVA::PrimerSet::LAMP;
use LLNL::LAVA::Core qw(generateDistancePenalties calculate_proportional_geometry countDegenerateBases);
use LLNL::LAVA::Validator qw(checkPrimerMismatchTolerance getPrimerTargetedSequences isIUPACCompatible rev_comp generateIUPACCode validateCompleteSignatureSpacing);
use LLNL::LAVA::PipelineUtils qw(buildNativeReversePool analyzeAll enumeratePairs buildMetricsArray reducePairInfosByPenalty reducePrimersByOverlap reduceSignaturesByOverlap flattenInfoData buildBigMerge calculateSignatureIntersection createPerSignatureFiles createAmplificationFiles analyzeSignatureCombinations generateCombinations calculateDynamicPairLengths); # buildReversePrimers retiré (DEPRECATED, remplacé par buildNativeReversePool)

# Activer l'auto-flush de STDOUT pour les logs temps réel via Flask / Enable STDOUT auto-flush for real-time logs via Flask
# Enable STDOUT autoflush for real-time log streaming via Flask
$| = 1;
# Autoflush STDERR pour que \r fonctionne en temps reel (comme tqdm)
# Autoflush STDERR so \r works in real-time (like tqdm)
use IO::Handle;
STDERR->autoflush(1);
# Detection du terminal interactif : barre en place si TTY, silencieuse si fichier
# Detect interactive terminal: in-place bar if TTY, silent if redirected to file
our $_LAVA_IS_TTY = -t STDERR ? 1 : 0;


################################################################################
# FONCTIONS DE VALIDATION ET D'ANALYSE DES SIGNATURES
# Ces fonctions sont désormais dans LLNL::LAVA::PipelineUtils (Phase 36). / These functions are now in LLNL::LAVA::PipelineUtils (Phase 36).
# Functions now in LLNL::LAVA::PipelineUtils (Phase 36 harmonization):
#   - calculateSignatureIntersection
#   - analyzeSignatureCombinations
#   - generateCombinations
#   - createPerSignatureFiles
#   - createAmplificationFiles
#   - calculateDynamicPairLengths
################################################################################

=head2 getOligosWithMismatchTolerance

Version améliorée de getOligos qui intègre la tolérance aux mismatches. / Improved version of getOligos that integrates mismatch tolerance.
Utilise d'abord Primer3 sur la première séquence, puis applique notre logique / First uses Primer3 on the first sequence, then applies our logic 
de tolérance aux mismatches pour valider et modifier les amorces candidates. / for mismatch tolerance to validate and modify candidate primers.

=cut

sub getOligosWithMismatchTolerance {
  my ($enumerator, $alignment, $min_match_percent, $min_iupac_percent, $min_primer_coverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, $label) = @_;
  
  $label //= "Forward";
  my $progress_label = "Validation $label";
  
  my $sequenceCount = $alignment->num_sequences();
  if ($sequenceCount <= 0) {
    confess("data error - MSA must contain at least one sequence");
  }
  
  print "INFO: Utilisation de la tolérance aux mismatches activée\n";
  print "  - Seuil de concordance stricte / Strict match threshold: ${min_match_percent}%\n";
  print "  - Seuil de couverture IUPAC / IUPAC coverage threshold: ${min_iupac_percent}%\n";
  
  # Extraire toutes les séquences du MSA / Extract all sequences from the MSA
  my @sequences = ();
  foreach my $sequence ($alignment->each_seq()) {
    my $seqContent = $sequence->seq();
    $seqContent = uc($seqContent);  # Convertir en majuscules d'abord / Convert to uppercase first
    $seqContent =~ s/[^ATCG]/N/g;  # Puis normaliser (remplacer caractères non-ADN par N) / Then normalize (replace non-DNA characters with N)
    push @sequences, $seqContent;
  }
  
  # Utiliser l'enumerator original pour obtenir les amorces candidates / Use the original enumerator to get candidate primers
  my @candidatePrimers = $enumerator->getOligos($alignment);
  my @validatedPrimers = ();
  my $strict_count = 0;
  my $degenerate_count = 0;
  my $rejected_count = 0;
  
  my $nb_fwd_candidates = scalar(@candidatePrimers);
  print "INFO: Analyse de $nb_fwd_candidates amorces candidates Forward avec tolerance aux mismatches...\n";

  # Barre de progression / Progress bar (auto-detect Term::ProgressBar ou ASCII fallback)
  my $_has_pb = eval { require Term::ProgressBar; 1 } || 0;
  my $_pb_obj = undef;
  my $_pb_t0  = time();
  if ($_has_pb && -t STDOUT) {
    $_pb_obj = Term::ProgressBar->new({
      name   => $progress_label,
      count  => $nb_fwd_candidates,
      ETA    => 'linear',
      remove => 0,
      fh     => \*STDERR,
    });
    $_pb_obj->minor(0);
  }
  my $_pb_done = 0;

  foreach my $primer (@candidatePrimers) {
    my $location = $primer->location();
    my $length = $primer->length();
    my $original_sequence = $primer->sequence();
    
    # Appliquer notre analyse de tolérance aux mismatches / Apply our mismatch tolerance analysis
    my ($final_sequence, $coverage_percent, $is_degenerate, $compatible_seq_ids) = 
      checkPrimerMismatchTolerance(\@sequences, $location, $length, $original_sequence,
                                  $min_match_percent, $min_iupac_percent, $min_primer_coverage,
                                  $maxTotalDegen, $maxConsecDegen, 
                                  $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);
    
    # Utiliser le même seuil que l'algorithme interne (paramètre passé à la fonction) / Use the same threshold as the internal algorithm (parameter passed to the function)
    my $min_primer_acceptance = $min_primer_coverage;
    if ($coverage_percent >= $min_primer_acceptance) {
      # Amorce acceptée - créer la version finale / Primer accepted - create the final version
      my $validatedPrimer = $primer->clone();
      
      if ($is_degenerate) {
        # Amorce dégénérée - modifier la séquence / Degenerate primer - modify the sequence
        $validatedPrimer->sequence($final_sequence);
        $validatedPrimer->setTag("is_degenerate", 1);
        $validatedPrimer->setTag("original_sequence", $original_sequence);
        $validatedPrimer->setTag("iupac_coverage", sprintf("%.1f", $coverage_percent));
        $validatedPrimer->setTag("compatible_sequence_ids", $compatible_seq_ids);
        $degenerate_count++;
        
        print "DEGENERATE PRIMER acceptée - Pos: $location, Couv: " . 
              sprintf("%.1f", $coverage_percent) . "%, Seq: $final_sequence\n";
        # Mise à jour barre / Update progress bar
        $_pb_done++;
        if ($_has_pb && $_pb_obj) { $_pb_obj->update($_pb_done); }
        elsif ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
          # Ligne de progression structuree pour Flask / Structured progress line for Flask
          if ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
            my $pct = int($_pb_done/$nb_fwd_candidates*100);
            my $eta = ($_pb_done > 0 && $_pb_done < $nb_fwd_candidates)
                      ? int(($nb_fwd_candidates-$_pb_done)/($_pb_done/(time()-$_pb_t0+0.001)))
                      : 0;
            my $rate = $_pb_done / (time()-$_pb_t0+0.001);
            printf("[LAVA-PROGRESS] %s|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                   $progress_label,$_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
          }
        }
      } else {
        # Amorce stricte - garder la séquence originale / Strict primer - keep the original sequence
        $validatedPrimer->setTag("is_degenerate", 0);
        $validatedPrimer->setTag("iupac_coverage", "100.0");
        $validatedPrimer->setTag("compatible_sequence_ids", $compatible_seq_ids);
        $strict_count++;
        # Mise à jour barre / Update progress bar
        $_pb_done++;
        if ($_has_pb && $_pb_obj) { $_pb_obj->update($_pb_done); }
        elsif ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
          # Ligne de progression structuree pour Flask / Structured progress line for Flask
          if ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
            my $pct = int($_pb_done/$nb_fwd_candidates*100);
            my $eta = ($_pb_done > 0 && $_pb_done < $nb_fwd_candidates)
                      ? int(($nb_fwd_candidates-$_pb_done)/($_pb_done/(time()-$_pb_t0+0.001)))
                      : 0;
            my $rate = $_pb_done / (time()-$_pb_t0+0.001);
            printf("[LAVA-PROGRESS] %s|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                   $progress_label,$_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
          }
        }
      }
      
      push @validatedPrimers, $validatedPrimer;
    } else {
      # Amorce rejetée / Primer rejected
      $rejected_count++;
      print "REJECTED PRIMER - Pos: $location, Couv: " . 
            sprintf("%.1f", $coverage_percent) . "% < ${min_primer_acceptance}%\n";
      # Mise à jour barre / Update progress bar
      $_pb_done++;
      if ($_has_pb && $_pb_obj) { $_pb_obj->update($_pb_done); }
      elsif ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
        # Ligne de progression structuree pour Flask / Structured progress line for Flask
        if ($_pb_done % 200 == 0 || $_pb_done == $nb_fwd_candidates) {
          my $pct = int($_pb_done/$nb_fwd_candidates*100);
          my $eta = ($_pb_done > 0 && $_pb_done < $nb_fwd_candidates)
                    ? int(($nb_fwd_candidates-$_pb_done)/($_pb_done/(time()-$_pb_t0+0.001)))
                    : 0;
          my $rate = $_pb_done / (time()-$_pb_t0+0.001);
          printf("[LAVA-PROGRESS] %s|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                 $progress_label,$_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
        }
      }
    }
  }
  
  # Finaliser la barre / Finalize progress bar
  if ($_has_pb && $_pb_obj) { $_pb_obj->update($nb_fwd_candidates); }
  elsif ($_LAVA_IS_TTY) {
    # Effacer la barre et passer a la ligne suivante / Clear bar and move to next line
    printf(STDERR "\r%-80s\n", "");
  }
  print "RÉSULTATS tolérance mismatches:\n";
  print "  - Amorces strictes acceptées / accepted: $strict_count\n";
  print "  - Amorces dégénérées acceptées / accepted: $degenerate_count\n";
  print "  - Amorces rejetées: $rejected_count\n";
  print "  - Total validé / Total validated: " . scalar(@validatedPrimers) . "/" . scalar(@candidatePrimers) . "\n\n";
  
  return @validatedPrimers;
}

################################################################################

{ # Fake main() to enforce scope
  my %options;
  my %optionMap =
    (
      "alignment_fasta=s" => \$options{"alignment_fasta"},
      "output_file=s" => \$options{"output_file"}, 
      "signature_max_length=i" => \$options{"signature_max_length"},
      "total_signature_length=i" => \$options{"total_signature_length"},

      "outer_primer_target_length=i" => \$options{"outer_primer_target_length"},
      "outer_primer_min_length=i" => \$options{"outer_primer_min_length"},
      "outer_primer_max_length=i" => \$options{"outer_primer_max_length"},
      "outer_primer_target_tm=f" => \$options{"outer_primer_target_tm"},
      "outer_primer_min_tm=f" => \$options{"outer_primer_min_tm"},
      "outer_primer_max_tm=f" => \$options{"outer_primer_max_tm"},

      "stem_primer_target_length=i" => \$options{"stem_primer_target_length"},
      "stem_primer_min_length=i" => \$options{"stem_primer_min_length"},
      "stem_primer_max_length=i" => \$options{"stem_primer_max_length"},
      "stem_primer_target_tm=f" => \$options{"stem_primer_target_tm"},
      "stem_primer_min_tm=f" => \$options{"stem_primer_min_tm"},
      "stem_primer_max_tm=f" => \$options{"stem_primer_max_tm"},

      "middle_primer_target_length=i" => \$options{"middle_primer_target_length"},
      "middle_primer_min_length=i" => \$options{"middle_primer_min_length"},
      "middle_primer_max_length=i" => \$options{"middle_primer_max_length"},
      "middle_primer_target_tm=f" => \$options{"middle_primer_target_tm"},
      "middle_primer_min_tm=f" => \$options{"middle_primer_min_tm"},
      "middle_primer_max_tm=f" => \$options{"middle_primer_max_tm"},

      "inner_primer_target_length=i" => \$options{"inner_primer_target_length"},
      "inner_primer_min_length=i" => \$options{"inner_primer_min_length"},
      "inner_primer_max_length=i" => \$options{"inner_primer_max_length"},
      "inner_primer_target_tm=f" => \$options{"inner_primer_target_tm"},
      "inner_primer_min_tm=f" => \$options{"inner_primer_min_tm"},
      "inner_primer_max_tm=f" => \$options{"inner_primer_max_tm"},
    
      "max_poly_bases=i" => \$options{"max_poly_bases"}, 
      
      "max_total_degenerate_bases=i" => \$options{"max_total_degenerate_bases"},
      "max_consecutive_degenerate_bases=i" => \$options{"max_consecutive_degenerate_bases"},
      "max_3prime_degenerate_bases=i" => \$options{"max_3prime_degenerate_bases"},
      "max_tolerated_mismatches=i" => \$options{"max_tolerated_mismatches"},
      "three_prime_zone_size=i" => \$options{"three_prime_zone_size"},
      "min_base_frequency=f" => \$options{"min_base_frequency"},
      "entropy_threshold=f" => \$options{"entropy_threshold"},

      "outer_pair_target_length=i" => \$options{"outer_pair_target_length"},
      "middle_pair_target_length=i" => \$options{"middle_pair_target_length"},
      "inner_pair_target_length=i" => \$options{"inner_pair_target_length"},

      "include_stem_primers=i" => \$options{"include_stem_primers"},
      "min_signatures_for_success=i" => \$options{"min_signatures_for_success"},
      "min_primer_spacing=i" => \$options{"min_primer_spacing"},
      "min_inner_pair_spacing=i" => \$options{"min_inner_pair_spacing"},
      "max_overlap_percent=f" => \$options{"max_overlap_percent"},
      "resolve_overlap_by=s" => \$options{"resolve_overlap_by"},
      # --- REDUCTION SPATIALE PAR FENETRE / SPATIAL WINDOW REDUCTION ---
      "window_size=i"    => \$options{"window_size"},    # largeur fenetre en nt (0=desactive)
      "max_per_window=i" => \$options{"max_per_window"}, # max candidats par fenetre

      # --- NOUVEAUX PARAMÈTRES D'ARCHITECTURE ---
      "max_dist_outer_middle=i" => \$options{"max_dist_outer_middle"},
      "max_dist_middle_inner=i" => \$options{"max_dist_middle_inner"},
      # -----------------------------------------

      "primer3_executable=s" => \$options{"primer3_executable"},
      "thermodynamic_path=s" => \$options{"thermodynamic_path"},
      "alignment_format=s" => \$options{"alignment_format"},
      # Ces paramètres sont nécessaires pour Primer3Conserved / These parameters are required for Primer3Conserved
      "dntp_conc=f" => \$options{"dntp_conc"},
      "salt_divalent=f" => \$options{"salt_divalent"},
      "salt_monovalent=f" => \$options{"salt_monovalent"},
      "dna_conc=f" => \$options{"dna_conc"}, # new
      "max_tm_diff=f" => \$options{"max_tm_diff"},
      "max_primer_gen=f" => \$options{"max_primer_gen"}, # new

      # Sigmoid Penalty Parameters
      "penalty_plateau=f" => \$options{"penalty_plateau"},
      "penalty_slope=f" => \$options{"penalty_slope"},

      # --- NOUVEAUX PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES ---
      "primer_min_match_percent=f" => \$options{"primer_min_match_percent"},
      "primer_iupac_min_percent=f" => \$options{"primer_iupac_min_percent"},
      "min_primer_coverage=f" => \$options{"min_primer_coverage"},
      "signature_common_target_min_percent=f" => \$options{"signature_common_target_min_percent"},
      # --------------------------------------------------------

      # TODO: Not sure if the pair target lengths should be exposed to the 
      # user, or adjusted based on other parameters
      #"outer_pair_target_length=i" => \$options{"outer_pair_target_length"}, 
      #"middle_pair_target_length=i" => \$options{"middle_pair_target_length"}, 
      #"inner_pair_target_length=i" => \$options{"inner_pair_target_length"}, 

      "option_file|options_file=s" => \$options{"option_file"},
    );

  my %optionDefaults =
    (
      "signature_max_length" => 400,
      "outer_primer_target_length" => 20,
      "outer_primer_min_length" => 18,
      "outer_primer_max_length" => 23,
      "outer_primer_target_tm" => "60.0",
      "stem_primer_target_length" => 20,
      "stem_primer_min_length" => 18,
      "stem_primer_max_length" => 23,
      "stem_primer_target_tm" => "60.0",
      "middle_primer_target_length" => 20,
      "middle_primer_min_length" => 18,
      "middle_primer_max_length" => 23,
      "middle_primer_target_tm" => "60.0",
      "inner_primer_target_length" => 23,
      "inner_primer_min_length" => 20,
      "inner_primer_max_length" => 26,
      "inner_primer_target_tm" => "62.0",
      "max_poly_bases" => 2,
      "include_stem_primers" => 1,
      "min_signatures_for_success" => 1, # Should probably never go lower / Ne devrait probablement jamais descendre plus bas
      "min_primer_spacing" => 1,
      "min_inner_pair_spacing" => 1,
      # --- NOUVEAUX PARAMÈTRES D'ARCHITECTURE (valeurs par défaut) / NEW ARCHITECTURE PARAMETERS (default values) ---
      "max_dist_outer_middle" => 30,
      "max_dist_middle_inner" => 30,
      # --- PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES ---
      "primer_min_match_percent" => 80,
      "primer_iupac_min_percent" => 98,
      "min_primer_coverage" => 80,
      "signature_common_target_min_percent" => 70,
      # -----------------------------------------
      # Some LAMP-specific approximate targets for a "minimum sized" signature
      # Currently, no penalty is assessed for lengths under the target size, so
      # these sizes are a little larger than they need to be.
      "outer_pair_target_length" => 250, 
      "middle_pair_target_length" => 180,
      "inner_pair_target_length" => 50, 
      "max_overlap_percent" => 0,
      "resolve_overlap_by" => "penalty",
      "dntp_conc" => 1.4,
      "salt_divalent" => 8,
      "salt_monovalent" => 50,
      "dna_conc" => 400,
      "dna_conc" => 400,
      "penalty_plateau" => 0.25,
      "penalty_slope" => 0.15,
      "max_primer_gen" => 10001, # primer3 rounding error off by 1?
      "primer3_executable" => "/usr/bin/primer3_core",
      "thermodynamic_path" => "/etc/primer3_config/",
      "alignment_format"   => "fasta",
      # Reduction spatiale : 0 = desactive (comportement par defaut) / Spatial reduction: 0 = disabled (default behavior)
      "window_size"        => 0,
      "max_per_window"     => 0,
    );

  my $usageString = "Usage:\n" .
    "./lava.pl \n" .
      "    --alignment_fasta <fasta_file>\n" .
      "    --output_file <output_file>\n" .
      "    [--signature_max_length <length, default=" .
        $optionDefaults{"signature_max_length"} .
	">]\n" .
      # Outer primer options
      "    [--outer_primer_target_length <length, default=" .
        $optionDefaults{"outer_primer_target_length"} .
	">]\n" .
      "    [--outer_primer_min_length <length, default=" .
        $optionDefaults{"outer_primer_min_length"} .
	">]\n" .
      "    [--outer_primer_max_length <length, default=" .
        $optionDefaults{"outer_primer_max_length"} .
	">]\n" .
      "    [--outer_primer_target_tm <tm, default=" .
        $optionDefaults{"outer_primer_target_tm"} .
        "C>]\n" .
      "    [--outer_primer_min_tm <tm, default=outer_primer_target_tm - 1.0>]\n" .
      "    [--outer_primer_max_tm <tm, default=outer_primer_target_tm + 1.0>]\n" .
      # STEM primer options
      "    [--stem_primer_target_length <length, default=" . 
        $optionDefaults{"stem_primer_target_length"} .
        ">]\n" .
      "    [--stem_primer_min_length <length, default=" .
        $optionDefaults{"stem_primer_min_length"} .
	">]\n" .
      "    [--stem_primer_max_length <length, default=" .
        $optionDefaults{"stem_primer_max_length"} .
	">]\n" .
      "    [--stem_primer_target_tm <tm, default=" .
        $optionDefaults{"stem_primer_target_tm"} .
	"C>]\n" .
      "    [--stem_primer_min_tm <tm, default=stem_primer_target_tm - 1.0>]\n" .
      "    [--stem_primer_max_tm <tm, default=stem_primer_target_tm + 1.0>]\n" .
      # Middle primer options
      "    [--middle_primer_target_length <length, default=" .
        $optionDefaults{"middle_primer_target_length"} .
	">]\n" .
      "    [--middle_primer_min_length <length, default=" .
        $optionDefaults{"middle_primer_min_length"} .
	">]\n" .
      "    [--middle_primer_max_length <length, default=" .
        $optionDefaults{"middle_primer_max_length"} .
	">]\n" .
      "    [--middle_primer_target_tm <tm, default=" .
        $optionDefaults{"middle_primer_target_tm"} .
	"C>]\n" .
      "    [--middle_primer_min_tm <tm, default=middle_primer_target_tm - 1.0>]\n" .
      "    [--middle_primer_max_tm <tm, default=middle_primer_target_tm + 1.0>]\n" .
      # Outer primer options
      "    [--inner_primer_target_length <length, default=" .
        $optionDefaults{"inner_primer_target_length"} .
	">]\n" .
      "    [--inner_primer_min_length <length, default=" .
        $optionDefaults{"inner_primer_min_length"} .
	">]\n" .
      "    [--inner_primer_max_length <length, default=" .
        $optionDefaults{"inner_primer_max_length"} .
	">]\n" .
      "    [--inner_primer_target_tm <tm, default=" .
        $optionDefaults{"inner_primer_target_tm"} .
	"C>]\n" .
      "    [--inner_primer_min_tm <tm, default=inner_primer_target_tm - 1.0>]\n" .
      "    [--inner_primer_max_tm <tm, default=inner_primer_target_tm + 1.0>]\n" .
      # Other kinds of options
      "    [--max_poly_bases <max, default=" .
        $optionDefaults{"max_poly_bases"} .
	">]\n" .
      "    [--min_primer_spacing <max, default=" .
        $optionDefaults{"min_primer_spacing"} .  
    ">]\n" .
      "    [--min_inner_pair_spacing <max, default=" .
        $optionDefaults{"min_inner_pair_spacing"} .
    ">]\n" .
      "    [--outer_pair_target_length <length, default=" .
        $optionDefaults{"outer_pair_target_length"} .
    ">]\n" .
      "    [--middle_pair_target_length <length, default=" .
      $optionDefaults{"middle_pair_target_length"} .
    ">]\n" .
      "    [--inner_pair_target_length <length, default=" .
      $optionDefaults{"inner_pair_target_length"} .
    ">]\n" .
      "    [--include_stem_primers <length, default=" .
        $optionDefaults{"include_stem_primers"} .
	">]\n" .
      "    [--min_signatures_for_success <length, default=" .
        $optionDefaults{"min_signatures_for_success"} .
  ">]\n" .
    "    [--max_overlap_percent <length, default=" .
      $optionDefaults{"max_overlap_percent"} .
  ">]\n" .
    "    [--dna_conc <length, default=" .
      $optionDefaults{"dna_conc"} .
  ">]\n" .
    "    [--dntp_conc <length, default=" .
      $optionDefaults{"dntp_conc"} .
  ">]\n" .
    "    [--salt_monovalent <length, default=" .
      $optionDefaults{"salt_monovalent"} .
  ">]\n" .
    "    [--salt_divalent <length, default=" .
      $optionDefaults{"salt_divalent"} .
 ">]\n" .
    "    [--max_primer_gen <length, default=" .
      $optionDefaults{"max_primer_gen"} .
	">]\n" .
      # --- PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES ---
      "    [--primer_min_match_percent <percent, default=" .
        $optionDefaults{"primer_min_match_percent"} .
	">]\n" .
      "    [--primer_iupac_min_percent <percent, default=" .
        $optionDefaults{"primer_iupac_min_percent"} .
	">]\n" .
      "    [--min_primer_coverage <percent, default=" .
        $optionDefaults{"min_primer_coverage"} .
	">]\n" .
      "    [--signature_common_target_min_percent <percent, default=" .
        $optionDefaults{"signature_common_target_min_percent"} .
	">]\n" .
      # -----------------------------------------
      "    [--primer3_executable <path_to_primer3, default=" .
        $optionDefaults{"primer3_executable"} .
	">]\n" .
      "    [--thermodynamic_path <path_to_primer3_configuration, default=" .
        $optionDefaults{"thermodynamic_path"} .
	">]\n" .
      "    [--alignment_format <file format of alignment, default=\"" .
        $optionDefaults{"alignment_format"} .
	"\">]\n" .
      "    [--alignment_format <file format of alignment, default=\"" .
        $optionDefaults{"alignment_format"} .
	"\">]\n" .
      "    [--penalty_plateau <float, default=" . $optionDefaults{"penalty_plateau"} . ">]\n" .
      "    [--penalty_slope <float, default=" . $optionDefaults{"penalty_slope"} . ">]\n" .
      "    [--option_file <options_xml> (cmd line options take precedence)]\n";

  # TODO: Probably want to be able to use multiple files for parameter
  # definition, so we can have the thermo parameters set, and separately have
  # the file IO parameters.
  GetOptions(%optionMap);
  loadOptionsFromFile(\%options);
  my $options_r = \%options;

  # TODO: perl-check for file existence cause BioPerl dump isn't useful
  my $alignmentFastaName = optionRequired($options_r, "alignment_fasta", $usageString);
  my $outputFileName = optionRequired($options_r, "output_file", $usageString);

  my $signatureMaxLength = 
    optionWithDefault($options_r, "signature_max_length", 
      $optionDefaults{"signature_max_length"});
  my $totalSignatureLength = 
    optionWithDefault($options_r, "total_signature_length",
      $signatureMaxLength); # Default to max length if not specified

  my $maxTotalDegen = optionWithDefault($options_r, "max_total_degenerate_bases", 2);
  my $maxConsecDegen = optionWithDefault($options_r, "max_consecutive_degenerate_bases", 2);
  my $max3PrimeDegen = optionWithDefault($options_r, "max_3prime_degenerate_bases", 2);
  my $maxToleratedMismatches = optionWithDefault($options_r, "max_tolerated_mismatches", 0);
  my $threePrimeZoneSize = optionWithDefault($options_r, "three_prime_zone_size", 6);
  my $minBaseFrequency = optionWithDefault($options_r, "min_base_frequency", 0.05);
  print "Config: Min Base Frequency = $minBaseFrequency\n";
  my $entropyThreshold = optionWithDefault($options_r, "entropy_threshold", 1.5);

  my $outerPrimerTargetLength =
    optionWithDefault($options_r, "outer_primer_target_length", 
      $optionDefaults{"outer_primer_target_length"});
  my $outerPrimerMinLength =
    optionWithDefault($options_r, "outer_primer_min_length", 
      $optionDefaults{"outer_primer_min_length"});
  if($outerPrimerMinLength > $outerPrimerTargetLength)
  {
    $outerPrimerMinLength = $outerPrimerTargetLength;
  }
  my $outerPrimerMaxLength =
    optionWithDefault($options_r, "outer_primer_max_length", 
      $optionDefaults{"outer_primer_max_length"});
  if($outerPrimerMaxLength < $outerPrimerTargetLength)
  {
    $outerPrimerMaxLength = $outerPrimerTargetLength;
  }

  my $outerPrimerTargetTM =
    optionWithDefault($options_r, "outer_primer_target_tm", 
      $optionDefaults{"outer_primer_target_tm"});
  my $outerPrimerMinTM =
    optionWithDefault($options_r, "outer_primer_min_tm", 
      ($outerPrimerTargetTM - 1.0));
  if($outerPrimerMinTM > $outerPrimerTargetTM)
  {
    $outerPrimerMinTM = $outerPrimerTargetTM;
  }
  my $outerPrimerMaxTM =
    optionWithDefault($options_r, "outer_primer_max_tm", 
      ($outerPrimerTargetTM + 1.0));
  if($outerPrimerMaxTM < $outerPrimerTargetTM)
  {
    $outerPrimerMaxTM = $outerPrimerTargetTM;
  }

  my $stemPrimerTargetLength =
    optionWithDefault($options_r, "stem_primer_target_length", 
      $optionDefaults{"stem_primer_target_length"});
  my $stemPrimerMinLength =
    optionWithDefault($options_r, "stem_primer_min_length", 
      $optionDefaults{"stem_primer_min_length"});
  if($stemPrimerMinLength > $stemPrimerTargetLength)
  {
    $stemPrimerMinLength = $stemPrimerTargetLength;
  }
  my $stemPrimerMaxLength =
    optionWithDefault($options_r, "stem_primer_max_length", 
      $optionDefaults{"stem_primer_max_length"});
  if($stemPrimerMaxLength < $stemPrimerTargetLength)
  {
    $stemPrimerMaxLength = $stemPrimerTargetLength;
  }

  my $stemPrimerTargetTM =
    optionWithDefault($options_r, "stem_primer_target_tm", 
      $optionDefaults{"stem_primer_target_tm"});
  my $stemPrimerMinTM =
    optionWithDefault($options_r, "stem_primer_min_tm", 
      ($stemPrimerTargetTM - 1.0));
  if($stemPrimerMinTM > $stemPrimerTargetTM)
  {
    $stemPrimerMinTM = $stemPrimerTargetTM;
  }
  my $stemPrimerMaxTM =
    optionWithDefault($options_r, "stem_primer_max_tm", 
      ($stemPrimerTargetTM + 1.0));
  if($stemPrimerMaxTM < $stemPrimerTargetTM)
  {
    $stemPrimerMaxTM = $stemPrimerTargetTM;
  }


  my $middlePrimerTargetLength =
    optionWithDefault($options_r, "middle_primer_target_length", 
      $optionDefaults{"middle_primer_target_length"});
  my $middlePrimerMinLength =
    optionWithDefault($options_r, "middle_primer_min_length", 
      $optionDefaults{"middle_primer_min_length"});
  if($middlePrimerMinLength > $middlePrimerTargetLength)
  {
    $middlePrimerMinLength = $middlePrimerTargetLength;
  }
  my $middlePrimerMaxLength =
    optionWithDefault($options_r, "middle_primer_max_length", 
      $optionDefaults{"middle_primer_max_length"});
  if($middlePrimerMaxLength < $middlePrimerTargetLength)
  {
    $middlePrimerMaxLength = $middlePrimerTargetLength;
  }

  my $middlePrimerTargetTM =
    optionWithDefault($options_r, "middle_primer_target_tm", 
      $optionDefaults{"middle_primer_target_tm"});
  my $middlePrimerMinTM =
    optionWithDefault($options_r, "middle_primer_min_tm", 
      ($middlePrimerTargetTM - 1.0));
  if($middlePrimerMinTM > $middlePrimerTargetTM)
  {
    $middlePrimerMinTM = $middlePrimerTargetTM;
  }
  my $middlePrimerMaxTM =
    optionWithDefault($options_r, "middle_primer_max_tm", 
      ($middlePrimerTargetTM + 1.0));
  if($middlePrimerMaxTM < $middlePrimerTargetTM)
  {
    $middlePrimerMaxTM = $middlePrimerTargetTM;
  }

  my $innerPrimerTargetLength =
    optionWithDefault($options_r, "inner_primer_target_length", 
      $optionDefaults{"inner_primer_target_length"});
  my $innerPrimerMinLength =
    optionWithDefault($options_r, "inner_primer_min_length", 
      $optionDefaults{"inner_primer_min_length"});
  if($innerPrimerMinLength > $innerPrimerTargetLength)
  {
    $innerPrimerMinLength = $innerPrimerTargetLength;
  }
  my $innerPrimerMaxLength =
    optionWithDefault($options_r, "inner_primer_max_length", 
      $optionDefaults{"inner_primer_max_length"});
  if($innerPrimerMaxLength < $innerPrimerTargetLength)
  {
    $innerPrimerMaxLength = $innerPrimerTargetLength;
  }

  my $innerPrimerTargetTM =
    optionWithDefault($options_r, "inner_primer_target_tm", 
      $optionDefaults{"inner_primer_target_tm"});
  my $innerPrimerMinTM =
    optionWithDefault($options_r, "inner_primer_min_tm", 
      ($innerPrimerTargetTM - 1.0));
  if($innerPrimerMinTM > $innerPrimerTargetTM)
  {
    $innerPrimerMinTM = $innerPrimerTargetTM;
  }
  my $innerPrimerMaxTM =
    optionWithDefault($options_r, "inner_primer_max_tm", 
      ($innerPrimerTargetTM + 1.0));
  if($innerPrimerMaxTM < $innerPrimerTargetTM)
  {
    $innerPrimerMaxTM = $innerPrimerTargetTM;
  }

  my $maxPolyBases = 
    optionWithDefault($options_r, "max_poly_bases", 
      $optionDefaults{"max_poly_bases"});
  
  my $includeStemPrimers = 
    optionWithDefault($options_r, "include_stem_primers", 
      $optionDefaults{"include_stem_primers"});

  my $maxDeltaTm = 
    optionWithDefault($options_r, "max_tm_diff", 5.0);
  my $signatureCommonTargetMinPercent =
    optionWithDefault($options_r, "signature_common_target_min_percent",
      optionWithDefault($options_r, "min_signatures_for_success",
        $optionDefaults{"signature_common_target_min_percent"}));
  # Lit signature_common_target_min_percent si fourni, sinon se replie sur min_signatures_for_success envoyé par l'IHM Flask
  # Reads signature_common_target_min_percent if provided, otherwise falls back to min_signatures_for_success sent by the Flask GUI
  my $maxSigOverlapPercent = 
    optionWithDefault($options_r, "max_overlap_percent",
      $optionDefaults{"max_overlap_percent"});
  my $resolveOverlapBy = 
    optionWithDefault($options_r, "resolve_overlap_by",
      $optionDefaults{"resolve_overlap_by"});

  my $dnaConc = 
    optionWithDefault($options_r, "dna_conc",
     $optionDefaults{"dna_conc"});
  # Redéclaration de ces paramètres pour éviter les erreurs, même si pas supportés par la version actuelle de Primer3 / Redeclaration of these parameters to prevent errors, even if not supported by the current Primer3 version
  my $dntpConc = optionWithDefault($options_r, "dntp_conc", $optionDefaults{"dntp_conc"});
  my $saltMonovalent = optionWithDefault($options_r, "salt_monovalent", $optionDefaults{"salt_monovalent"});

  my $saltDivalent = optionWithDefault($options_r, "salt_divalent", $optionDefaults{"salt_divalent"});

  my $penaltyPlateau = optionWithDefault($options_r, "penalty_plateau", $optionDefaults{"penalty_plateau"});
  my $penaltySlope = optionWithDefault($options_r, "penalty_slope", $optionDefaults{"penalty_slope"});
  my $maxEnumeratedPrimers = int(
    optionWithDefault($options_r, "max_primer_gen",
    $optionDefaults{"max_primer_gen"}));
    
  my $minPrimerSpacing = 
    optionWithDefault($options_r, "min_primer_spacing", 
      $optionDefaults{"min_primer_spacing"});
  my $minInnerPairSpacing =
    optionWithDefault($options_r, "min_inner_pair_spacing", 
      $optionDefaults{"min_inner_pair_spacing"});
  #print "Max poly: $maxPolyBases\n";

  my $outerPairTargetLength = 
    optionWithDefault($options_r, "outer_pair_target_length", 
      $optionDefaults{"outer_pair_target_length"});
  my $middlePairTargetLength = 
    optionWithDefault($options_r, "middle_pair_target_length",
      $optionDefaults{"middle_pair_target_length"});
  my $innerPairTargetLength =
    optionWithDefault($options_r, "inner_pair_target_length", 
      $optionDefaults{"inner_pair_target_length"});

  # --- CALCUL DYNAMIQUE DES LONGUEURS CIBLES (PipelineUtils) ---
  # --- DYNAMIC TARGET LENGTH CALCULATION (PipelineUtils) ---
  if (exists $options_r->{"max_dist_outer_middle"} || exists $options_r->{"max_dist_middle_inner"})
  {
    my $maxDistOuterMiddle = 
      optionWithDefault($options_r, "max_dist_outer_middle",
        $optionDefaults{"max_dist_outer_middle"});
    my $maxDistMiddleInner =
      optionWithDefault($options_r, "max_dist_middle_inner",
        $optionDefaults{"max_dist_middle_inner"});

    ($middlePairTargetLength, $innerPairTargetLength) = calculateDynamicPairLengths(
      $outerPairTargetLength, $maxDistOuterMiddle, $maxDistMiddleInner, $minInnerPairSpacing
    );
  }
  # --- FIN DU CALCUL DYNAMIQUE ---

  # Eventually want to let the user specify which penalty method
  # is used to calculate the spacing penalty, making the objective function
  # more customizable

  my $primer3ExecutablePath = optionWithDefault($options_r, "primer3_executable",
    $optionDefaults{"primer3_executable"});
  my $thermo_path = optionWithDefault($options_r, "thermodynamic_path",
    $optionDefaults{"thermodynamic_path"});
  my $alignmentFormat = optionWithDefault($options_r, "alignment_format",
    $optionDefaults{"alignment_format"});

  # --- RÉCUPÉRATION DES PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES ---
  my $primerMinMatchPercent = optionWithDefault($options_r, "primer_min_match_percent",
    $optionDefaults{"primer_min_match_percent"});
  my $primerIupacMinPercent = optionWithDefault($options_r, "primer_iupac_min_percent", 
    $optionDefaults{"primer_iupac_min_percent"});
  my $minPrimerCoverage = optionWithDefault($options_r, "min_primer_coverage", 
    $optionDefaults{"min_primer_coverage"});
  # $signatureCommonTargetMinPercent deja declare via min_signatures_for_success (GUI)
  # Already declared above via min_signatures_for_success (GUI parameter)
  
  print "Configuration tolérance mismatches:\n";
  print "  - Match strict minimum: ${primerMinMatchPercent}%\n";
  print "  - Couverture IUPAC minimum: ${primerIupacMinPercent}%\n";
  print "  - Seuil élimination primer: ${minPrimerCoverage}%\n";
  print "  - Intersection commune minimum: ${signatureCommonTargetMinPercent}%\n\n";

  # In theory, the overall score logic belongs in a PrimerSetAnalyzer, 
  # but I hope this helps me optimize the inner loop implementing it
  # here, and only instantiating LAMP signatures for the best combinations
  my $innerPenaltyWeight = "1.2";
  my $stemPenaltyWeight = ".7";
  my $middlePenaltyWeight = "1.1";
  my $outerPenaltyWeight = "1.0";

  my $innerToStemPenaltyWeight = 0.5; # Reduced LAVA 2026
  my $innerToLoopPenaltyWeight = 0.5; 
  my $loopToMiddlePenaltyWeight = 0.5; 
  my $innerToMiddlePenaltyWeight = 0.5;
  my $middleToOuterPenaltyWeight = 0.5;
  my $innerForwardToReversePenaltyWeight = 0.5;

  # Let the games begin...

  # Load the input alignment, could be a single sequence
  # TODO: # Make sure the alignment format option suggestion is working
  my $alignIN = Bio::AlignIO->new(-file => "< $alignmentFastaName", -format => $alignmentFormat);
  my $inputMSA = $alignIN->next_aln();

  if (!$inputMSA || $inputMSA->num_sequences() < 1) {
    print STDERR "ERROR: INPUT_EMPTY - Le fichier ne contient aucune sequence valide.\n";
    exit(2);
  }

  if ($inputMSA->num_sequences() >= 2) {
    if (!$inputMSA->is_flush()) {
      print STDERR "ERROR: INPUT_NOT_ALIGNED - Les sequences n'ont pas toutes la meme longueur. Le fichier doit etre un alignement multiple (MSA), pas des sequences brutes.\n";
      exit(2);
    }

    # Verification supplementaire : Bio::AlignIO (BioPerl) ajoute automatiquement des tirets (-) 
    # a la fin des sequences plus courtes lors de la lecture, ce qui fait que is_flush() renvoie 1 
    # meme sur un fichier FASTA non aligne. Nous verifions donc que les sequences brutes dans le fichier 
    # ont bien la meme longueur avant padding.
    if (open(my $fh_check, '<', $alignmentFastaName)) {
      my %raw_lengths;
      my $cur_len = 0;
      my $cur_id = "";
      while (my $line = <$fh_check>) {
        chomp $line;
        if ($line =~ /^>/) {
          if ($cur_id ne "" && $cur_len > 0) {
            $raw_lengths{$cur_len} = 1;
          }
          $cur_id = $line;
          $cur_len = 0;
        } else {
          $line =~ s/\r//g;
          $cur_len += length($line);
        }
      }
      if ($cur_id ne "" && $cur_len > 0) {
        $raw_lengths{$cur_len} = 1;
      }
      close($fh_check);
      if (scalar(keys %raw_lengths) > 1) {
        print STDERR "ERROR: INPUT_NOT_ALIGNED - Les sequences n'ont pas toutes la meme longueur (" . join(", ", keys %raw_lengths) . " bp). Le fichier doit etre un alignement multiple (MSA), pas des sequences brutes.\n";
        exit(2);
      }
    }
  }

  my $sequenceLength = $inputMSA->length;

  # Extraire les séquences du MSA pour la validation d'intersection commune / Extract MSA sequences for common intersection validation
  my @sequences = ();
  my @sequence_names = ();
  my @sequence_objects = ();
  
  foreach my $sequence ($inputMSA->each_seq()) {
    my $seqContent = $sequence->seq();
    $seqContent = uc($seqContent);  # Convertir en majuscules d'abord / Convert to uppercase first
    $seqContent =~ s/[^ATCG]/N/g;  # Puis normaliser (remplacer caractères non-ADN par N) / Then normalize (replace non-DNA characters with N)
    push @sequences, $seqContent;
    push @sequence_names, $sequence->display_id();
    push @sequence_objects, $sequence;
  }

  # Ideally we would  have separate forward and reverse primer generation,
  # But since Primer3 doesn't accept "PRIMER_INTERNAL_OLIGO_MAX_STABILITY", 
  # we're going to have to filter that out ourselves, but it does mean that we can
  # cheat and just reverse complement the forward primers to get the reverse primers.

  # Enumerate outer primers
  my $outerEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $outerEnumerator->setPrimer3Targets(
    {
      "target_length" => $outerPrimerTargetLength,
      "min_length" => $outerPrimerMinLength,
      "max_length" => $outerPrimerMaxLength,
      "target_tm" => $outerPrimerTargetTM,
      "min_tm" => $outerPrimerMinTM,
      "max_tm" => $outerPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
      "dna_conc" => $dnaConc,
      "dntp_conc" => $dntpConc,
      "salt_monovalent" => $saltMonovalent,
      "salt_divalent" => $saltDivalent,
      "entropy_threshold" => $entropyThreshold,
    });

  print "Enumerating outer forward primers\n";
  my @outerForwardPrimers = getOligosWithMismatchTolerance($outerEnumerator, $inputMSA, 
                                                          $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Outer Forward (F3)");

  print "  Generated \"" .
    scalar(@outerForwardPrimers) .
    "\" outer forward primers (avec tolérance mismatches)\n";


  # Option B : Generation NATIVE des Reverse Outer via Primer3 sur RC(MSA)
  print "Enumerating outer NATIVE reverse primers (Option B)\n";
  my @outerReversePrimers = buildNativeReversePool(
    $outerEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Outer Reverse (B3)"
  );
  print "  Generated \"" . scalar(@outerReversePrimers) . "\" outer native reverse primers\n";


  # Enumerate STEM primers, since the STEM primers extend in the opposite 
  # direction of the other LAMP primers, the back-STEM primers are 
  # generated on the as-is sequence, and the forward-STEM primers are 
  # built in the opposite orientation
  my $stemEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $stemEnumerator->setPrimer3Targets(
    {
      "target_length" => $stemPrimerTargetLength,
      "min_length" => $stemPrimerMinLength,
      "max_length" => $stemPrimerMaxLength,
      "target_tm" => $stemPrimerTargetTM,
      "min_tm" => $stemPrimerMinTM,
      "max_tm" => $stemPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
      "dna_conc" => $dnaConc,
      "dntp_conc" => $dntpConc,
      "salt_monovalent" => $saltMonovalent,
      "salt_divalent" => $saltDivalent,
      "entropy_threshold" => $entropyThreshold,
    });

  # This difference in naming is intentional for now (stemBackPrimers instead of 
  # stemReversePrimers), to serve as a reminder that
  # STEM primers extend the other direction, and that their locations need to be 
  # with the opposite orientation
  
  my @stemBackPrimers = ();
  my @stemForwardPrimers = ();
  
  if($includeStemPrimers == $TRUE) {
    # BSTEM : genere nativement sur le brin + (Back Stem = sens du brin +)
    # BSTEM: natively generated on plus strand (Back Stem = sense of plus strand)
    print "Enumerating STEM BACK (BSTEM) primers on plus strand\n";
    @stemBackPrimers = getOligosWithMismatchTolerance($stemEnumerator, $inputMSA,
                                                        $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                        $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Stem Back (BSTEM)");

    print "  Generated \"" .
      scalar(@stemBackPrimers) .
      "\" STEM BACK (BSTEM) primers\n";

    # FSTEM : Option B - genere nativement sur RC(MSA) pour garantir la protection 3'
    # FSTEM: Option B - natively generated on RC(MSA) to guarantee 3-prime protection
    print "Enumerating STEM FORWARD (FSTEM) NATIVE reverse primers (Option B)\n";
    @stemForwardPrimers = buildNativeReversePool(
      $stemEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Stem Forward (FSTEM)"
    );
    print "  Generated \"" . scalar(@stemForwardPrimers) . "\" STEM FORWARD (FSTEM) native primers\n";
  } else {
    print "STEM primers désactivés - génération ignorée\n";
  }

  # Enumerate middle primers
  my $middleEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $middleEnumerator->setPrimer3Targets(
    {
      "target_length" => $middlePrimerTargetLength,
      "min_length" => $middlePrimerMinLength,
      "max_length" => $middlePrimerMaxLength,
      "target_tm" => $middlePrimerTargetTM,
      "min_tm" => $middlePrimerMinTM,
      "max_tm" => $middlePrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
      "dna_conc" => $dnaConc,
      "dntp_conc" => $dntpConc,
      "salt_monovalent" => $saltMonovalent,
      "salt_divalent" => $saltDivalent,
      "entropy_threshold" => $entropyThreshold,
    });

  print "Enumerating middle forward primers\n";
  my @middleForwardPrimers = getOligosWithMismatchTolerance($middleEnumerator, $inputMSA,
                                                           $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                           $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Middle Forward (F2)");

  print "  Generated \"" .
    scalar(@middleForwardPrimers) .
    "\" middle primers (avec tolérance mismatches)\n";

  # Option B : Generation NATIVE des Reverse Middle via Primer3 sur RC(MSA)
  print "Enumerating middle NATIVE reverse primers (Option B)\n";
  my @middleReversePrimers = buildNativeReversePool(
    $middleEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Middle Reverse (B2)"
  );
  print "  Generated \"" . scalar(@middleReversePrimers) . "\" middle native reverse primers\n";

  # Enumerate inner primers 
  my $innerEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $innerEnumerator->setPrimer3Targets(
    {
      "target_length" => $innerPrimerTargetLength,
      "min_length" => $innerPrimerMinLength,
      "max_length" => $innerPrimerMaxLength,
      "target_tm" => $innerPrimerTargetTM,
      "min_tm" => $innerPrimerMinTM,
      "max_tm" => $innerPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
      "dna_conc" => $dnaConc,
      "dntp_conc" => $dntpConc,
      "salt_monovalent" => $saltMonovalent,
      "salt_divalent" => $saltDivalent,
      "entropy_threshold" => $entropyThreshold,
    });

  print "Enumerating inner forward primers\n";
  my @innerForwardPrimers = getOligosWithMismatchTolerance($innerEnumerator, $inputMSA,
                                                          $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Inner Forward (F1c)");

  print "  Generated \"" .
    scalar(@innerForwardPrimers) .
    "\" inner primers (avec tolérance mismatches)\n";

  # Option B : Generation NATIVE des Reverse Inner via Primer3 sur RC(MSA)
  print "Enumerating inner NATIVE reverse primers (Option B)\n";
  my @innerReversePrimers = buildNativeReversePool(
    $innerEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Inner Reverse (B1c)"
  );
  print "  Generated \"" . scalar(@innerReversePrimers) . "\" inner native reverse primers\n";

  # TODO: want to flip any primer locations to reflect the standard
  # positive strand 5' location notation if they were generated
  # on an anti-sense strand, so all the length-based calculations
  # are handled only here, and the locations are standardized for the
  # rest of the process.

  # Analyze every oligo to get oligo penalty scores
  # Currently sharing one default analyzer for all the primers
  my $outerPrimerAnalyzer = LLNL::LAVA::PrimerAnalyzer::PCRPrimer->new();
  my $middlePrimerAnalyzer = $outerPrimerAnalyzer;
  my $innerPrimerAnalyzer = $outerPrimerAnalyzer;
  my $stemPrimerAnalyzer = $outerPrimerAnalyzer;

  print "Analyzing outer forward primers\n";
  my $outerForwardPrimerMeasurements_r =
    analyzeAll(\@outerForwardPrimers, $outerPrimerAnalyzer);
  print "Analyzing outer reverse primers\n";
  my $outerReversePrimerMeasurements_r =
    analyzeAll(\@outerReversePrimers, $outerPrimerAnalyzer);

  my $stemForwardPrimerMeasurements_r = [];
  my $stemBackPrimerMeasurements_r = [];
  
  if($includeStemPrimers == $TRUE) {
    print "Analyzing STEM \"forward\" primers\n";
    $stemForwardPrimerMeasurements_r =
      analyzeAll(\@stemForwardPrimers, $stemPrimerAnalyzer);
    print "Analyzing STEM \"back\" primers\n";
    $stemBackPrimerMeasurements_r =
      analyzeAll(\@stemBackPrimers, $stemPrimerAnalyzer);
  } else {
    print "Analyse de / Analysis ofs STEM primers ignorée\n";
  }

  print "Analyzing middle forward primers\n";
  my $middleForwardPrimerMeasurements_r = 
    analyzeAll(\@middleForwardPrimers, $middlePrimerAnalyzer);
  print "Analyzing middle reverse primers\n";
  my $middleReversePrimerMeasurements_r = 
    analyzeAll(\@middleReversePrimers, $middlePrimerAnalyzer);

  print "Analyzing inner forward primers\n";
  my $innerForwardPrimerMeasurements_r = 
    analyzeAll(\@innerForwardPrimers, $innerPrimerAnalyzer);
  print "Analyzing inner reverse primers\n";
  my $innerReversePrimerMeasurements_r = 
    analyzeAll(\@innerReversePrimers, $innerPrimerAnalyzer);

  # Sort all primers by 5' start location, and separately by score
  # It's tempting to rely on their current order, but I want to make
  # double sure we get increasing penalty sorting, so I'll do it explicitly.
  print "Sorting primer sets\n";

  # Not using an identifier to cross-reference between the sets, because
  # each location+length pair should be unique
  
  # Outer primers sorted 2 ways
  my @outerForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$outerReversePrimerMeasurements_r};

  my @outerForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$outerReversePrimerMeasurements_r};

  # STEM primers sorted 2 ways (seulement si activés / only if enabled)
  my @stemForwardInfoByLocation = ();
  my @stemBackInfoByLocation = ();
  my @stemForwardInfoByPenalty = ();
  my @stemBackInfoByPenalty = ();
  
  if($includeStemPrimers) {
    @stemForwardInfoByLocation =
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getLocation()] } 
      @{$stemForwardPrimerMeasurements_r};
    @stemBackInfoByLocation =
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getLocation()] } 
      @{$stemBackPrimerMeasurements_r};

    @stemForwardInfoByPenalty =
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getPenalty()] } 
      @{$stemForwardPrimerMeasurements_r};
    @stemBackInfoByPenalty =
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getPenalty()] } 
      @{$stemBackPrimerMeasurements_r};
  }

  # Middle primers sorted 2 ways
  my @middleForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$middleReversePrimerMeasurements_r};

  my @middleForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$middleReversePrimerMeasurements_r};

  # Inner primers sorted 2 ways
  my @innerForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$innerReversePrimerMeasurements_r};

  my @innerForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$innerReversePrimerMeasurements_r};

  print "Enumerating signatures\n";

  #=============================================================================
  # THE BIG MERGE (Single-Pass Architecture - LAVA 2026)
  #=============================================================================
  # Remplace l'ancien système multi-passes (combinationPlan x12 itérations).
  # Toutes les listes maîtresses sont construites UNE SEULE FOIS avec
  # max_overlap_percent = 100, préservant toute la diversité combinatoire.
  # La réduction par chevauchement n'est appliquée qu'en fin de pipeline,
  # sur les signatures complètes (via $maxSigOverlapPercent).
  # Replaces the old multi-pass system (combinationPlan x12 iterations).
  # All master lists are built ONCE with max_overlap_percent = 100,
  # preserving full combinatorial diversity for the nested loops.
  #=============================================================================

  my $allFoundSignatures_r = []; # Collect ALL signatures before any reduction

  print "Building Master Primer Lists (The Big Merge)...\n";

  my $bigMerge = buildBigMerge({
    inner_f_loc     => \@innerForwardInfoByLocation,
    inner_f_pen     => \@innerForwardInfoByPenalty,
    inner_r_loc     => \@innerReverseInfoByLocation,
    inner_r_pen     => \@innerReverseInfoByPenalty,
    special_f_loc   => \@stemForwardInfoByLocation,
    special_f_pen   => \@stemForwardInfoByPenalty,
    special_r_loc   => \@stemBackInfoByLocation,
    special_r_pen   => \@stemBackInfoByPenalty,
    include_special => $includeStemPrimers,
    middle_f_loc    => \@middleForwardInfoByLocation,
    middle_f_pen    => \@middleForwardInfoByPenalty,
    middle_r_loc    => \@middleReverseInfoByLocation,
    middle_r_pen    => \@middleReverseInfoByPenalty,
    outer_f_loc     => \@outerForwardInfoByLocation,
    outer_f_pen     => \@outerForwardInfoByPenalty,
    outer_r_loc     => \@outerReverseInfoByLocation,
    outer_r_pen     => \@outerReverseInfoByPenalty,
    # Reduction spatiale : garder les K meilleurs par fenetre de W nt / Spatial reduction: keep K best per W-nt window
    window_size     => optionWithDefault($options_r, "window_size",    $optionDefaults{"window_size"}),
    max_per_window  => optionWithDefault($options_r, "max_per_window", $optionDefaults{"max_per_window"}),
  });

  my $masterInnerF_r       = $bigMerge->{inner_f};
  my $masterInnerF_data_r  = $bigMerge->{inner_f_data};
  my $masterInnerR_r       = $bigMerge->{inner_r};
  my $masterInnerR_data_r  = $bigMerge->{inner_r_data};
  my $masterStemF_r        = $bigMerge->{special_f};
  my $masterStemF_data_r   = $bigMerge->{special_f_data};
  my $masterStemR_r        = $bigMerge->{special_r};
  my $masterStemR_data_r   = $bigMerge->{special_r_data};
  my $masterMiddleF_r      = $bigMerge->{middle_f};
  my $masterMiddleF_data_r = $bigMerge->{middle_f_data};
  my $masterMiddleR_r      = $bigMerge->{middle_r};
  my $masterMiddleR_data_r = $bigMerge->{middle_r_data};
  my $masterOuterF_r       = $bigMerge->{outer_f};
  my $masterOuterF_data_r  = $bigMerge->{outer_f_data};
  my $masterOuterR_r       = $bigMerge->{outer_r};
  my $masterOuterR_data_r  = $bigMerge->{outer_r_data};

  # Pointeurs vers les listes maîtresses pour la boucle de combinaison
  # Pointers to master lists for the combination loop
  my $innerForwardSubset_r     = $masterInnerF_r;
  my $innerForwardSubsetData_r = $masterInnerF_data_r;
  my $innerReverseSubset_r     = $masterInnerR_r;
  my $innerReverseSubsetData_r = $masterInnerR_data_r;
  my $stemForwardSubset_r      = $masterStemF_r;
  my $stemForwardSubsetData_r  = $masterStemF_data_r;
  my $stemReverseSubset_r      = $masterStemR_r;
  my $stemReverseSubsetData_r  = $masterStemR_data_r;
  my $middleForwardSubset_r    = $masterMiddleF_r;
  my $middleForwardSubsetData_r= $masterMiddleF_data_r;
  my $middleReverseSubset_r    = $masterMiddleR_r;
  my $middleReverseSubsetData_r= $masterMiddleR_data_r;
  my $outerForwardSubset_r     = $masterOuterF_r;
  my $outerForwardSubsetData_r = $masterOuterF_data_r;
  my $outerReverseSubset_r     = $masterOuterR_r;
  my $outerReverseSubsetData_r = $masterOuterR_data_r;

  if ($includeStemPrimers == $TRUE) {
    printf "Master Lists Counts:\n  Inner F: %d, R: %d\n  STEM F: %d, R: %d\n  Middle F: %d, R: %d\n  Outer F: %d, R: %d\n",
      scalar(@{$masterInnerF_r}), scalar(@{$masterInnerR_r}),
      scalar(@{$masterStemF_r}),  scalar(@{$masterStemR_r}),
      scalar(@{$masterMiddleF_r}), scalar(@{$masterMiddleR_r}),
      scalar(@{$masterOuterF_r}), scalar(@{$masterOuterR_r});
  } else {
    printf "Master Lists Counts:\n  Inner F: %d, R: %d\n  Middle F: %d, R: %d\n  Outer F: %d, R: %d\n",
      scalar(@{$masterInnerF_r}), scalar(@{$masterInnerR_r}),
      scalar(@{$masterMiddleF_r}), scalar(@{$masterMiddleR_r}),
      scalar(@{$masterOuterF_r}), scalar(@{$masterOuterR_r});
  }

################################################################################
# Now that the subgroups are picked, try to create signatures for this set of 
# combination of subgroups
################################################################################

    # Vérification de l'état des primers STEM / Checking STEM primers state
    
    if($includeStemPrimers == $TRUE) {
      print "Primer counts used for this plan iteration (WITH STEMS):\n" .
        "  " . scalar(@{$innerForwardSubset_r}) . " Inner Forward\n" .
        "  " . scalar(@{$innerReverseSubset_r}) . " Inner Reverse\n" .
        "  " . scalar(@{$stemForwardSubset_r}) . " STEM Forward\n" .
        "  " . scalar(@{$stemReverseSubset_r}) . " STEM Reverse\n" .
        "  " . scalar(@{$middleForwardSubset_r}) . " Middle Forward\n" .
        "  " . scalar(@{$middleReverseSubset_r}) . " Middle Reverse\n" .
        "  " . scalar(@{$outerForwardSubset_r}) . " Outer Forward\n" .
        "  " . scalar(@{$outerReverseSubset_r}) . " Outer Reverse\n\n";
    } else {
      print "Primer counts used for this plan iteration (LAMP 6-primer mode - NO STEMS):\n" .
        "  " . scalar(@{$innerForwardSubset_r}) . " Inner Forward\n" .
        "  " . scalar(@{$innerReverseSubset_r}) . " Inner Reverse\n" .
        "  0 STEM Forward (disabled)\n" .
        "  0 STEM Reverse (disabled)\n" .
        "  " . scalar(@{$middleForwardSubset_r}) . " Middle Forward\n" .
        "  " . scalar(@{$middleReverseSubset_r}) . " Middle Reverse\n" .
        "  " . scalar(@{$outerForwardSubset_r}) . " Outer Forward\n" .
        "  " . scalar(@{$outerReverseSubset_r}) . " Outer Reverse\n\n";
    }

    
    my $innerForwardCount = scalar(@{$innerForwardSubset_r});
    my $innerReverseCount = scalar(@{$innerReverseSubset_r});
    my $stemForwardCount = scalar(@{$stemForwardSubset_r});
    my $stemReverseCount = scalar(@{$stemReverseSubset_r});
    my $middleForwardCount = scalar(@{$middleForwardSubset_r});
    my $middleReverseCount = scalar(@{$middleReverseSubset_r});
    my $outerForwardCount = scalar(@{$outerForwardSubset_r});
    my $outerReverseCount = scalar(@{$outerReverseSubset_r});


    # TODO: HERE! (should)
    # Pre-compute top Middle->Outer pairings, so we don't keep iterating over them
    #
    # Can't pre-compute others though, because STEM primer may or may not be needed?
    # (Maybe conditionally compute other sets...?)

    #my $innerPairCount = scalar(@{$innerPairs_r});
    my @bestSignatureForInnerForward = ();
    #my $bestSignatureCount = 0;

    # Pre-compute a set of distance penalties for faster use
    # ------------------------------------------------------
    # GENERATION PROPORTIONNELLE SIGMOÏDE (LAVA 2026)
    my $geometry = calculate_proportional_geometry($totalSignatureLength);
    
    # Générer des pénalités spécifiques pour chaque distance / Generate specific penalties for each distance
    my $f2_f1_target = $geometry->{'f2_f1_target'};
    my $loop_target = int($f2_f1_target / 2);
    
    my $innerToLoopPenalties_r = generateDistancePenalties($signatureMaxLength, $loop_target, $penaltyPlateau, $penaltySlope);
    my $loopToMiddlePenalties_r = generateDistancePenalties($signatureMaxLength, $loop_target, $penaltyPlateau, $penaltySlope);
    my $middleToOuterPenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'f3_f2_target'}, $penaltyPlateau, $penaltySlope);
    my $innerToMiddlePenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'f2_f1_target'}, $penaltyPlateau, $penaltySlope);
    my $innerToInnerPenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'inner_target'}, $penaltyPlateau, $penaltySlope);

    print "Generating Sigmoid Penalties (Core.pm)...\n";

    # To remember the optimum combination with
    # 3 columns: STEM, middle, outer
    my @bestForwardInfos = (); 
    # 2 columns: spacing_penalty, primer3_penalty
    my @bestForwardPenalties = ();

    # To help short-cut when no possibilities are found
    my $forwardSetCount = 0;

    print "Scanning Forward Primer Combinations...\n";
    # Barre de progression Flask pour les signatures Forward / Flask progress bar for Forward signatures
    my $_sig_fwd_t0   = time();
    my $_sig_fwd_done = 0;
    my $_sig_fwd_hits = 0;  # Nombre de signatures Forward trouvees / Forward signatures found
    print STDERR "  Recherche combinatoire Stem Forward: $innerForwardCount amorces F1c...\n";

    for(my $innerIndex = 0; $innerIndex < $innerForwardCount; $innerIndex++)
    {
      # Emission LAVA-PROGRESS toutes les 50 iterations / Emit LAVA-PROGRESS every 50 iterations
      $_sig_fwd_done = $innerIndex + 1;
      if ($_sig_fwd_done % 50 == 0 || $_sig_fwd_done == $innerForwardCount) {
        if ($_LAVA_IS_TTY || 1) {
          my $elapsed = time() - $_sig_fwd_t0 + 0.001;
          my $eta = ($_sig_fwd_done < $innerForwardCount)
                    ? int(($innerForwardCount - $_sig_fwd_done) / ($_sig_fwd_done / $elapsed))
                    : 0;
          my $rate = $_sig_fwd_done / $elapsed;
          printf("[LAVA-PROGRESS] Signatures Stem Forward|%d|%d|Sig: %d|%.1f it/s|%d\n",
                 $_sig_fwd_done, $innerForwardCount, $_sig_fwd_hits, $rate, $eta);
        }
      }
      my $innerInfo = $innerForwardSubset_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty, $innerTm) = 
        @{$innerForwardSubsetData_r->[$innerIndex]};

      my $bestSetPenalty = 1000000; # Riduculously large starting value

      # Calculate the first and last base locations to consider for the 
      # STEM forward primer (on the plus strand darn it... need an inversion)
      my $searchStartAt = $innerLocation - $signatureMaxLength +
        $innerLength + 20; # 20 represents 2 other primer min lengths...?
      if($searchStartAt < 0)
      {
	$searchStartAt = 0;
      }

      # STEM ARCHITECTURE: STEM primers are placed AFTER inner forward primer (F1c)
      # et AVANT le B1c. La borne max = F1c_end + innerPairTargetLength/2
      # STEM Start At and End At are indexes where the STEM SEARCH starts and ends
      # FSTEM must lie between F1c (innerLocation+innerLength) and the midpoint F1-B1
      my $stemStartAt = $innerLocation + $innerLength + $minPrimerSpacing;
      my $stemEndAt   = $innerLocation + $innerLength + int($innerPairTargetLength / 2);
      if($stemEndAt < $stemStartAt) { $stemEndAt = $stemStartAt + 50; }  # fallback
      if($stemEndAt < 0) { $stemEndAt = 0; }

      # Progression : amorce inner suivante (debug retiré / debug removed)
      #print "(PA $sequenceLength long, start at $searchStartAt, maxLen $signatureMaxLength, $innerLocation->$stemStartAt-$stemEndAt)";
      #print "(PA* $stemStartAt-$stemEndAt)";

      # If no STEM primers sought, then overwrite the STEM primer list with the
      # single placeholder, one-length (but ideally zero-length), zero-penalty 
      # STEM primer, placed at the end of the inner primer, to make sure it 
      # appears to fit within the acceptable locations
      if($includeStemPrimers == $FALSE)
      {
	#print "\n\n\nNO STEM?!\n\n\n";
        my $placeHolderPrimer = LLNL::LAVA::Oligo->new(
	  {
            "sequence" => "N",
	    "location" => $stemEndAt + 1, # Extra position to un-do length of 1
	    "strand" => "minus",
	  });
        $placeHolderPrimer->setTag("primer3_penalty", 0);
        $placeHolderPrimer->setTag("primer3_tm", 0);

        my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new(
	  {
            "penalty" => 0,
	    "sequence" => $placeHolderPrimer->sequence(),
	    "location" =>$placeHolderPrimer->location(),
	    "length" => $placeHolderPrimer->length(),
	    "analyzed_primer" => $placeHolderPrimer,
	  });

        $stemForwardSubset_r = [$placeHolderInfo];
        $stemForwardSubsetData_r = [[$stemEndAt + 1, 1, 0]]; # [location, length, penalty]
	$stemForwardCount = 1; 
      }

      # Start of the 3-level nested loop for forward primers.
      # Should  exhaustively iterate over STEM, middle, outer 
      # combinations based on the inner pair
      for(my $i = 0; $i < $stemForwardCount; $i++)
      {
	my $stemInfo = $stemForwardSubset_r->[$i];
	my ($stemLocation, $stemLength, $stemPenalty, $stemTm) = 
	  @{$stemForwardSubsetData_r->[$i]};
        #my $stemLocation = $stemInfo->getLocation();
        #my $stemLength = $stemInfo->getLength();

	# Special inversion for STEM primer, because forward STEM primer was
	# designed on the minus strand.
	
        # Seek to the first STEM primer within range
	# but, accept placeholder STEM primer
        if($stemLocation < $stemStartAt &&
	   $stemLength != 1)
	{
	  next;
	}

	# Stop when STEM primer goes out of range	
	if($stemLocation > $stemEndAt &&
	   $stemLength != 1)
	{
	  last;
	}

        # Progression : itération STEM (debug retiré / debug removed)

        # No check for STEM->inner overlap, because stemEndAt is the location limiter

	# STEM ARCHITECTURE: Middle primers are positioned independently of STEM primers
	# since STEM primers are now between inner forward and reverse primers
        my $middleStartAt = $searchStartAt;
        my $middleEndAt = $innerLocation - $minPrimerSpacing;  
	if($middleEndAt < 0)
	{
	  $middleEndAt = 0;
	}
        #print "(LA $middleStartAt-$middleEndAt)";
        # STEM ARCHITECTURE: Calculate distance from inner forward to STEM forward
        my $innerToStemDistance = $stemLocation - ($innerLocation + $innerLength);
        # Ensure distance is non-negative
        if($innerToStemDistance < 0) { $innerToStemDistance = 0; }

        # --- DYNAMIC THERMAL FILTER (Stem vs Inner) ---
        # Only check if stem exists (not placeholder which has 0 tm)
        if ($includeStemPrimers == $TRUE && $stemTm > 0) {
            next if (abs($innerTm - $stemTm) > $maxDeltaTm);
        }

        for(my $j = 0; $j < $middleForwardCount; $j++)
	{
	  my $middleInfo = $middleForwardSubset_r->[$j];
	  my ($middleLocation, $middleLength, $middlePenalty, $midTm) = 
	    @{$middleForwardSubsetData_r->[$j]};

          #my $middleLocation = $middleInfo->getLocation();
          #my $middleLength = $middleInfo->getLength();

	  # Seek to the first middle primer within range
	  if($middleLocation < $middleStartAt)
	  {
	    next;
	  }

	  # Stop when middle primer goes out of range
	  if($middleLocation > $middleEndAt)
	  {
	    last;
	  }

          # Progression : itération Middle (debug retiré / debug removed)
          
          # --- DYNAMIC THERMAL FILTER ---
          # Check Middle vs Stem (if Stem exists) or Middle vs Inner (if no Stem)
          if ($includeStemPrimers == $TRUE && $stemTm > 0) {
              next if (abs($stemTm - $midTm) > $maxDeltaTm); 
          } else {
             # If no stem, Middle follows Inner
             next if (abs($innerTm - $midTm) > $maxDeltaTm);
          }

	  # STEM ARCHITECTURE: Only check spacing between middle and inner primers
	  # STEM primers are no longer between middle and inner
	  if($middleLocation + $middleLength + $minPrimerSpacing > $innerLocation)
	  {
	    next;
	  }

	  my $outerStartAt = $searchStartAt;
          my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;

         #print "(MA $outerStartAt-$outerEndAt)";

         # STEM ARCHITECTURE: Calculate direct distance between inner and middle primers
         # STEM primers are no longer between them
          my $innerToMiddleDistance = $innerLocation - ($middleLocation + $middleLength);
          # Ensure distance is non-negative
          if($innerToMiddleDistance < 0) { $innerToMiddleDistance = 0; }

               

          for(my $k = 0; $k < $outerForwardCount; $k++)
	  {
	    my $outerInfo = $outerForwardSubset_r->[$k];
	    my ($outerLocation, $outerLength, $outerPenalty, $outTm) = 
	      @{$outerForwardSubsetData_r->[$k]};
	    #my $outerLocation = $outerInfo->getLocation();
            #my $outerLength = $outerInfo->getLength();
 
            # Seek to first outer primer within range
	    if($outerLocation < $outerStartAt)
	    {
	      next;
	    }

	    # Stop when outer primer goes out of range
	    if($outerLocation > $outerEndAt)
	    {
              last;
	    }

            #print "O";
            # Next primer if this outer doesn't leave enough spacing to the middle primer
            if($outerLocation + $outerLength + $minPrimerSpacing >
	      $middleLocation)
	    {
	      next;
	    }
            #print "(OA)";
            
            # --- DYNAMIC THERMAL FILTER (Middle vs Outer) ---
            next if (abs($midTm - $outTm) > $maxDeltaTm);

	    # Inter-primer distance used for calculating spacing penalty
            # Calculate middle to outer distance
            my $middleToOuterDistance = $middleLocation - ($outerLocation + $outerLength);
            if($middleToOuterDistance < 0) { $middleToOuterDistance = 0; }
            
            my $spacingPenalty = 0;
            my $primer3Penalty = 0;
            my $detailStr = "";
            # Clamper les distances au max des tableaux de penalites (evite index OOB silencieux → undef)
            # Clamp distances to max penalty array index (prevents silent OOB → undef penalty)
            my $maxPenIdx = $signatureMaxLength - 1;
            my $d_stem   = ($innerToStemDistance   < $maxPenIdx) ? $innerToStemDistance   : $maxPenIdx;
            my $d_middle = ($innerToMiddleDistance < $maxPenIdx) ? $innerToMiddleDistance : $maxPenIdx;
            my $d_outer  = ($middleToOuterDistance < $maxPenIdx) ? $middleToOuterDistance : $maxPenIdx;

            if($includeStemPrimers == $TRUE)
            {
              # STEM ARCHITECTURE: Calculate spacing between inner-STEM and inner-middle
              $spacingPenalty = 
                ($innerToLoopPenalties_r->[$d_stem]    * $innerToStemPenaltyWeight) +
                ($innerToMiddlePenalties_r->[$d_middle] * $innerToMiddlePenaltyWeight) +
                ($middleToOuterPenalties_r->[$d_outer]  * $middleToOuterPenaltyWeight);
              $primer3Penalty = 
                $innerPenalty * $innerPenaltyWeight +
                $stemPenalty  * $stemPenaltyWeight +
                $middlePenalty * $middlePenaltyWeight +
                $outerPenalty * $outerPenaltyWeight;
              
              $detailStr = sprintf("Spc[I_S:%.1f I_M:%.1f M_O:%.1f] Thm[I:%.1f S:%.1f M:%.1f O:%.1f]",
                    ($innerToLoopPenalties_r->[$d_stem]    * $innerToStemPenaltyWeight),
                    ($innerToMiddlePenalties_r->[$d_middle] * $innerToMiddlePenaltyWeight),
                    ($middleToOuterPenalties_r->[$d_outer]  * $middleToOuterPenaltyWeight),
                    ($innerPenalty * $innerPenaltyWeight),
                    ($stemPenalty  * $stemPenaltyWeight),
                    ($middlePenalty * $middlePenaltyWeight),
                    ($outerPenalty * $outerPenaltyWeight));
            }
            else
            {
              $spacingPenalty = 
                ($innerToMiddlePenalties_r->[$d_middle] * $innerToMiddlePenaltyWeight) +
                ($middleToOuterPenalties_r->[$d_outer]  * $middleToOuterPenaltyWeight);
              $primer3Penalty = 
                $innerPenalty * $innerPenaltyWeight +
                $middlePenalty * $middlePenaltyWeight +
                $outerPenalty * $outerPenaltyWeight;
              
              $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]",
                    ($innerToMiddlePenalties_r->[$d_middle] * $innerToMiddlePenaltyWeight),
                    ($middleToOuterPenalties_r->[$d_outer]  * $middleToOuterPenaltyWeight),
                    ($innerPenalty * $innerPenaltyWeight),
                    ($middlePenalty * $middlePenaltyWeight),
                    ($outerPenalty * $outerPenaltyWeight));
            }
           
            my $forwardSetPenalty = $spacingPenalty + $primer3Penalty;
            if($forwardSetPenalty < $bestSetPenalty)
            {
              $bestForwardInfos[$innerIndex] = [$stemInfo, $middleInfo, $outerInfo];
              $bestSetPenalty = $forwardSetPenalty;
              $forwardSetCount++;
              $_sig_fwd_hits++;  # Compteur de signatures Fwd / Fwd signature counter
              $bestForwardPenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty, $detailStr];
            }
          } # End forward outer iteration
        } # End forward middle iteration
      } # End forward STEM iteration
    } # End forward inner iteration

    # Finaliser la barre Forward / Finalize Forward bar
    printf(STDERR "\r%-80s\n", "") if $_LAVA_IS_TTY;
    print "  [Stem Fwd] $forwardSetCount combinaisons Forward trouvees sur $innerForwardCount amorces F1c.\n";

    # Stop trying if no forward primer sets were found
    if($forwardSetCount == 0)
    {
      print "No valid forward primer combinations found.\n";
      next;
    }

    print "Scanning Reverse Primer Combinations...\n";

    my @bestSigantureForInnerReverse = ();

    # To remember the optimum combination with
    # 3 columns: STEM, middle, outer
    my @bestReverseInfos = (); 
    # 2 columns: spacing_penalty, primer3_penalty
    my @bestReversePenalties = ();

    ## To help short-cut when no possibilities are found
    #my $reverseSetCount = 0;

    # Barre de progression Flask pour les signatures Reverse / Flask progress bar for Reverse signatures
    my $_sig_rev_t0   = time();
    my $_sig_rev_done = 0;
    my $_sig_rev_hits = 0;  # Nombre de signatures Reverse trouvees / Reverse signatures found
    print STDERR "  Recherche combinatoire Stem Reverse: $innerReverseCount amorces B1c...\n";

    for(my $innerIndex = 0; $innerIndex < $innerReverseCount; $innerIndex++)
    {
      # Emission LAVA-PROGRESS toutes les 50 iterations / Emit LAVA-PROGRESS every 50 iterations
      $_sig_rev_done = $innerIndex + 1;
      if ($_sig_rev_done % 50 == 0 || $_sig_rev_done == $innerReverseCount) {
        if ($_LAVA_IS_TTY || 1) {
          my $elapsed = time() - $_sig_rev_t0 + 0.001;
          my $eta = ($_sig_rev_done < $innerReverseCount)
                    ? int(($innerReverseCount - $_sig_rev_done) / ($_sig_rev_done / $elapsed))
                    : 0;
          my $rate = $_sig_rev_done / $elapsed;
          printf("[LAVA-PROGRESS] Signatures Stem Reverse|%d|%d|Sig: %d|%.1f it/s|%d\n",
                 $_sig_rev_done, $innerReverseCount, $_sig_rev_hits, $rate, $eta);
        }
      }
      my $innerInfo = $innerReverseSubset_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty, $innerTm) = 
        @{$innerReverseSubsetData_r->[$innerIndex]};

      #my $innerLocation = $innerInfo->getLocation();
      #my $innerLength = $innerInfo->getLength();

      my $bestSetPenalty = 1000000; # Riduculously large starting value

      # Calculate the first and last base locations to consider for the 
      # STEM reverse primer (on the plus strand, so need inversion again?)
      my $searchEndAt = $innerLocation + $signatureMaxLength - 
	$innerLength - 20; # -20 represents 2 other primer min lengths.
	
      # STEM ARCHITECTURE: BSTEM primers are positioned BEFORE inner reverse primer (B1c)
      # et APRES le F1c. La borne min = B1c_start - innerPairTargetLength/2
      # BSTEM must lie between the midpoint F1-B1 and B1c
      my $stemStartAt = $innerLocation - int($innerPairTargetLength / 2);
      my $stemEndAt   = $innerLocation - $minPrimerSpacing;
      if($stemStartAt < 0) { $stemStartAt = 0; }
      if($stemStartAt > $stemEndAt) { $stemStartAt = ($stemEndAt > 50) ? $stemEndAt - 50 : 0; }  # fallback

      # If no STEM primers sought, then overwrite the STEM primer list with the
      # single placeholder, one-length (but ideally zero-length), zero-penalty
      # STEM primer, placed at the end of the inner primer, to make sure it 
      # appears to fit within the acceptable locations
      if($includeStemPrimers == $FALSE)
      {
        my $placeHolderPrimer = LLNL::LAVA::Oligo->new(
	  {
            "sequence" => "N",
	    "location" => $stemStartAt - 1, # Extra position to un-do length of 1
	    "strand" => "plus",
	  });
        $placeHolderPrimer->setTag("primer3_penalty", 0);
        $placeHolderPrimer->setTag("primer3_tm", 0);

        my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new(
	  { 
            "penalty" => 0,
	    "sequence" => $placeHolderPrimer->sequence(),
	    "location" =>$placeHolderPrimer->location(),
	    "length" => $placeHolderPrimer->length(),
	    "analyzed_primer" => $placeHolderPrimer,
	  });

        $stemReverseSubset_r = [$placeHolderInfo];
        $stemReverseSubsetData_r = [[$stemEndAt + 1, 1, 0]]; # [location, length, penalty]
	$stemReverseCount = 1; 
      }

      # Start of the 3-level nested loop for reverse primers.
      # Should  exhaustively iterate over STEM, middle, outer 
      # combinations based on the inner pair
      for(my $i = 0; $i < $stemReverseCount; $i++)
      {
	my $stemInfo = $stemReverseSubset_r->[$i];
	my ($stemLocation, $stemLength, $stemPenalty, $stemTm) = 
	  @{$stemReverseSubsetData_r->[$i]};

        #my $stemLocation = $stemInfo->getLocation();
        #my $stemLength = $stemInfo->getLength();

	# Special inversion for STEM primer, because reverse STEM primer was
	# designed on the minus strand.
	
        # Seek to the first STEM primer within range 
	# but, accept placeholder STEM primer
        if($stemLocation < $stemStartAt &&
	   $stemLength != 1)
	{
	  next;
	}

	# Stop when STEM primer goes out of range	
	if($stemLocation > $stemEndAt &&
	   $stemLength != 1)
	{
	  last;
	}

        #print "L";

        # No check for STEM->inner overlap because stemStartAt is the location limiter

	# STEM ARCHITECTURE: Middle primers are positioned independently of STEM primers  
	# since STEM primers are now between inner forward and reverse primers
        my $middleStartAt = $innerLocation + $minPrimerSpacing;
        my $middleEndAt = $searchEndAt;

        # STEM ARCHITECTURE: Calculate distance from STEM back to inner reverse  
	my $innerToStemDistance = ($innerLocation - $innerLength) - ($stemLocation + $stemLength);
        # Ensure distance is non-negative
        if($innerToStemDistance < 0) { $innerToStemDistance = 0; }

        # --- DYNAMIC THERMAL FILTER (Stem vs Inner) ---
        if ($includeStemPrimers == $TRUE && $stemTm > 0) {
            next if (abs($innerTm - $stemTm) > $maxDeltaTm);
        }

        for(my $j = 0; $j < $middleReverseCount; $j++)
	{
	  my $middleInfo = $middleReverseSubset_r->[$j];
	  my ($middleLocation, $middleLength, $middlePenalty, $midTm) = 
	    @{$middleReverseSubsetData_r->[$j]};

          #my $middleLocation = $middleInfo->getLocation();
          #my $middleLength = $middleInfo->getLength();

	  # Seek to the first middle primer within range
	  if($middleLocation < $middleStartAt)
	  {
	    next;
	  }

	  # Stop when middle primer goes out of range
	  if($middleLocation > $middleEndAt)
	  {
	    last;
	  }

          # Progression : itération Middle reverse (debug retiré / debug removed)
          
          # --- DYNAMIC THERMAL FILTER ---
          if ($includeStemPrimers == $TRUE && $stemTm > 0) {
              next if (abs($stemTm - $midTm) > $maxDeltaTm);
          } else {
              next if (abs($innerTm - $midTm) > $maxDeltaTm);
          }

          # STEM ARCHITECTURE: Only check spacing between middle and inner primers
	  # STEM primers are no longer between middle and inner  
          if($middleLocation - $middleLength - $minPrimerSpacing < $innerLocation)
          {
            next;
	  }

          my $outerStartAt = $middleLocation + 1 + $minPrimerSpacing;
          my $outerEndAt = $searchEndAt;


          for(my $k = 0; $k < $outerReverseCount; $k++)
	  {
	    my $outerInfo = $outerReverseSubset_r->[$k];
	    my ($outerLocation, $outerLength, $outerPenalty, $outTm) = 
	      @{$outerReverseSubsetData_r->[$k]};

	    #my $outerLocation = $outerInfo->getLocation();
            #my $outerLength = $outerInfo->getLength();
 
            # Seek to first outer primer within range
	    if($outerLocation < $outerStartAt)
	    {
	      next;
	    }

	    # Stop when outer primer goes out of range
	    if($outerLocation > $outerEndAt)
	    {
              last;
	    }

            #print "O";
            
            # --- DYNAMIC THERMAL FILTER ---
            next if (abs($midTm - $outTm) > $maxDeltaTm);

            # Next primer if this outer doesn't leave enough spacing to the middle primer
	    if($outerLocation - $outerLength - $minPrimerSpacing <
	       $middleLocation)
	    {
	      next;
	    }

	    # Inter-primer distance used for calculating spacing penalty
            # Calculate distances for reverse primers
            my $middleToOuterDistance = ($outerLocation - $outerLength) - $middleLocation;
            if($middleToOuterDistance < 0) { $middleToOuterDistance = 0; }
            my $innerToMiddleDistance = ($middleLocation - $middleLength) - ($innerLocation + 1);
            if($innerToMiddleDistance < 0) { $innerToMiddleDistance = 0; }
            
            my $spacingPenalty = 0;
            my $primer3Penalty = 0;
            my $detailStr = "";
            # Clamper les distances au max des tableaux de penalites / Clamp distances to max penalty array index
            my $maxPenIdx_r = $signatureMaxLength - 1;
            my $d_stem_r   = ($innerToStemDistance   < $maxPenIdx_r) ? $innerToStemDistance   : $maxPenIdx_r;
            my $d_middle_r = ($innerToMiddleDistance < $maxPenIdx_r) ? $innerToMiddleDistance : $maxPenIdx_r;
            my $d_outer_r  = ($middleToOuterDistance < $maxPenIdx_r) ? $middleToOuterDistance : $maxPenIdx_r;

            if($includeStemPrimers == $TRUE)
            {
              # STEM ARCHITECTURE: Calculate spacing between inner-STEM and inner-middle
              $spacingPenalty = 
                ($innerToLoopPenalties_r->[$d_stem_r]    * $innerToStemPenaltyWeight) +
                ($innerToMiddlePenalties_r->[$d_middle_r] * $innerToMiddlePenaltyWeight) +
                ($middleToOuterPenalties_r->[$d_outer_r]  * $middleToOuterPenaltyWeight);
              $primer3Penalty = 
                $innerPenalty * $innerPenaltyWeight +
                $stemPenalty  * $stemPenaltyWeight +
                $middlePenalty * $middlePenaltyWeight +
                $outerPenalty * $outerPenaltyWeight;
              
              $detailStr = sprintf("Spc[I_S:%.1f I_M:%.1f M_O:%.1f] Thm[I:%.1f S:%.1f M:%.1f O:%.1f]",
                    ($innerToLoopPenalties_r->[$d_stem_r]    * $innerToStemPenaltyWeight),
                    ($innerToMiddlePenalties_r->[$d_middle_r] * $innerToMiddlePenaltyWeight),
                    ($middleToOuterPenalties_r->[$d_outer_r]  * $middleToOuterPenaltyWeight),
                    ($innerPenalty * $innerPenaltyWeight),
                    ($stemPenalty  * $stemPenaltyWeight),
                    ($middlePenalty * $middlePenaltyWeight),
                    ($outerPenalty * $outerPenaltyWeight));
            }
            else
            {
              $spacingPenalty = 
                ($innerToMiddlePenalties_r->[$d_middle_r] * $innerToMiddlePenaltyWeight) +
                ($middleToOuterPenalties_r->[$d_outer_r]  * $middleToOuterPenaltyWeight);
              $primer3Penalty = 
                $innerPenalty * $innerPenaltyWeight +
                $middlePenalty * $middlePenaltyWeight +
                $outerPenalty * $outerPenaltyWeight;
              
              $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]",
                    ($innerToMiddlePenalties_r->[$d_middle_r] * $innerToMiddlePenaltyWeight),
                    ($middleToOuterPenalties_r->[$d_outer_r]  * $middleToOuterPenaltyWeight),
                    ($innerPenalty * $innerPenaltyWeight),
                    ($middlePenalty * $middlePenaltyWeight),
                    ($outerPenalty * $outerPenaltyWeight));
            }
 
            my $reverseSetPenalty = $spacingPenalty + $primer3Penalty;
            if($reverseSetPenalty < $bestSetPenalty)
            {
              $bestReverseInfos[$innerIndex] = [$stemInfo, $middleInfo, $outerInfo];
              $bestSetPenalty = $reverseSetPenalty;
              $_sig_rev_hits++;  # Compteur de signatures Rev / Rev signature counter
              $bestReversePenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty, $detailStr];
            }
          } # End reverse outer iteration
        } # End reverse middle iteration
      } # End reverse STEM iteration
    } # End reverse inner iteration

    # Finaliser la barre Reverse / Finalize Reverse bar
    printf(STDERR "\r%-80s\n", "") if $_LAVA_IS_TTY;
    print "  [Stem Rev] $_sig_rev_hits combinaisons Reverse trouvees sur $innerReverseCount amorces B1c.\n";

    ## Stop trying if no reverse primer sets were found (probably an un-needed optimization)
    #if($reverseSetCount == 0)
    #{
    #  print "R";
    #  next;
    #}

    # Now, try to combine forward and reverse primer sets into full signatures
    print "Combining Best F/R Halves to create LAMP Signatures...\n";
    my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner iteration
    for(my $i = 0; $i < $innerForwardCount; $i++)
    {
      # Skip inner primers without primer sets
      if(! exists($bestForwardInfos[$i]))
      {
	next;
      }

      my $finnerInfo = $innerForwardSubset_r->[$i];
      my ($fstemInfo, $fmiddleInfo, $fouterInfo) = @{$bestForwardInfos[$i]};
      my ($forwardSpacingPenalty, $forwardPrimer3Penalty, $forwardDetailStr) = 
        @{$bestForwardPenalties[$i]};

      my $forwardStart = $fouterInfo->getLocation();
      my $forwardEnd = $finnerInfo->getLocation() + $finnerInfo->getLength() - 1;

      # Used to bound the upper end of the inner iteration search
      my $maxReverseLocation = $forwardStart + $signatureMaxLength - 1;
    
      # Used to help bound the lower end of the inner iteration search
      my $previousCompatibleIndexFound = $FALSE;
      
      for(my $j = $previousFirstCompatibleIndex; $j < $innerReverseCount; $j++)
      {
        # Skip reverse primers without primer sets
        if(! exists($bestReverseInfos[$j]))
        {
	  next;
        }
	my $binnerInfo = $innerReverseSubset_r->[$j];
	my ($bstemInfo, $bmiddleInfo, $bouterInfo) = @{$bestReverseInfos[$j]};
        my ($reverseSpacingPenalty, $reversePrimer3Penalty, $reverseDetailStr) = 
	  @{$bestReversePenalties[$j]};

        my $reverseEnd = $bouterInfo->getLocation();
        my $reverseStart = $binnerInfo->getLocation() - $binnerInfo->getLength() + 1;
        
        #print "\n  Outer $reverseStart -> $reverseEnd";

        # Advance to the next compatible reverse primer by skipping all the
        # primers located too far 5' with respect to the forward primer
        if($previousCompatibleIndexFound == $FALSE)
        {
          if($reverseStart <= $forwardEnd)
          {
            next;
          }
          else
          {
            $previousFirstCompatibleIndex = $j;
            $previousCompatibleIndexFound = $TRUE;
          }
        }

	# Stop searching if the inner iteration bounds are exceeded
        if($reverseStart > $maxReverseLocation)
        {
          last;
        }
 
        # Enforce minimum inner spacing distance
	my $innerSpacing = $reverseStart - ($forwardEnd + 1);
        if($innerSpacing < $minInnerPairSpacing)
	{
	  next;
	}

        # VALIDATION COMPLÈTE D'ESPACEMENT POUR TOUS LES PRIMERS
        # Nouvelle logique qui vérifie TOUS les primers de la signature / New logic that verifies ALL primers of the signature
        my @forwardPrimers = ();
        my @reversePrimers = ();
        
        # Collecter tous les primers forward avec leurs noms
        $fouterInfo->{name} = 'F3';
        $fmiddleInfo->{name} = 'F2';
        $finnerInfo->{name} = 'F1';
        push @forwardPrimers, $fouterInfo;
        push @forwardPrimers, $fmiddleInfo;  
        push @forwardPrimers, $finnerInfo;
        if($includeStemPrimers == $TRUE) {
          $fstemInfo->{name} = 'FSTEM';
          push @forwardPrimers, $fstemInfo;
        }
        
        # Collecter tous les primers reverse avec leurs noms
        if($includeStemPrimers == $TRUE) {
          $bstemInfo->{name} = 'BSTEM';
          push @reversePrimers, $bstemInfo;
        }
        $binnerInfo->{name} = 'B1';
        $bmiddleInfo->{name} = 'B2';
        $bouterInfo->{name} = 'B3';
        push @reversePrimers, $binnerInfo;
        push @reversePrimers, $bmiddleInfo;
        push @reversePrimers, $bouterInfo;
        
        # Utiliser la validation complète d'espacement / Use complete spacing validation
        if (!validateCompleteSignatureSpacing(\@forwardPrimers, \@reversePrimers, $minPrimerSpacing)) {
            next;
        }
       
        
       
        # Enforce max signature length
        if($reverseEnd - ($forwardStart + 1) > $signatureMaxLength)
	{
	  next;
	}

        #print "*";
    
        # TODO: Spacing penalty should probably exclude minimum required spacings?
        my $innerSpacingPenalty = 
          ($innerToInnerPenalties_r->[$innerSpacing] *
	   $innerForwardToReversePenaltyWeight);

        my $totalPenalty = $forwardSpacingPenalty +
	  $innerSpacingPenalty +
	  $reverseSpacingPenalty +
	  $forwardPrimer3Penalty +
	  $reversePrimer3Penalty;

        my $innerPair = LLNL::LAVA::PrimerSet::PCRPair->new(
	  {
	    "forward_info" => $finnerInfo,
	    "reverse_info" => $binnerInfo,
	  });
        my $middlePair = LLNL::LAVA::PrimerSet::PCRPair->new(
	  {
	    "forward_info" => $fmiddleInfo,
	    "reverse_info" => $bmiddleInfo,
	  });
        my $outerPair = LLNL::LAVA::PrimerSet::PCRPair->new(
          {
            "forward_info" => $fouterInfo,
            "reverse_info" => $bouterInfo,
	  });
  
        my $innerPairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $innerPair,
          });
        my $middlePairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $middlePair,
          });
        my $outerPairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $outerPair,
          });
  
        my $signature = LLNL::LAVA::PrimerSet::LAMP->new(
          {
            "inner_info" => $innerPairInfo,
            "middle_info" => $middlePairInfo,
            "outer_info" => $outerPairInfo,
          });
  
        $signature->setTag("lamp_penalty", $totalPenalty);
  
        # Just for fine-tuning and debugging reports
        my $f_penalty_sum = $forwardSpacingPenalty + $forwardPrimer3Penalty;
        my $r_penalty_sum = $reverseSpacingPenalty + $reversePrimer3Penalty;
        $signature->setTag("penalty_notes", sprintf("Total F:%.1f R:%.1f | F{%s} | R{%s}", $f_penalty_sum, $r_penalty_sum, $forwardDetailStr, $reverseDetailStr));
  
        if($includeStemPrimers == $TRUE)
        {
          $signature->setTag("fstem_info", $fstemInfo);
          $signature->setTag("bstem_info", $bstemInfo);
          $signature->setTag("has_stem_primers", $TRUE);
        }
        
        # Tags STEM définis, la validation sera effectuée en bloc post-collecte
        # STEM tags set, validation will be done in the post-collection batch
         
        push(@{$allFoundSignatures_r}, $signature);
      } # End forward sets iteration
    } # End reverse sets iteration

  print "Found " .
    scalar(@{$allFoundSignatures_r}) .
    " total signatures across all iterations\n";

  if(scalar(@{$allFoundSignatures_r}) == 0)
  {
    print "Failed to identify signatures - exiting normally\n";
    exit(0);
  }

  # --- VALIDATION PAR SIGNATURE (Essential for correct tagging) ---
  # Validation individuelle de chaque signature avant reduction
  # Individual per-signature validation before reduction (mirrors LOOP behaviour)
  print "Validating and calculating coverage for " . scalar(@{$allFoundSignatures_r}) . " signatures...\n";
  my $validated_count = 0;
  foreach my $signature (@{$allFoundSignatures_r}) {
      # Recalcule l'intersection et met a jour les tags de couverture
      # Recalculate intersection and update coverage tags
      my ($amplified_seqs, $coverage, $status) = calculateSignatureIntersection(
          $signature,
          scalar(@sequences),
          $signatureCommonTargetMinPercent,
          $includeStemPrimers,
          "stem"
      );
      $validated_count++;
  }
  print "Validation complete.\n";

  # Sort signatures by score BEFORE reduction for the complete file
  my @allSignatures = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1] }
    map {[$_, $_->getTag("lamp_penalty")]}
    @{$allFoundSignatures_r};

  # Write ALL signatures to a separate file before reduction
  my $allSignaturesFileName = "$outputFileName.all_signatures";
  open(OUTALLSIGS, "> $allSignaturesFileName") ||
    confess("file error - failed to open output file \"$allSignaturesFileName\" " .
      "for writing: $!");

  print "Writing all " . scalar(@allSignatures) . " signatures to $allSignaturesFileName\n";

  for(my $i = 0; $i < scalar(@allSignatures); $i++)
  {
    my $signature = $allSignatures[$i];
    my $signatureName = "$i";
 
    # Set the linker to the dash, but restore the exsting linker afterwards
    my $originalLinker = $signature->linker();  
    $signature->linker("");

    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    my $coverage_percent = $signature->getTag("signature_coverage_percent") || "N/A";
    my $target_count = $signature->getTag("signature_target_count") || "N/A";

    my $penaltyNotes = $signature->getTagExists("penalty_notes") ? $signature->getTag("penalty_notes") : "";
    my $sigNum = $signatureName + 1;
    my $degenerate_bases = $signature->getTagExists("degenerate_bases") ? $signature->getTag("degenerate_bases") : 0;
    my $sigLength = $signature->getLength();
    my $headerLine = "Signature ${sigNum} (length: ${sigLength}bp) (penalty: $penalty) $penaltyNotes (coverage: ${target_count}seqs/${coverage_percent}%) (degenerate: ${degenerate_bases} bases) (locations: $locationSummary)";
    if($includeStemPrimers && $signature->getTagExists("fstem_info")) {
      my $stemLocationSummary = $signature->getStemLocationSummary();
      $headerLine .= " STEM (locations: $stemLocationSummary)";
    }
    print OUTALLSIGS "$headerLine\n";
    print OUTALLSIGS ">${sigNum}_F3\n";
    print OUTALLSIGS $signature->getF3() . "\n";
    print OUTALLSIGS ">${sigNum}_B3\n";
    print OUTALLSIGS $signature->getB3() . "\n";
    print OUTALLSIGS ">${sigNum}_F2\n";
    print OUTALLSIGS $signature->getF2() . "\n";
    print OUTALLSIGS ">${sigNum}_B2\n";
    print OUTALLSIGS $signature->getB2() . "\n";
    print OUTALLSIGS ">${sigNum}_F1\n";
    print OUTALLSIGS $signature->getF1() . "\n";
    print OUTALLSIGS ">${sigNum}_B1\n";
    print OUTALLSIGS $signature->getB1() . "\n";
            if($includeStemPrimers == $TRUE)
        {
          my $fstemSequence = ($signature->getTag("fstem_info"))->getSequence();
          my $bstemSequence = ($signature->getTag("bstem_info"))->getSequence();
              my $stemLocationSummary = $signature->getStemLocationSummary();
      print OUTALLSIGS ">${sigNum}_FSTEM\n";
          print OUTALLSIGS $fstemSequence . "\n";
          print OUTALLSIGS ">${sigNum}_BSTEM\n";
          print OUTALLSIGS $bstemSequence . "\n";
        }

    # Return the linker back to its original state
    $signature->linker($originalLinker);
  }

  close(OUTALLSIGS) ||
    confess("file error - failed to close output file \"$allSignaturesFileName\": $!");

  # Filtrer pour ne garder que les signatures validées (couverture >= seuil) pour le fichier principal .primers
  # Filter to only keep validated signatures (coverage >= threshold) for the main .primers file
  my @valid_sigs = grep { $_->getTag("validation_status") eq "VALIDEE" } @{$allFoundSignatures_r};
  $allFoundSignatures_r = \@valid_sigs;

  # NOW apply the overlap reduction for the main output files
  # Réduction finale des signatures par chevauchement / Final signature overlap reduction
  my $possibleSignatures_r = reduceSignaturesByOverlap(
    {
      "signatures" => $allFoundSignatures_r,
      "max_overlap_percent" => $maxSigOverlapPercent,
      "resolve_overlap_by" => $resolveOverlapBy,
    });

  print "After reduction: " .
    scalar(@{$possibleSignatures_r}) .
    " final signatures\n";

  # Sort signatures by Coverage (Desc) -> Degeneracy (Asc) -> Penalty (Asc)
  my @possibleSignatures = 
    map {$_->[0]}
    sort {
        # 1. Coverage (Descending)
        # Note: coverage percent is stored as a string "XX.XX", numerical comparison needed
        my $covA = $a->[0]->getTag("signature_coverage_percent") || 0;
        my $covB = $b->[0]->getTag("signature_coverage_percent") || 0;
        
        # 2. Degenerate Bases (Ascending)
        my $degA = $a->[2];
        my $degB = $b->[2];
        
        # 3. Penalty (Ascending)
        my $penA = $a->[1];
        my $penB = $b->[1];

        if ($covB != $covA) {
            return $covB <=> $covA;
        } elsif ($degA != $degB) {
            return $degA <=> $degB;
        } else {
            return $penA <=> $penB;
        }
    }
    map {
        my $sig = $_;
        my $penalty = $sig->getTag("lamp_penalty");
        
        # Calculate total degenerate bases for the full signature
        my $seqs = $sig->getF3() . $sig->getB3() . $sig->getFIP() . $sig->getBIP();
        if ($includeStemPrimers) {
            $seqs .= ($sig->getTag("fstem_info"))->getSequence();
            $seqs .= ($sig->getTag("bstem_info"))->getSequence();
        }
        my $degenerateCount = countDegenerateBases($seqs);
        $sig->setTag("degenerate_bases", $degenerateCount); # Store for output
        
        [$sig, $penalty, $degenerateCount]
    }
    @{$possibleSignatures_r};

  # Remplacer la référence globale par la nouvelle liste triée / Replace global reference with sorted list
  $possibleSignatures_r = \@possibleSignatures;

  # Analyser les combinaisons de signatures (SUR LES SIGNATURES RÉDUITES ET VALIDÉES)
  if (scalar(@possibleSignatures) > 0) {
    my $num_signatures = scalar(@possibleSignatures);
    print "\n🔍 Analyse de / Analysis ofs combinaisons sur les $num_signatures signatures finales après réduction...\n";
    
    # Vérifier d'abord si une signature atteint déjà 100% de couverture / First check if a signature already reaches 100% coverage
    my $has_perfect_signature = 0;
    my $max_coverage = 0;
    
    foreach my $signature (@possibleSignatures) {
      my $coverage = $signature->getTag("signature_coverage_percent") || 0;
      $max_coverage = $coverage if $coverage > $max_coverage;
      
      if ($coverage >= 100.0) {
        $has_perfect_signature = 1;
        last;
      }
    }
    
    print "🔍 Couverture maximale des signatures individuelles: ${max_coverage}%\n";
    
    if ($has_perfect_signature) {
      print "✅ Une ou plusieurs signatures atteignent déjà 100% de couverture.\n";
      print "   L'analyse des combinaisons n'est pas nécessaire.\n\n";
    } else {
      print "📊 Aucune signature n'atteint 100% - Lancement de l'analyse des combinaisons...\n\n";
      
      # Limiter l'analyse des combinaisons pour éviter les calculs trop longs / Limit combination analysis to avoid excessive computation time
      my $max_signatures_for_analysis = 15;  # Limite raisonnable
      my @signatures_to_analyze = @possibleSignatures;
      
      if ($num_signatures > $max_signatures_for_analysis) {
        print "⚠️  Trop de signatures ($num_signatures) pour l'analyse complète des combinaisons.\n";
        print "   Analyse limitée aux $max_signatures_for_analysis meilleures signatures (triées par couverture/dégénérescence).\n";
        
        # ELLES SONT DÉJÀ TRIÉES PAR NOTRE NOUVEAU CRITÈRE
        @signatures_to_analyze = @signatures_to_analyze[0 .. $max_signatures_for_analysis - 1];
      }
      
      my $combination_results = analyzeSignatureCombinations(\@signatures_to_analyze, scalar(@sequences));
      
      # Sauvegarder les résultats de combinaisons dans un fichier / Save combination results to a file
      my $outputFileBase = $outputFileName;
      $outputFileBase =~ s/\.(txt|fasta|fa)$//;  # Enlever l'extension si présente / Remove extension if present
      my $combinations_file = "${outputFileBase}_combinations.txt";
      open(my $comb_fh, '>', $combinations_file) or die "Cannot open $combinations_file: $!";
      
      print $comb_fh "ANALYSE DES COMBINAISONS DE SIGNATURES\n";
      print $comb_fh "=====================================\n\n";
      print $comb_fh "Nombre total de signatures: " . scalar(@signatures_to_analyze) . "\n";
      print $comb_fh "Nombre total de séquences: " . scalar(@sequences) . "\n\n";
      
      for my $size (sort {$a <=> $b} keys %{$combination_results}) {
        print $comb_fh "COMBINAISONS $size par $size:\n";
        print $comb_fh "=" x 30 . "\n";
        
        my $results = $combination_results->{$size};
        my $max_index = ($#{$results} < 9) ? $#{$results} : 9;
        for my $i (0 .. $max_index) {  # Top 10 pour chaque taille
          my $result = $results->[$i];
          my $names_str = join(" + ", @{$result->{signature_names}});
          printf $comb_fh "%2d. %s: %d séquences (%.2f%%)\n", 
                 $i + 1, $names_str, $result->{union_count}, $result->{union_coverage};
        }
        print $comb_fh "\n";
      }
      
      close($comb_fh);
      print "Analyse de / Analysis ofs combinaisons sauvegardée dans: $combinations_file\n\n";
    }
  }

  # Write the output fasta
  #TODO: watch for stompping?
  open(OUTANSWER, "> $outputFileName") ||
    confess("file error - failed to open output file \"$outputFileName\" " .
      "for writing: $!");

  my $possibleSignatureCount = scalar(@possibleSignatures);
  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
    my $signatureName = "$i";
 
    # Set the linker to the dash, but restore the exsting linker afterwards
    my $originalLinker = $signature->linker();  
    $signature->linker("");

    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();

    #print OUTANSWER ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTANSWER ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTagExists("penalty_notes") ? $signature->getTag("penalty_notes") : "";
    my $target_count = $signature->getTagExists("signature_target_count") ? $signature->getTag("signature_target_count") : 0;
    my $coverage_percent = $signature->getTagExists("signature_coverage_percent") ? $signature->getTag("signature_coverage_percent") : "0.00";
    my $degenerate_bases = $signature->getTagExists("degenerate_bases") ? $signature->getTag("degenerate_bases") : 0;
    
    my $sigNum = $signatureName + 1;
    my $sigLength = $signature->getLength();
    my $headerLine = "Signature ${sigNum} (length: ${sigLength}bp) (penalty: $penalty) $penaltyNotes (coverage: ${target_count}seqs/${coverage_percent}%) (degenerate: ${degenerate_bases} bases) (locations: $locationSummary)";
    if($includeStemPrimers && $signature->getTagExists("fstem_info")) {
      my $stemLocationSummary = $signature->getStemLocationSummary();
      $headerLine .= " STEM (locations: $stemLocationSummary)";
    }
    print OUTANSWER "$headerLine\n";
    print OUTANSWER ">${sigNum}_F3\n";
    print OUTANSWER $signature->getF3() . "\n";
    print OUTANSWER ">${sigNum}_B3\n";
    print OUTANSWER $signature->getB3() . "\n";
    print OUTANSWER ">${sigNum}_FIP\n";
    print OUTANSWER $signature->getFIP() . "\n";
    print OUTANSWER ">${sigNum}_BIP\n";
    print OUTANSWER $signature->getBIP() . "\n";
    if($includeStemPrimers == $TRUE)
    {
      my $fstemSequence = ($signature->getTag("fstem_info"))->getSequence();
      my $bstemSequence = ($signature->getTag("bstem_info"))->getSequence();
      my $stemLocationSummary = $signature->getStemLocationSummary();
      print OUTANSWER ">${sigNum}_FSTEM\n";
      print OUTANSWER $fstemSequence . "\n";
      print OUTANSWER ">${sigNum}_BSTEM\n";
      print OUTANSWER $bstemSequence . "\n";
    }

    # Return the linker back to its original state
    $signature->linker($originalLinker);
  }

  close(OUTANSWER) ||
    confess("file error - failed to cose output file \"$outputFileName\": $!");

  # Write the linked-marker file, using the dash-linker in context
  #TODO: watch for stompping?
  my $dashFileName = "$outputFileName.dash";
  open(OUTDASH, "> $dashFileName") || 
    confess("file error - failed to open output file \"$dashFileName\" " .
      "for writing: $!");

  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
    my $signatureName = "$i";
 
    # Set the linker to the dash, but restore the exsting linker afterwards
    my $originalLinker = $signature->linker();  
    $signature->linker("-");

    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();
     
    #print OUTDASH ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTDASH ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTagExists("penalty_notes") ? $signature->getTag("penalty_notes") : "";
    my $target_count = $signature->getTagExists("signature_target_count") ? $signature->getTag("signature_target_count") : 0;
    my $coverage_percent = $signature->getTagExists("signature_coverage_percent") ? $signature->getTag("signature_coverage_percent") : "0.00";
    my $degenerate_bases = $signature->getTagExists("degenerate_bases") ? $signature->getTag("degenerate_bases") : 0;

    my $sigNum = $signatureName + 1;
    my $sigLength = $signature->getLength();
    my $headerLine = "Signature ${sigNum} (length: ${sigLength}bp) (penalty: $penalty) $penaltyNotes (coverage: ${target_count}seqs/${coverage_percent}%) (degenerate: ${degenerate_bases} bases) (locations: $locationSummary)";
    if($includeStemPrimers && $signature->getTagExists("fstem_info")) {
      my $stemLocationSummary = $signature->getStemLocationSummary();
      $headerLine .= " STEM (locations: $stemLocationSummary)";
    }
    print OUTDASH "$headerLine\n";
    print OUTDASH ">${sigNum}_F3\n";
    print OUTDASH $signature->getF3() . "\n";
    print OUTDASH ">${sigNum}_B3\n"; 
    print OUTDASH $signature->getB3() . "\n";
    print OUTDASH ">${sigNum}_FIP\n"; 
    print OUTDASH $signature->getFIP() . "\n";
    print OUTDASH ">${sigNum}_BIP\n"; 
    print OUTDASH $signature->getBIP() . "\n";
    if($includeStemPrimers == $TRUE)
    {
      my $fstemSequence = ($signature->getTag("fstem_info"))->getSequence();
      my $bstemSequence = ($signature->getTag("bstem_info"))->getSequence();
      my $stemLocationSummary = $signature->getStemLocationSummary();
      print OUTDASH ">${sigNum}_FSTEM\n";
      print OUTDASH $fstemSequence . "\n";
      print OUTDASH ">${sigNum}_BSTEM\n";
      print OUTDASH $bstemSequence . "\n";
    }
      
    # Return the linker back to its original state
    $signature->linker($originalLinker);
  }

  close(OUTDASH) ||
    confess("file error - failed to cose output file \"$dashFileName\": $!");

  # Write the individual primers (in extension orientation) as an answer file
  #TODO: watch for stompping?
  my $primersFileName = "$outputFileName.primers";
  open(OUTPRIMERS, "> $primersFileName") || 
    confess("file error - failed to open output file \"$primersFileName\" " .
      "for writing: $!");

  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
 
    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();
     
    #print OUTPRIMERS ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTPRIMERS ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTagExists("penalty_notes") ? $signature->getTag("penalty_notes") : "";
    my $target_count = $signature->getTagExists("signature_target_count") ? $signature->getTag("signature_target_count") : 0;
    my $coverage_percent = $signature->getTagExists("signature_coverage_percent") ? $signature->getTag("signature_coverage_percent") : "0.00";
    my $degenerate_bases = $signature->getTagExists("degenerate_bases") ? $signature->getTag("degenerate_bases") : 0;

    my $sigNum = $i + 1;
    my $sigLength = $signature->getLength();
    my $headerLine = "Signature ${sigNum} (length: ${sigLength}bp) (penalty: $penalty) $penaltyNotes (coverage: ${target_count}seqs/${coverage_percent}%) (degenerate: ${degenerate_bases} bases) (locations: $locationSummary)";
    if($includeStemPrimers && $signature->getTagExists("fstem_info")) {
      my $stemLocationSummary = $signature->getStemLocationSummary();
      $headerLine .= " STEM (locations: $stemLocationSummary)";
    }
    print OUTPRIMERS "$headerLine\n";
    print OUTPRIMERS ">${sigNum}_F3\n";
    print OUTPRIMERS $signature->getF3() . "\n";
    print OUTPRIMERS ">${sigNum}_B3\n"; 
    print OUTPRIMERS $signature->getB3() . "\n";
    print OUTPRIMERS ">${sigNum}_F2\n";
    print OUTPRIMERS $signature->getF2() . "\n";
    print OUTPRIMERS ">${sigNum}_B2\n";
    print OUTPRIMERS $signature->getB2() . "\n";
    print OUTPRIMERS ">${sigNum}_F1\n";
    print OUTPRIMERS $signature->getF1() . "\n";
    print OUTPRIMERS ">${sigNum}_B1\n";
    print OUTPRIMERS $signature->getB1() . "\n";
    if($includeStemPrimers == $TRUE)
    {
      my $fstemSequence = ($signature->getTag("fstem_info"))->getSequence();
      my $bstemSequence = ($signature->getTag("bstem_info"))->getSequence();
      my $stemLocationSummary = $signature->getStemLocationSummary();
      print OUTPRIMERS ">${sigNum}_FSTEM\n";
      print OUTPRIMERS $fstemSequence . "\n";
      print OUTPRIMERS ">${sigNum}_BSTEM\n";
      print OUTPRIMERS $bstemSequence . "\n";
    }
  }

  close(OUTPRIMERS) ||
    confess("file error - failed to cose output file \"$primersFileName\": $!");

  # Générer les fichiers FASTA des séquences amplifiées et exclues / Generate FASTA files of amplified and excluded sequences
  my $output_base = $outputFileName;
  $output_base =~ s/\.(txt|fasta|fa)$//;
  createAmplificationFiles($possibleSignatures_r, \@sequence_objects, \@sequence_names, $output_base);

  # Creer les fichiers par signature individuelle (PipelineUtils unifie)
  # Create per-signature files (unified PipelineUtils)
  createPerSignatureFiles($possibleSignatures_r, \@sequence_names, $output_base, "STEM");

  print "Exiting normally\n";
}

# Les fonctions utilitaires partagées / The shared utility functions (buildReversePrimers, analyzeAll, enumeratePairs,
# buildMetricsArray, reducePairInfosByPenalty, reducePrimersByOverlap,
# reduceSignaturesByOverlap, flattenInfoData) sont désormais dans: / are now in:
# lib/LLNL/LAVA/PipelineUtils.pm
# Shared utility functions are now in: lib/LLNL/LAVA/PipelineUtils.pm

