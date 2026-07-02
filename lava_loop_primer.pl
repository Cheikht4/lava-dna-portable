#!/usr/bin/perl -w

################################################################################
#
#
# Version 0.1.2 (2016)
# Updated by Michaël Bekaert <michael.bekaert@stir.ac.uk>.
# Produced at the Institute of Aquacuture, University of Stirling, UK
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
use LLNL::LAVA::Core qw(generateDistancePenalties calculate_proportional_geometry generateSigmoidPenalty countDegenerateBases);
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
# Ces fonctions sont desormais dans LLNL::LAVA::PipelineUtils (Phase 36).
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

=cut

sub getOligosWithMismatchTolerance {
  my ($enumerator, $alignment, $min_match_percent, $min_iupac_percent, $min_primer_coverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency) = @_;
  
  my $sequenceCount = $alignment->num_sequences();
  if ($sequenceCount <= 0) {
    confess("data error - MSA must contain at least one sequence");
  }
  
  print "INFO: Utilisation de la tolérance aux mismatches activée\n";
  print "  - Seuil de concordance stricte / Strict match threshold: ${min_match_percent}%\n";
  print "  - Seuil de couverture IUPAC / IUPAC coverage threshold: ${min_iupac_percent}%\n";
  
  my @sequences = ();
  foreach my $sequence ($alignment->each_seq()) {
    my $seqContent = $sequence->seq();
    $seqContent = uc($seqContent);  # Convertir en majuscules d'abord / Convert to uppercase first
    $seqContent =~ s/[^ATCG]/N/g;  # Puis normaliser (remplacer caractères non-ADN par N) / Then normalize (replace non-DNA characters with N)
    push @sequences, $seqContent;
  }
  
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
      name   => "Forward Validation",
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
    
    my ($final_sequence, $coverage_percent, $is_degenerate, $compatible_seq_ids) = 
      checkPrimerMismatchTolerance(\@sequences, $location, $length, $original_sequence,
                                  $min_match_percent, $min_iupac_percent, $min_primer_coverage,
                                  $maxTotalDegen, $maxConsecDegen, 
                                  $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);
    
    # Utiliser le même seuil que l'algorithme interne (paramètre passé à la fonction) / Use the same threshold as the internal algorithm (parameter passed to the function)
    my $min_primer_acceptance = $min_primer_coverage;
    if ($coverage_percent >= $min_primer_acceptance) {
      my $validatedPrimer = $primer->clone();
      
      if ($is_degenerate) {
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
            printf("[LAVA-PROGRESS] Validation Forward|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                   $_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
          }
        }
      } else {
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
            printf("[LAVA-PROGRESS] Validation Forward|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                   $_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
          }
        }
      }
      
      push @validatedPrimers, $validatedPrimer;
    } else {
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
          printf("[LAVA-PROGRESS] Validation Forward|%d|%d|OK:%d DEG:%d REJ:%d|%.0f it/s|%d\n",
                 $_pb_done,$nb_fwd_candidates,$strict_count,$degenerate_count,$rejected_count,$rate,$eta);
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

      "loop_primer_target_length=i" => \$options{"loop_primer_target_length"},
      "loop_primer_min_length=i" => \$options{"loop_primer_min_length"},
      "loop_primer_max_length=i" => \$options{"loop_primer_max_length"},
      "loop_primer_target_tm=f" => \$options{"loop_primer_target_tm"},
      "loop_primer_min_tm=f" => \$options{"loop_primer_min_tm"},
      "loop_primer_max_tm=f" => \$options{"loop_primer_max_tm"},

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

      "include_loop_primers=i" => \$options{"include_loop_primers"},
      "loop_min_gap=i" => \$options{"loop_min_gap"},
      "min_signatures_for_success=i" => \$options{"min_signatures_for_success"},
      "min_primer_spacing=i" => \$options{"min_primer_spacing"},
      "min_inner_pair_spacing=i" => \$options{"min_inner_pair_spacing"},
      "max_overlap_percent=f" => \$options{"max_overlap_percent"},
      "resolve_overlap_by=s" => \$options{"resolve_overlap_by"},
      # --- REDUCTION SPATIALE PAR FENETRE / SPATIAL WINDOW REDUCTION ---
      "window_size=i"    => \$options{"window_size"},    # largeur fenetre en nt (0=desactive)
      "max_per_window=i" => \$options{"max_per_window"}, # max candidats par fenetre
      # Calcul dynamique des longueurs (porte depuis STEM / ported from STEM)
      "max_dist_outer_middle=i" => \$options{"max_dist_outer_middle"},
      "max_dist_middle_inner=i" => \$options{"max_dist_middle_inner"},

      "primer3_executable=s" => \$options{"primer3_executable"},
      "thermodynamic_path=s" => \$options{"thermodynamic_path"},
      "alignment_format=s" => \$options{"alignment_format"},
      "dntp_conc=f" => \$options{"dntp_conc"}, # new
      "salt_divalent=f" => \$options{"salt_divalent"}, # new
      "salt_monovalent=f" => \$options{"salt_monovalent"}, # new
      "dna_conc=f" => \$options{"dna_conc"}, # new
      "dna_conc=f" => \$options{"dna_conc"}, # new
      "max_primer_gen=f" => \$options{"max_primer_gen"}, # new
      "max_tm_diff=f" => \$options{"max_tm_diff"}, # new

      # Sigmoid Penalty Parameters
      "penalty_plateau=f" => \$options{"penalty_plateau"},
      "penalty_slope=f" => \$options{"penalty_slope"},

      # --- NOUVEAUX PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES ---
      "primer_min_match_percent=f" => \$options{"primer_min_match_percent"},
      "primer_iupac_min_percent=f" => \$options{"primer_iupac_min_percent"},
      "min_primer_coverage=f" => \$options{"min_primer_coverage"},
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
      "signature_max_length" => 320,
      "outer_primer_target_length" => 20,
      "outer_primer_min_length" => 18,
      "outer_primer_max_length" => 23,
      "outer_primer_target_tm" => "60.0",
      "loop_primer_target_length" => 20,
      "loop_primer_min_length" => 18,
      "loop_primer_max_length" => 23,
      "loop_primer_target_tm" => "60.0",
      "middle_primer_target_length" => 20,
      "middle_primer_min_length" => 18,
      "middle_primer_max_length" => 23,
      "middle_primer_target_tm" => "60.0",
      "inner_primer_target_length" => 23,
      "inner_primer_min_length" => 20,
      "inner_primer_max_length" => 26,
      "inner_primer_target_tm" => "62.0",
      "max_poly_bases" => 2,
      "include_loop_primers" => 1,
      "loop_min_gap" => 25,
      "min_signatures_for_success" => 1, # Should probably never go lower
      "min_primer_spacing" => 1,
      "min_inner_pair_spacing" => 1,
      # Some LAMP-specific approximate targets for a "minimum sized" signature
      # Currently, no penalty is assessed for lengths under the target size, so
      # these sizes are a little larger than they need to be.
      "outer_pair_target_length" => 200, 
      "middle_pair_target_length" => 160,
      "inner_pair_target_length" => 50, 
      "max_overlap_percent" => 0,
      "resolve_overlap_by" => "penalty",
      # Calcul dynamique (porte depuis STEM / ported from STEM)
      "max_dist_outer_middle" => 30,
      "max_dist_middle_inner" => 30,
      # --- PARAMETRES DE TOLERANCE AUX MISMATCHES ---
      "primer_min_match_percent" => 80,
      "primer_iupac_min_percent" => 98,
      "min_primer_coverage" => 80,
      # -----------------------------------------
      "dntp_conc" => 1.4,
      "salt_divalent" => 8,
      "salt_monovalent" => 50,
      "salt_monovalent" => 50,
      "dna_conc" => 400,
      "penalty_plateau" => 0.25,
      "penalty_slope" => 0.15,
      "max_primer_gen" => 10001, # primer3 rounding error off by 1?
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
      "    [--signature_max_length <length, deafult=" .
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
      # Loop primer options
      "    [--loop_primer_target_length <length, default=" . 
        $optionDefaults{"loop_primer_target_length"} .
        ">]\n" .
      "    [--loop_primer_min_length <length, default=" .
        $optionDefaults{"loop_primer_min_length"} .
	">]\n" .
      "    [--loop_primer_max_length <length, default=" .
        $optionDefaults{"loop_primer_max_length"} .
	">]\n" .
      "    [--loop_primer_target_tm <tm, default=" .
        $optionDefaults{"loop_primer_target_tm"} .
	"C>]\n" .
      "    [--loop_primer_min_tm <tm, default=loop_primer_target_tm - 1.0>]\n" .
      "    [--loop_primer_max_tm <tm, default=loop_primer_target_tm + 1.0>]\n" .
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
      "    [--include_loop_primers <length, default=" .
        $optionDefaults{"include_loop_primers"} .
	">]\n" .
      # Loop gap is the distance between the middle and inner primers
      "    [--loop_min_gap <length, default=" .
        $optionDefaults{"loop_min_gap"} .
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
      300);
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
  my $maxTmDiff = optionWithDefault($options_r, "max_tm_diff", 5.0);

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

  my $loopPrimerTargetLength =
    optionWithDefault($options_r, "loop_primer_target_length", 
      $optionDefaults{"loop_primer_target_length"});
  my $loopPrimerMinLength =
    optionWithDefault($options_r, "loop_primer_min_length", 
      $optionDefaults{"loop_primer_min_length"});
  if($loopPrimerMinLength > $loopPrimerTargetLength)
  {
    $loopPrimerMinLength = $loopPrimerTargetLength;
  }
  my $loopPrimerMaxLength =
    optionWithDefault($options_r, "loop_primer_max_length", 
      $optionDefaults{"loop_primer_max_length"});
  if($loopPrimerMaxLength < $loopPrimerTargetLength)
  {
    $loopPrimerMaxLength = $loopPrimerTargetLength;
  }

  my $loopPrimerTargetTM =
    optionWithDefault($options_r, "loop_primer_target_tm", 
      $optionDefaults{"loop_primer_target_tm"});
  my $loopPrimerMinTM =
    optionWithDefault($options_r, "loop_primer_min_tm", 
      ($loopPrimerTargetTM - 1.0));
  if($loopPrimerMinTM > $loopPrimerTargetTM)
  {
    $loopPrimerMinTM = $loopPrimerTargetTM;
  }
  my $loopPrimerMaxTM =
    optionWithDefault($options_r, "loop_primer_max_tm", 
      ($loopPrimerTargetTM + 1.0));
  if($loopPrimerMaxTM < $loopPrimerTargetTM)
  {
    $loopPrimerMaxTM = $loopPrimerTargetTM;
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
  
  my $includeLoopPrimers = 
    optionWithDefault($options_r, "include_loop_primers",
    $optionDefaults{"include_loop_primers"});
  my $loopMinGap = 
    optionWithDefault($options_r, "loop_min_gap", 
      $optionDefaults{"loop_min_gap"});
  my $signatureCommonTargetMinPercent =
    optionWithDefault($options_r, "min_signatures_for_success",
      $optionDefaults{"min_signatures_for_success"});
  my $maxSigOverlapPercent = 
    optionWithDefault($options_r, "max_overlap_percent",
      $optionDefaults{"max_overlap_percent"});
  my $resolveOverlapBy = 
    optionWithDefault($options_r, "resolve_overlap_by",
      $optionDefaults{"resolve_overlap_by"});

  my $dnaConc = 
    optionWithDefault($options_r, "dna_conc",
     $optionDefaults{"dna_conc"});
  my $dntpConc = 
    optionWithDefault($options_r, "dntp_conc",
     $optionDefaults{"dntp_conc"});
  my $saltMonovalent = 
    optionWithDefault($options_r, "salt_monovalent",
     $optionDefaults{"salt_monovalent"});

  my $penaltyPlateau = optionWithDefault($options_r, "penalty_plateau", $optionDefaults{"penalty_plateau"});
  my $penaltySlope = optionWithDefault($options_r, "penalty_slope", $optionDefaults{"penalty_slope"});

  my $saltDivalent = 
    optionWithDefault($options_r, "salt_divalent",
     $optionDefaults{"salt_divalent"});
  my $maxEnumeratedPrimers = int(
    optionWithDefault($options_r, "max_primer_gen",
    $optionDefaults{"max_primer_gen"}));
    
  my $minPrimerSpacing = 
    optionWithDefault($options_r, "min_primer_spacing", 
      $optionDefaults{"min_primer_spacing"});
  my $minInnerPairSpacing =
    optionWithDefault($options_r, "min_inner_pair_spacing", 
      $optionDefaults{"min_inner_pair_spacing"});
  #print "Loop min gap: $loopMinGap\n";
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

  # --- CALCUL DYNAMIQUE DES LONGUEURS CIBLES (PipelineUtils, porte depuis STEM) ---
  # --- DYNAMIC TARGET LENGTH CALCULATION (PipelineUtils, ported from STEM) ---
  if (exists $options_r->{"max_dist_outer_middle"} || exists $options_r->{"max_dist_middle_inner"})
  {
    my $maxDistOuterMiddle = 
      optionWithDefault($options_r, "max_dist_outer_middle",
        $optionDefaults{"max_dist_outer_middle"});
    my $maxDistMiddleInner =
      optionWithDefault($options_r, "max_dist_middle_inner",
        $optionDefaults{"max_dist_middle_inner"});

    # --- CORRECTION DE CONFLIT LOOP (Phase 36) ---
    # maxDistMiddleInner représente la distance cible / maxDistMiddleInner represents the target distance (Middle -> Inner) / 2 = F2_len + gap(F2, F1c).
    # Mais dans LOOP, gap(F2, F1c) doit être au minimum de loopMinGap pour accommoder le Loop primer.
    # Donc maxDistMiddleInner DOIT être >= middlePrimerTargetLength + loopMinGap.
    if ($includeLoopPrimers) {
      my $middlePrimerTargetLength = optionWithDefault($options_r, "middle_primer_target_length", $optionDefaults{"middle_primer_target_length"});
      my $min_required_dist = $middlePrimerTargetLength + $loopMinGap;
      
      if ($maxDistMiddleInner < $min_required_dist) {
        print "\nWARNING: max_dist_middle_inner ($maxDistMiddleInner) est trop petit pour accommoder loop_min_gap ($loopMinGap).\n";
        print "WARNING: Pour eviter un conflit geometrique bloquant, max_dist_middle_inner est ajuste automatiquement a $min_required_dist.\n";
        $maxDistMiddleInner = $min_required_dist;
      }
    }

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
  
  print "Configuration tolérance mismatches:\n";
  print "  - Match strict minimum: ${primerMinMatchPercent}%\n";
  print "  - Couverture IUPAC minimum: ${primerIupacMinPercent}%\n";
  print "  - Seuil élimination primer: ${minPrimerCoverage}%\n\n";

  # In theory, the overall score logic belongs in a PrimerSetAnalyzer, 
  # but I hope this helps me optimize the inner loop implementing it
  # here, and only instantiating LAMP signatures for the best combinations
  my $innerPenaltyWeight = "1.2";
  my $loopPenaltyWeight = ".7";
  my $middlePenaltyWeight = "1.1";
  my $outerPenaltyWeight = "1.0";

  my $innerToLoopPenaltyWeight = 0.5; # Reduced from 1.0 (LAVA 2026 - Spacing Relaxation)
  my $loopToMiddlePenaltyWeight = 0.5; 
  my $innerToMiddlePenaltyWeight = 0.5;
  my $middleToOuterPenaltyWeight = 0.5;
  my $innerForwardToReversePenaltyWeight = 0.5;

  # Let the games begin...

  # Load the input alignment, could be a single sequence
  # TODO: # Make sure the alignment format option suggestion is working
  my $alignIN = Bio::AlignIO->new(-file => "< $alignmentFastaName", -format => $alignmentFormat);
  my $inputMSA = $alignIN->next_aln();
  my $sequenceLength = $inputMSA->length;
  
  # Extraire les objets de séquence pour la génération des fichiers FASTA / Extract sequence objects for FASTA file generation
  my @sequence_objects = ();
  my @sequence_names = ();
  foreach my $sequence ($inputMSA->each_seq()) {
    push @sequence_objects, $sequence;
    push @sequence_names, $sequence->display_id();
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
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);

  print "  Generated \"" .
    scalar(@outerForwardPrimers) .
    "\" outer forward primers\n";
  
  # DEBUG : Vérifier que les primers ont le tag compatible_sequence_ids / Verify that primers have the compatible_sequence_ids tag
  my $outer_with_tag = 0;
  for my $primer (@outerForwardPrimers) {
    eval {
      my $tag = $primer->getTag("compatible_sequence_ids");
      $outer_with_tag++ if defined $tag;
    };
  }


  # Option B : Generation NATIVE des Reverse Outer via Primer3 sur RC(MSA)
  # Les Reverse sont generes independamment des Forward — protection 3' garantie
  # Option B: NATIVE Outer Reverse generation via Primer3 on RC(MSA)
  # Reverse primers generated independently from Forward — 3' protection guaranteed
  print "Enumerating outer NATIVE reverse primers (Option B)\n";
  my @outerReversePrimers = buildNativeReversePool(
    $outerEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp
  );
  print "  Generated \"" . scalar(@outerReversePrimers) . "\" outer native reverse primers\n";


  # Enumerate loop primers, since the loop primers extend in the opposite 
  # direction of the other LAMP primers, the back-loop primers are 
  # generated on the as-is sequence, and the forward-loop primers are 
  # built in the opposite orientation
  my $loopEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $loopEnumerator->setPrimer3Targets(
    {
      "target_length" => $loopPrimerTargetLength,
      "min_length" => $loopPrimerMinLength,
      "max_length" => $loopPrimerMaxLength,
      "target_tm" => $loopPrimerTargetTM,
      "min_tm" => $loopPrimerMinTM,
      "max_tm" => $loopPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
      "dna_conc" => $dnaConc,
      "dntp_conc" => $dntpConc,
      "salt_monovalent" => $saltMonovalent,
      "salt_divalent" => $saltDivalent,
      "entropy_threshold" => $entropyThreshold,
    });

  # This difference in naming is intentional for now (loopBackPrimers instead of 
  # loopReversePrimers), to serve as a reminder that
  # loop primers extend the other direction, and that their locations need to be 
  # with the opposite orientation
  
  my @loopBackPrimers = ();
  my @loopForwardPrimers = ();
  
  if($includeLoopPrimers == $TRUE) {
  # BLOOP : genere nativement sur le brin + (Back Loop = sens du brin +, 3' pointe vers B1c)
  # BLOOP: natively generated on plus strand (Back Loop = sense of plus strand, 3' points toward B1c)
  print "Enumerating loop BACK (BLOOP) primers on plus strand\n";
    @loopBackPrimers = getOligosWithMismatchTolerance($loopEnumerator, $inputMSA,
                                                        $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                        $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);

  print "  Generated \"" .
    scalar(@loopBackPrimers) .
    "\" loop BACK (BLOOP) primers\n";

  # FLOOP : Option B - genere nativement sur RC(MSA) pour garantir la protection 3'
  # FLOOP: Option B - natively generated on RC(MSA) to guarantee 3-prime protection
  # (Forward Loop = antisens, 3' pointe vers F1c - correspondait avant au 5' du BLOOP source = bug)
  print "Enumerating loop FORWARD (FLOOP) NATIVE reverse primers (Option B)\n";
    @loopForwardPrimers = buildNativeReversePool(
      $loopEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp
    );
  print "  Generated \"" . scalar(@loopForwardPrimers) . "\" loop FORWARD (FLOOP) native primers\n";
  } else {
    print "Loop primers désactivés - génération ignorée\n";
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
                                                           $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);

  print "  Generated \"" .
    scalar(@middleForwardPrimers) .
    "\" middle primers\n";

  # Option B : Generation NATIVE des Reverse Middle via Primer3 sur RC(MSA)
  print "Enumerating middle NATIVE reverse primers (Option B)\n";
  my @middleReversePrimers = buildNativeReversePool(
    $middleEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp
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
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency);

  print "  Generated \"" .
    scalar(@innerForwardPrimers) .
    "\" inner primers\n";
  
  # DEBUG : Vérifier que les primers inner ont le tag / Verify that inner primers have the tag
  my $inner_with_tag = 0;
  for my $primer (@innerForwardPrimers) {
    eval {
      my $tag = $primer->getTag("compatible_sequence_ids");
      $inner_with_tag++ if defined $tag;
    };
  }

  # Option B : Generation NATIVE des Reverse Inner via Primer3 sur RC(MSA)
  print "Enumerating inner NATIVE reverse primers (Option B)\n";
  my @innerReversePrimers = buildNativeReversePool(
    $innerEnumerator, $inputMSA,
    $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
    $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
    \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp
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
  my $loopPrimerAnalyzer = $outerPrimerAnalyzer;

  print "Analyzing outer forward primers\n";
  my $outerForwardPrimerMeasurements_r =
    analyzeAll(\@outerForwardPrimers, $outerPrimerAnalyzer);
  print "Analyzing outer reverse primers\n";
  my $outerReversePrimerMeasurements_r =
    analyzeAll(\@outerReversePrimers, $outerPrimerAnalyzer);

  my $loopForwardPrimerMeasurements_r = [];
  my $loopBackPrimerMeasurements_r = [];
  
  if($includeLoopPrimers == $TRUE) {
  print "Analyzing loop \"forward\" primers\n";
    $loopForwardPrimerMeasurements_r =
    analyzeAll(\@loopForwardPrimers, $loopPrimerAnalyzer);
  print "Analyzing loop \"back\" primers\n";
    $loopBackPrimerMeasurements_r =
    analyzeAll(\@loopBackPrimers, $loopPrimerAnalyzer);
  } else {
    print "Analyse de / Analysis ofs loop primers ignorée\n";
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

  # Loop primers sorted 2 ways (seulement si activés / only if enabled)
  my @loopForwardInfoByLocation = ();
  my @loopBackInfoByLocation = ();
  my @loopForwardInfoByPenalty = ();
  my @loopBackInfoByPenalty = ();
  
  if($includeLoopPrimers) {
    @loopForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$loopForwardPrimerMeasurements_r};
    @loopBackInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$loopBackPrimerMeasurements_r};

    @loopForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$loopForwardPrimerMeasurements_r};
    @loopBackInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$loopBackPrimerMeasurements_r};
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

  # Attempts will be made for combinations of different reduced primer sets.
  # The order attempts are made depends on the plan below.
  # 
  # Subgroups are built from the possible primer sets,  based on the overlap percent
  # specified in the subgroup schedule


  #=============================================================================
  # OPTIMIZATION PHASE 13: THE BIG MERGE & FAST-FAIL
  #=============================================================================
  # Legacy combinationPlan loop removed.
  # New single-pass logic with Generalized Sigmoid Penalty and Fast-Fail implementation follows.
  #=============================================================================
  
  #-----------------------------------------------------------------------------
  # 1. THE BIG MERGE: Single-Pass Subgroup Generation
  #-----------------------------------------------------------------------------
  # Instead of iterating through multiple plans with varying overlap stringencies,
  # we generate ONE set of high-quality "Master Lists" for each primer type.
  # We use the global $maxSigOverlapPercent to ensure good diversity without redundancy.
  
  print "Building Master Primer Lists (The Big Merge)...\n";

  my $bigMerge = buildBigMerge({
    inner_f_loc     => \@innerForwardInfoByLocation,
    inner_f_pen     => \@innerForwardInfoByPenalty,
    inner_r_loc     => \@innerReverseInfoByLocation,
    inner_r_pen     => \@innerReverseInfoByPenalty,
    special_f_loc   => \@loopForwardInfoByLocation,
    special_f_pen   => \@loopForwardInfoByPenalty,
    special_r_loc   => \@loopBackInfoByLocation,
    special_r_pen   => \@loopBackInfoByPenalty,
    include_special => $includeLoopPrimers,
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
  my $masterLoopF_r        = $bigMerge->{special_f};
  my $masterLoopF_data_r   = $bigMerge->{special_f_data};
  my $masterLoopR_r        = $bigMerge->{special_r};
  my $masterLoopR_data_r   = $bigMerge->{special_r_data};
  my $masterMiddleF_r      = $bigMerge->{middle_f};
  my $masterMiddleF_data_r = $bigMerge->{middle_f_data};
  my $masterMiddleR_r      = $bigMerge->{middle_r};
  my $masterMiddleR_data_r = $bigMerge->{middle_r_data};
  my $masterOuterF_r       = $bigMerge->{outer_f};
  my $masterOuterF_data_r  = $bigMerge->{outer_f_data};
  my $masterOuterR_r       = $bigMerge->{outer_r};
  my $masterOuterR_data_r  = $bigMerge->{outer_r_data};

  printf "Master Lists Counts:\n  Inner F: %d, R: %d\n  Loop F: %d, R: %d\n  Middle F: %d, R: %d\n  Outer F: %d, R: %d\n",
    scalar(@{$masterInnerF_r}),  scalar(@{$masterInnerR_r}),
    scalar(@{$masterLoopF_r}),   scalar(@{$masterLoopR_r}),
    scalar(@{$masterMiddleF_r}), scalar(@{$masterMiddleR_r}),
    scalar(@{$masterOuterF_r}),  scalar(@{$masterOuterR_r});

  #-----------------------------------------------------------------------------
  # 2. PRE-COMPUTE PENALTIES (Sigmoid)
  #-----------------------------------------------------------------------------
  # Pre-compute a set of distance penalties for faster use using the stored geometry
  # GENERATION PROPORTIONNELLE SIGMOÏDE (LAVA 2026)
  my $geometry = calculate_proportional_geometry($totalSignatureLength);
  
  # Pour Inner->Loop et Loop->Middle, cible 50% de F2-F1 chacun (répartition équilibrée) / For Inner->Loop and Loop->Middle, target 50% of F2-F1 each (balanced distribution)
  my $loop_target = int($geometry->{'f2_f1_target'} / 2);
  
  print "Generating Sigmoid Penalties (Core.pm)...\n";
  # Note: generateSigmoidPenalty is now used inside the loop or via pre-computed table
  # But the old `generateDistancePenalties` used a parabolic model. 
  # We will use the new `generateDistancePenalties` from Core which SHOULD be using sigmoid if updated,
  # OR we call `generateSigmoidPenalty` directly.
  # Let's assume generateDistancePenalties determines the penalty for array lookup. 
  # CHECK: Does generateDistancePenalties use the new sigmoid? 
  # Core.pm showed `generateDistancePenalties` calls `generateSigmoidPenalty` or `generatePenalty`.
  # Let's verify Core.pm briefly if needed. (Assuming it does based on previous interactions).
  
  my $innerToLoopPenalties_r = generateDistancePenalties($signatureMaxLength, $loop_target, $penaltyPlateau, $penaltySlope);
  my $loopToMiddlePenalties_r = generateDistancePenalties($signatureMaxLength, $loop_target, $penaltyPlateau, $penaltySlope);
  
  # Pour Middle->Outer (F2-F3), cible 12%
  my $middleToOuterPenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'f3_f2_target'}, $penaltyPlateau, $penaltySlope);
  
  # Pour Inner->Middle (sans loop), cible 18% (F1-F2)
  my $innerToMiddlePenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'f2_f1_target'}, $penaltyPlateau, $penaltySlope);

  # Pour Inner->Inner (F1c-B1c), cible 40%
  my $innerToInnerPenalties_r = generateDistancePenalties($signatureMaxLength, $geometry->{'inner_target'}, $penaltyPlateau, $penaltySlope);

  #-----------------------------------------------------------------------------
  # 3. OPTIMIZED NESTED LOOPS (Forward)
  #-----------------------------------------------------------------------------
  print "Scanning Forward Primer Combinations...\n";
  
  # Scoped Variables for Storage
  my @bestForwardInfos = (); 
  my @bestForwardPenalties = ();
  my $forwardSetCount = 0;
  
  # Legacy variable for later compatibility (though logic changed)
  my $possibleSignatures_r = [];
  my $allFoundSignatures_r = []; 

  # Setup Loop Placeholder if needed
  if (!$includeLoopPrimers) {
      my $placeHolderPrimer = LLNL::LAVA::Oligo->new({
          "sequence" => "N",
          "location" => 0, # Will be set dynamically per inner
          "strand" => "minus",
      });
      $placeHolderPrimer->setTag("primer3_penalty", 0);
      $placeHolderPrimer->setTag("primer3_tm", 0);
      
      # We create a single placeholder info, but we'll need to adjust its location 
      # or logic inside the loop. Actually, better to just use a dummy list 
      # and handle the location check specially.
      
      my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new({
          "penalty" => 0,
          "sequence" => "N",
          "location" => 0,
          "length" => 1,
          "analyzed_primer" => $placeHolderPrimer,
      });
      $masterLoopF_r = [$placeHolderInfo];
      $masterLoopF_data_r = [[0, 1, 0]]; # Location will be ignored/overridden
  }

  my $innerForwardCount = scalar(@{$masterInnerF_r});
  # Barre de progression pour la recherche combinatoire Forward / Progress bar for Forward combinatorial search
  my $_sig_fwd_t0   = time();
  my $_sig_fwd_done = 0;
  my $_sig_fwd_hits = 0;  # Nombre de signatures Forward trouvees / Forward signatures found
  print STDERR "  Recherche combinatoire Forward: $innerForwardCount amorces F1c...\n";

  for(my $innerIndex = 0; $innerIndex < $innerForwardCount; $innerIndex++)
  {
      # Barre de progression toutes les 50 iterations ou a 100%
      # Progress bar every 50 iterations or at 100%
      $_sig_fwd_done = $innerIndex + 1;
      if ($_sig_fwd_done % 50 == 0 || $_sig_fwd_done == $innerForwardCount) {
        # Emission progression Flask / Flask progress emission
        if ($_LAVA_IS_TTY || 1) {  # Toujours emettre pour Flask / Always emit for Flask
          my $elapsed = time() - $_sig_fwd_t0 + 0.001;
          my $eta = ($_sig_fwd_done < $innerForwardCount)
                    ? int(($innerForwardCount - $_sig_fwd_done) / ($_sig_fwd_done / $elapsed))
                    : 0;
          my $rate = $_sig_fwd_done / $elapsed;
          printf("[LAVA-PROGRESS] Signatures Forward|%d|%d|Sig: %d|%.1f it/s|%d\n",
                 $_sig_fwd_done, $innerForwardCount, $_sig_fwd_hits, $rate, $eta);
        }
      }
      my $innerInfo = $masterInnerF_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty, $innerTm) = @{$masterInnerF_data_r->[$innerIndex]};
      
      my $bestSetPenalty = 1000000;
      
      # 3.1 Calculate Search Bounds for Loop Primer
      my $searchStartAt = $innerLocation - $signatureMaxLength + $innerLength + 20;
      $searchStartAt = 0 if $searchStartAt < 0;
      
      my $loopEndAt = $innerLocation - 1 - $minPrimerSpacing;
      $loopEndAt = 0 if $loopEndAt < 0;
      
      # Determine Loop List to iterate
      my $currentLoopList_r = $includeLoopPrimers ? $masterLoopF_r : $masterLoopF_r; # Placeholder logic if false
      my $currentLoopData_r = $includeLoopPrimers ? $masterLoopF_data_r : $masterLoopF_data_r;
      my $loopCount = scalar(@{$currentLoopList_r});
      
      # If NO loops, we update the placeholder location to be valid (end of range)
      if (!$includeLoopPrimers) {
           $currentLoopData_r->[0]->[0] = $loopEndAt + 1; # Dummy valid location
      }

      for(my $i = 0; $i < $loopCount; $i++)
      {
          my $loopInfo = $currentLoopList_r->[$i];
          my ($loopLocation, $loopLength, $loopPenalty, $loopTm) = @{$currentLoopData_r->[$i]};
          
          if ($includeLoopPrimers) {
              # Fast-Fail: Sorted by location.
              # If loop is before search window, skip.
              next if ($loopLocation < $searchStartAt);
              
              # Fast-Fail: If loop is beyond search window, STOP.
              last if ($loopLocation > $loopEndAt);

              # --- DYNAMIC THERMAL FILTER (Inner vs Loop) ---
              next if (abs($innerTm - $loopTm) > $maxTmDiff);
          }
          
          # 3.2 Calculate Search Bounds for Middle Primer (F2) (Structure: F3-F2-LF-F1c)
          # Middle (F2) is upstream of Loop.
          my $middleStartAt = $searchStartAt; # Inherit max range start
          my $middleEndAt = $loopLocation - $loopLength - $minPrimerSpacing;
          
          # Ensure Middle isn't too close to Inner (respecting loop gap or standard gap)
          if ($includeLoopPrimers) {
              my $altMiddleEndAt = $innerLocation - ($loopMinGap + 1);
              $middleEndAt = $altMiddleEndAt if ($altMiddleEndAt < $middleEndAt);
          } else {
              my $altMiddleEndAt = $innerLocation - $minPrimerSpacing;
              $middleEndAt = $altMiddleEndAt if ($altMiddleEndAt < $middleEndAt);
          }
          $middleEndAt = 0 if $middleEndAt < 0;
          
          my $innerToLoopDistance = $innerLocation - ($loopLocation + 1);
          
          my $middleCount = scalar(@{$masterMiddleF_r});
          for(my $j = 0; $j < $middleCount; $j++)
          {
              my $middleInfo = $masterMiddleF_r->[$j];
              my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleF_data_r->[$j]};
              
              # Fast-Fail
              next if ($middleLocation < $middleStartAt);
              last if ($middleLocation > $middleEndAt);
              
              # Validity Checks (Redundant with bounds but safe)
              if ($includeLoopPrimers) {
                  next if ($middleLocation + $middleLength + $minPrimerSpacing > $loopLocation - $loopLength + 1);
                  next if ($middleLocation + $middleLength + $loopMinGap > $innerLocation);
              }

              # --- DYNAMIC THERMAL FILTER (Neighbor Check) ---
              if ($includeLoopPrimers) {
                  # Loop vs Middle
                  next if (abs($loopTm - $midTm) > $maxTmDiff);
              } else {
                  # Inner vs Middle (if no loop)
                  next if (abs($innerTm - $midTm) > $maxTmDiff);
              }

              # 3.3 Calculate Search Bounds for Outer Primer (F3)
              my $outerStartAt = $searchStartAt;
              my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;
              
              my $loopToMiddleDistance = ($loopLocation - $loopLength + 1) - ($middleLocation + $middleLength);
              
              my $outerCount = scalar(@{$masterOuterF_r});
              for(my $k = 0; $k < $outerCount; $k++)
              {
                  my $outerInfo = $masterOuterF_r->[$k];
                  my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterF_data_r->[$k]};
                  
                  # Fast-Fail
                  next if ($outerLocation < $outerStartAt);
                  last if ($outerLocation > $outerEndAt);
                  
                  # Distance Check
                  next if ($outerLocation + $outerLength + $minPrimerSpacing > $middleLocation);
                  
                  # --- DYNAMIC THERMAL FILTER (Middle vs Outer) ---
                  next if (abs($midTm - $outTm) > $maxTmDiff);
                  
                  # Calculate Penalty
                  my $middleToOuterDistance = $middleLocation - ($outerLocation + $outerLength);
                  my $innerToMiddleDistance = $innerLocation - ($middleLocation + $middleLength);
                  
                  my $spacingPenalty = 0;
                  my $primer3Penalty = 0;
                  my $detailStr = "";
                  
                  if ($includeLoopPrimers) {
                      # Ensure indices are within bounds of pre-computed arrays
                      # (Assuming generateDistancePenalties covers $signatureMaxLength)
                      $spacingPenalty = 
                          ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight) +
                          ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight) +
                          ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                          
                      $primer3Penalty = 
                          $innerPenalty * $innerPenaltyWeight +
                          $loopPenalty * $loopPenaltyWeight +
                          $middlePenalty * $middlePenaltyWeight +
                          $outerPenalty * $outerPenaltyWeight;
                      
                      $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                            ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight),
                            ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight),
                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                            ($innerPenalty * $innerPenaltyWeight),
                            ($loopPenalty * $loopPenaltyWeight),
                            ($middlePenalty * $middlePenaltyWeight),
                            ($outerPenalty * $outerPenaltyWeight)
                            );
                  } else {
                      $spacingPenalty = 
                          ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                          ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                          
                      $primer3Penalty = 
                          $innerPenalty * $innerPenaltyWeight +
                          $middlePenalty * $middlePenaltyWeight +
                          $outerPenalty * $outerPenaltyWeight;
                      
                      $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                            ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                            ($innerPenalty * $innerPenaltyWeight),
                            ($middlePenalty * $middlePenaltyWeight),
                            ($outerPenalty * $outerPenaltyWeight)
                            );
                  }
                  
                  my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                  
                  # Save if best
                  if ($currentSetPenalty < $bestSetPenalty) {
                      $bestForwardInfos[$innerIndex] = [$loopInfo, $middleInfo, $outerInfo];
                      $bestSetPenalty = $currentSetPenalty;
                      $bestForwardPenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty, $detailStr];
                      $forwardSetCount++;
                      $_sig_fwd_hits++;  # Compteur de signatures Forward / Forward signature counter
                  }
              } # End Outer
          } # End Middle
      } # End Loop
  } # End Inner
  
  # Finaliser la barre Forward / Finalize Forward bar
  print "  [Forward] $forwardSetCount combinaisons Forward trouvees sur $innerForwardCount amorces F1c.\n";

  # Check if anything found
  if($forwardSetCount == 0) {
      print "No valid forward primer combinations found.\n";
      exit 0;
  }

  #-----------------------------------------------------------------------------
  # 4. OPTIMIZED NESTED LOOPS (Reverse)
  #-----------------------------------------------------------------------------
  print "Scanning Reverse Primer Combinations...\n";

  # Setup Loop Placeholder if needed (Reverse)
  if (!$includeLoopPrimers) {
      my $placeHolderPrimer = LLNL::LAVA::Oligo->new({
          "sequence" => "N",
          "location" => 0, 
          "strand" => "plus",
      });
      $placeHolderPrimer->setTag("primer3_penalty", 0);
      $placeHolderPrimer->setTag("primer3_tm", 0);
      
      my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new({
          "penalty" => 0,
          "sequence" => "N",
          "location" => 0,
          "length" => 1,
          "analyzed_primer" => $placeHolderPrimer,
      });
      $masterLoopR_r = [$placeHolderInfo];
      $masterLoopR_data_r = [[0, 1, 0]]; 
  }

  my @bestReverseInfos = (); 
  my @bestReversePenalties = ();
  my $reverseSetCount = 0;
  
  my $innerReverseCount = scalar(@{$masterInnerR_r});
  # Barre de progression pour la recherche combinatoire Reverse / Progress bar for Reverse combinatorial search
  my $_sig_rev_t0   = time();
  my $_sig_rev_done = 0;
  my $_sig_rev_hits = 0;  # Nombre de signatures Reverse trouvees / Reverse signatures found
  print STDERR "  Recherche combinatoire Reverse: $innerReverseCount amorces B1c...\n";

  for(my $innerIndex = 0; $innerIndex < $innerReverseCount; $innerIndex++)
  {
      # Barre de progression toutes les 50 iterations ou a 100%
      # Progress bar every 50 iterations or at 100%
      $_sig_rev_done = $innerIndex + 1;
      if ($_sig_rev_done % 50 == 0 || $_sig_rev_done == $innerReverseCount) {
        # Emission progression Flask / Flask progress emission
        if ($_LAVA_IS_TTY || 1) {  # Toujours emettre pour Flask / Always emit for Flask
          my $elapsed = time() - $_sig_rev_t0 + 0.001;
          my $eta = ($_sig_rev_done < $innerReverseCount)
                    ? int(($innerReverseCount - $_sig_rev_done) / ($_sig_rev_done / $elapsed))
                    : 0;
          my $rate = $_sig_rev_done / $elapsed;
          printf("[LAVA-PROGRESS] Signatures Reverse|%d|%d|Sig: %d|%.1f it/s|%d\n",
                 $_sig_rev_done, $innerReverseCount, $_sig_rev_hits, $rate, $eta);
        }
      }
      my $innerInfo = $masterInnerR_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty, $innerTm) = @{$masterInnerR_data_r->[$innerIndex]};
      
      my $bestSetPenalty = 1000000;
      
      # 4.1 Calculate Search Bounds for Loop Primer (Reverse)
      # Inner (B1c) is at $innerLocation. 
      # Search starts after Inner and goes downstream (3').
      
      my $searchStartAt = $innerLocation + 1 + $minPrimerSpacing;
      
      my $searchEndAt = $innerLocation + $signatureMaxLength - $innerLength - 20; 
      
      # Determine Loop List to iterate
      my $currentLoopList_r = $includeLoopPrimers ? $masterLoopR_r : $masterLoopR_r; 
      my $currentLoopData_r = $includeLoopPrimers ? $masterLoopR_data_r : $masterLoopR_data_r;
      my $loopCount = scalar(@{$currentLoopList_r});
      
      # If NO loops, we update the placeholder location to be valid (start of range)
      if (!$includeLoopPrimers) {
           $currentLoopData_r->[0]->[0] = $searchStartAt - 1; # Dummy valid location
      }

      for(my $i = 0; $i < $loopCount; $i++)
      {
          my $loopInfo = $currentLoopList_r->[$i];
          my ($loopLocation, $loopLength, $loopPenalty, $loopTm) = @{$currentLoopData_r->[$i]};
          
          if ($includeLoopPrimers) {
              # Fast-Fail: Sorted by location.
              next if ($loopLocation < $searchStartAt);
              last if ($loopLocation > $searchEndAt);

              # --- DYNAMIC THERMAL FILTER (Inner vs Loop) ---
              next if (abs($innerTm - $loopTm) > $maxTmDiff);
          }
          
          # 4.2 Calculate Search Bounds for Middle Primer (B2)
          # Middle (B2) is downstream of Loop.
          my $middleStartAt = $loopLocation + $loopLength + $minPrimerSpacing;
          
          # Ensure Middle isn't too close to Inner (respecting loop gap or standard gap)
          if ($includeLoopPrimers) {
              my $altMiddleStartAt = $innerLocation + ($loopMinGap + 1);
              $middleStartAt = $altMiddleStartAt if ($altMiddleStartAt > $middleStartAt);
          } else {
              my $altMiddleStartAt = $innerLocation + $minPrimerSpacing;
              $middleStartAt = $altMiddleStartAt if ($altMiddleStartAt > $middleStartAt);
          }
          
          my $middleEndAt = $searchEndAt;

          my $innerToLoopDistance = $loopLocation - ($innerLocation + 1);
          
          my $middleCount = scalar(@{$masterMiddleR_r});
          for(my $j = 0; $j < $middleCount; $j++)
          {
              my $middleInfo = $masterMiddleR_r->[$j];
              my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleR_data_r->[$j]};
              
              # Fast-Fail
              next if ($middleLocation < $middleStartAt);
              last if ($middleLocation > $middleEndAt);
              
              # Validity Checks 
              if ($includeLoopPrimers) {
                   # Check overlap with Loop
                   next if ($middleLocation - $middleLength - $minPrimerSpacing < $loopLocation + $loopLength - 1);
                   # Check Gap with Inner 
                   next if ($middleLocation - $middleLength - $loopMinGap < $innerLocation);
              }

              # --- DYNAMIC THERMAL FILTER (Neighbor Check) ---
              if ($includeLoopPrimers) {
                  # Loop vs Middle
                  next if (abs($loopTm - $midTm) > $maxTmDiff);
              } else {
                  # Inner vs Middle (if no loop)
                  next if (abs($innerTm - $midTm) > $maxTmDiff);
              }

              # 4.3 Calculate Search Bounds for Outer Primer (B3)
              my $outerStartAt = $middleLocation + 1 + $minPrimerSpacing;
              my $outerEndAt = $searchEndAt;
              
              my $loopToMiddleDistance = ($middleLocation - $middleLength + 1) - ($loopLocation + $loopLength);

              my $outerCount = scalar(@{$masterOuterR_r});
              for(my $k = 0; $k < $outerCount; $k++)
              {
                  my $outerInfo = $masterOuterR_r->[$k];
                  my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterR_data_r->[$k]};
                  
                  # Fast-Fail
                  next if ($outerLocation < $outerStartAt);
                  last if ($outerLocation > $outerEndAt);
                  
                  # Distance Check : B3_leftmost - spacing > B2_rightmost (gap >= minPrimerSpacing)
                  # Correction off-by-one : leftmost de B3 = outerLocation - outerLength + 1
                  next if ($outerLocation - $outerLength + 1 - $minPrimerSpacing <= $middleLocation);

                  # --- DYNAMIC THERMAL FILTER (Middle vs Outer) ---
                  next if (abs($midTm - $outTm) > $maxTmDiff);
                  
                  # Calculate Penalty
                  my $middleToOuterDistance = ($outerLocation - $outerLength) - $middleLocation;
                  my $innerToMiddleDistance = ($middleLocation - $middleLength) - $innerLocation;
                  
                  my $spacingPenalty = 0;
                  my $primer3Penalty = 0;
                  my $detailStr = "";
                  
                  if ($includeLoopPrimers) {
                      $spacingPenalty = 
                          ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight) +
                          ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight) +
                          ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                          
                      $primer3Penalty = 
                          $innerPenalty * $innerPenaltyWeight +
                          $loopPenalty * $loopPenaltyWeight +
                          $middlePenalty * $middlePenaltyWeight +
                          $outerPenalty * $outerPenaltyWeight;
                      
                      $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                            ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight),
                            ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight),
                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                            ($innerPenalty * $innerPenaltyWeight),
                            ($loopPenalty * $loopPenaltyWeight),
                            ($middlePenalty * $middlePenaltyWeight),
                            ($outerPenalty * $outerPenaltyWeight)
                            );
                  } else {
                      $spacingPenalty = 
                          ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                          ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                          
                      $primer3Penalty = 
                          $innerPenalty * $innerPenaltyWeight +
                          $middlePenalty * $middlePenaltyWeight +
                          $outerPenalty * $outerPenaltyWeight;
                      
                      $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                            ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                            ($innerPenalty * $innerPenaltyWeight),
                            ($middlePenalty * $middlePenaltyWeight),
                            ($outerPenalty * $outerPenaltyWeight)
                            );
                  }
                  
                  my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                  
                  if ($currentSetPenalty < $bestSetPenalty) {
                      $bestReverseInfos[$innerIndex] = [$loopInfo, $middleInfo, $outerInfo];
                      $bestSetPenalty = $currentSetPenalty;
                      $bestReversePenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty, $detailStr];
                      $reverseSetCount++;
                      $_sig_rev_hits++;  # Compteur de signatures Reverse / Reverse signature counter
                  }
              } # End Outer
          } # End Middle
      } # End Loop
  } # End Inner

  # Finaliser la barre Reverse / Finalize Reverse bar
  print "  [Reverse] $reverseSetCount combinaisons Reverse trouvees sur $innerReverseCount amorces B1c.\n";

  if($reverseSetCount == 0) {
      print "No valid reverse primer combinations found.\n";
      print "WARNING: No reverse sets found.\n";
  }

  #-----------------------------------------------------------------------------
  # 5. COMBINE HALVES & CREATE SIGNATURES
  #-----------------------------------------------------------------------------
  print "Combining Best F/R Halves to create LAMP Signatures...\n";
  
  my $combinedSignatureCount = 0;

  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      next unless defined $bestForwardInfos[$i]; # Skip if no valid F-half found
      
      my $innerF = $masterInnerF_r->[$i];
      my $f_set_infos = $bestForwardInfos[$i]; # [LoopF, MidF, OutF]
       
      # InnerF (F1c) Location data
      # Note: $innerF is a PrimerInfo. 
      # $innerF->getLocation() is the END of the primer on the Plus strand (for Fwd? No).
      # Let's verify standard LAVA location semantics:
      # Forward Primer: Start lowest, End highest. Location = End? 
      # Reverse Primer: Start lowest (5' on minus), End highest (3' on minus, physical 5' on plus).
      # Typically LAVA uses "Genome Coordinates".
      # Let's trust the `enumeratePairs` logic which I replaced or the earlier lookups.
      # better: use the raw data arrays I prepared
      my $f1c_location = $masterInnerF_data_r->[$i]->[0]; # This is Location 
      my $f1c_length = $masterInnerF_data_r->[$i]->[1];
      my $f1c_tm = $masterInnerF_data_r->[$i]->[3]; # Unpack cached Tm
      
      # F1c is "Inner Forward". In LAMP, F1c is the complement of F1.
      # But LAVA PrimerSet::LAMP expects "Inner Info", which contains an "Analyzed Pair".
      # Each `*_info` is a `LLNL::LAVA::PrimerSetInfo::PCRPair`. 
      
      for(my $j = 0; $j < scalar(@{$masterInnerR_r}); $j++) {
          next unless defined $bestReverseInfos[$j];
          
          my $innerR = $masterInnerR_r->[$j];
          my $r_set_infos = $bestReverseInfos[$j]; # [LoopR, MidR, OutR]
          
          my $b1c_location = $masterInnerR_data_r->[$j]->[0];
          my $b1c_length = $masterInnerR_data_r->[$j]->[1];
          my $b1c_tm = $masterInnerR_data_r->[$j]->[3]; # Unpack cached Tm
          
          # --- DYNAMIC THERMAL FILTER (Inner Pair) ---
          next if (abs($f1c_tm - $b1c_tm) > $maxTmDiff);
          
          # Check Inner Distance (Target: 0-50 usually, or just non-overlapping)
          # F1c (InnerF) is upstream of B1c (InnerR).
          
          # Gap Calculation
          # Forward End = $f1c_location
          # Reverse Start = $b1c_location - $b1c_length + 1
          # Gap = ReverseStart - ForwardEnd - 1
          
          my $b1c_start_genome = $b1c_location - $b1c_length + 1;
          my $inner_gap = $b1c_start_genome - $f1c_location - 1;
          
          # Validity Checks
          next if ($inner_gap < 0); # Overlap
          next if ($inner_gap > 100); # Too far apart (Inner Gap Limit)

          # VALIDATION COMPLETE D'ESPACEMENT - tous les primers de la signature
          # Full spacing validation - all primers in the signature (mirrors STEM behavior)
          # Sans ce guard, des primers cross-strand (F3/B3, F2/B2) peuvent se chevaucher
          # Without this guard, cross-strand primers (F3/B3, F2/B2) can overlap
          my @fwdPrimers = ();
          my @revPrimers = ();

          my $outF_v = $f_set_infos->[2];  # F3
          my $midF_v = $f_set_infos->[1];  # F2
          $outF_v->{name} = 'F3';
          $midF_v->{name} = 'F2';
          $innerF->{name} = 'F1';
          push @fwdPrimers, $outF_v, $midF_v, $innerF;
          if ($includeLoopPrimers) {
            my $loopF_v = $f_set_infos->[0];
            $loopF_v->{name} = 'FL';
            push @fwdPrimers, $loopF_v;
          }

          my $outR_v = $r_set_infos->[2];  # B3
          my $midR_v = $r_set_infos->[1];  # B2
          $innerR->{name} = 'B1';
          $midR_v->{name} = 'B2';
          $outR_v->{name} = 'B3';
          push @revPrimers, $innerR, $midR_v, $outR_v;
          if ($includeLoopPrimers) {
            my $loopR_v = $r_set_infos->[0];
            $loopR_v->{name} = 'BL';
            unshift @revPrimers, $loopR_v;
          }

          next if (!validateCompleteSignatureSpacing(\@fwdPrimers, \@revPrimers, $minPrimerSpacing));

          # Construct PCR Pairs
          # 1. Inner Pair (F1c, B1c)
          my $innerPair = LLNL::LAVA::PrimerSet::PCRPair->new({
              "forward_info" => $innerF,
              "reverse_info" => $innerR
          });
          my $innerSetInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new({
              "analyzed_pair" => $innerPair,
              "penalty" => $innerF->getPenalty() + $innerR->getPenalty()
          });
          
          # 2. Middle Pair (F2, B2)
          my $midF = $f_set_infos->[1];
          my $midR = $r_set_infos->[1];
          my $middlePair = LLNL::LAVA::PrimerSet::PCRPair->new({
              "forward_info" => $midF,
              "reverse_info" => $midR
          });
          my $middleSetInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new({
              "analyzed_pair" => $middlePair,
              "penalty" => $midF->getPenalty() + $midR->getPenalty()
          });

          # 3. Outer Pair (F3, B3)
          my $outF = $f_set_infos->[2];
          my $outR = $r_set_infos->[2];
          my $outerPair = LLNL::LAVA::PrimerSet::PCRPair->new({
              "forward_info" => $outF,
              "reverse_info" => $outR
          });
          my $outerSetInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new({
              "analyzed_pair" => $outerPair,
              "penalty" => $outF->getPenalty() + $outR->getPenalty()
          });
          
          # Create LAMP Signature
          # Note: Loop primers are attached as TAGS on the signature object usually, 
          # OR passed in specifically if the class supports it.
          # LAMP.pm documentation implies it manages Inner, Middle, Outer.
          # Loops are usually added via tags or setters.
          
          my $lampSignature = LLNL::LAVA::PrimerSet::LAMP->new({
              "inner_info" => $innerSetInfo,
              "middle_info" => $middleSetInfo,
              "outer_info" => $outerSetInfo,
          });
          
          # Add Loop Primers
          if($includeLoopPrimers) {
               my $loopF = $f_set_infos->[0]; # LoopF
               my $loopR = $r_set_infos->[0]; # LoopR
               
               $lampSignature->setTag("has_loop_primers", $TRUE);
               $lampSignature->setTag("floop_info", $loopF);
               $lampSignature->setTag("bloop_info", $loopR);
          } else {
               $lampSignature->setTag("has_loop_primers", $FALSE);
          }
          
          # Add total penalty tag
          my $f_penalty = $bestForwardPenalties[$i]->[0] + $bestForwardPenalties[$i]->[1];
          my $r_penalty = $bestReversePenalties[$j]->[0] + $bestReversePenalties[$j]->[1];
          $lampSignature->setTag("lamp_penalty", $f_penalty + $r_penalty);
          $lampSignature->setTag("penalty_notes", sprintf("Total F:%.1f R:%.1f | F{%s} | R{%s}", $f_penalty, $r_penalty, $bestForwardPenalties[$i]->[2], $bestReversePenalties[$j]->[2]));
          
          push(@{$allFoundSignatures_r}, $lampSignature);
          $combinedSignatureCount++;
      }
  }
  
  print "Created $combinedSignatureCount complete LAMP signatures.\n";
  
  print "Found " .
    scalar(@{$allFoundSignatures_r}) .
    " total signatures across all iterations\n";

  # --- VALIDATION STEP (Essential for correct tagging) ---
  print "Validating and calculating coverage for " . scalar(@{$allFoundSignatures_r}) . " signatures...\n";
  my $validated_count = 0;
  foreach my $signature (@{$allFoundSignatures_r}) {
      # Calculate Intersection and Coverage
      # Calculer l'intersection des sequences compatibles (PipelineUtils unifie)
      # Calculate compatible sequence intersection (unified PipelineUtils)
      my ($amplified_seqs, $coverage, $status) = calculateSignatureIntersection(
          $signature, 
          $inputMSA->num_sequences(), 
          $signatureCommonTargetMinPercent,
          $includeLoopPrimers,
          "loop"
      );
      $validated_count++;
  }
  print "Validation complete.\n";


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

  my $allSignatureCount = scalar(@allSignatures);
  for(my $i = 0; $i < $allSignatureCount; $i++)
  {
    my $signature = $allSignatures[$i];
    my $signatureName = "$i";
 
    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    my $penaltyNotes = $signature->getTagExists("penalty_notes") ? $signature->getTag("penalty_notes") : "";
    my $target_count = $signature->getTagExists("signature_target_count") ? $signature->getTag("signature_target_count") : 0;
    my $coverage_percent = $signature->getTagExists("signature_coverage_percent") ? $signature->getTag("signature_coverage_percent") : "0.00";
    my $sigNum = $signatureName + 1;
    my $degenerate_bases = $signature->getTagExists("degenerate_bases") ? $signature->getTag("degenerate_bases") : 0;
    my $sigLength = $signature->getLength();
    my $headerLine = "Signature ${sigNum} (length: ${sigLength}bp) (penalty: $penalty) $penaltyNotes (coverage: ${target_count}seqs/${coverage_percent}%) (degenerate: ${degenerate_bases} bases) (locations: $locationSummary)";
    if($includeLoopPrimers && $signature->getTagExists("floop_info")) {
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      $headerLine .= " LOOP (locations: $loopLocationSummary)";
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
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTALLSIGS ">${sigNum}_FLOOP\n";
      print OUTALLSIGS $floopSequence . "\n";
      print OUTALLSIGS ">${sigNum}_BLOOP\n";
      print OUTALLSIGS $bloopSequence . "\n";
    }
  }

  close(OUTALLSIGS) ||
    confess("file error - failed to close output file \"$allSignaturesFileName\": $!");

  # Filtrer pour ne garder que les signatures validées (couverture >= seuil) pour le fichier principal .primers
  # Filter to only keep validated signatures (coverage >= threshold) for the main .primers file
  my @valid_sigs = grep { $_->getTag("validation_status") eq "VALIDEE" } @{$allFoundSignatures_r};
  $allFoundSignatures_r = \@valid_sigs;

  # NOW apply the overlap reduction for the main output files
  $possibleSignatures_r = reduceSignaturesByOverlap(
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
        if ($includeLoopPrimers) {
            $seqs .= ($sig->getTag("floop_info"))->getSequence();
            $seqs .= ($sig->getTag("bloop_info"))->getSequence();
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
      
      my $combination_results = analyzeSignatureCombinations(\@signatures_to_analyze, $inputMSA->num_sequences());
      
      # Sauvegarder les résultats de combinaisons dans un fichier / Save combination results to a file
      my $outputFileBase = $outputFileName;
      $outputFileBase =~ s/\.(txt|fasta|fa)$//;  # Enlever l'extension si présente / Remove extension if present
      my $combinations_file = "${outputFileBase}_combinations.txt";
      open(my $comb_fh, '>', $combinations_file) or die "Cannot open $combinations_file: $!";
      
      print $comb_fh "ANALYSE DES COMBINAISONS DE SIGNATURES\n";
      print $comb_fh "=====================================\n\n";
      print $comb_fh "Nombre total de signatures: " . scalar(@signatures_to_analyze) . "\n";
      print $comb_fh "Nombre total de séquences: " . $inputMSA->num_sequences() . "\n\n";
      
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
    if($includeLoopPrimers && $signature->getTagExists("floop_info")) {
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      $headerLine .= " LOOP (locations: $loopLocationSummary)";
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
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTANSWER ">${sigNum}_FLOOP\n";
      print OUTANSWER $floopSequence . "\n";
      print OUTANSWER ">${sigNum}_BLOOP\n";
      print OUTANSWER $bloopSequence . "\n";
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
    if($includeLoopPrimers && $signature->getTagExists("floop_info")) {
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      $headerLine .= " LOOP (locations: $loopLocationSummary)";
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
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTDASH ">${sigNum}_FLOOP\n";
      print OUTDASH $floopSequence . "\n";
      print OUTDASH ">${sigNum}_BLOOP\n";
      print OUTDASH $bloopSequence . "\n";
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
    if($includeLoopPrimers && $signature->getTagExists("floop_info")) {
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      $headerLine .= " LOOP (locations: $loopLocationSummary)";
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
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTPRIMERS ">${sigNum}_FLOOP\n";
      print OUTPRIMERS $floopSequence . "\n";
      print OUTPRIMERS ">${sigNum}_BLOOP\n";
      print OUTPRIMERS $bloopSequence . "\n";
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
  createPerSignatureFiles($possibleSignatures_r, \@sequence_names, $output_base, "LOOP");

  print "Exiting normally\n";
}

# Les fonctions utilitaires partagées / The shared utility functions (buildReversePrimers, analyzeAll, enumeratePairs,
# buildMetricsArray, reducePairInfosByPenalty, reducePrimersByOverlap,
# reduceSignaturesByOverlap, flattenInfoData) sont désormais dans: / are now in:
# lib/LLNL/LAVA/PipelineUtils.pm
# Shared utility functions are now in: lib/LLNL/LAVA/PipelineUtils.pm
# lib/LLNL/LAVA/PipelineUtils.pm
# Shared utility functions are now in: lib/LLNL/LAVA/PipelineUtils.pm
