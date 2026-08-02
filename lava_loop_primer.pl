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
$| = 1;  # Autoflush STDOUT pour l'envoi en temps réel vers Flask
use lib 'lib';

use Getopt::Long;

our $penalty_guard_innerToLoop_neg = 0;
our $penalty_guard_innerToLoop_oob = 0;
our $penalty_guard_loopToMiddle_neg = 0;
our $penalty_guard_loopToMiddle_oob = 0;
our $penalty_guard_innerToMiddle_neg = 0;
our $penalty_guard_innerToMiddle_oob = 0;
our $penalty_guard_middleToOuter_neg = 0;
our $penalty_guard_middleToOuter_oob = 0;

sub penaltyAt {
    my ($table_r, $distance, $label) = @_;
    if ($distance < 0) {
        if ($label eq 'loopToMiddle') { $penalty_guard_loopToMiddle_neg++; }
        elsif ($label eq 'innerToLoop') { $penalty_guard_innerToLoop_neg++; }
        elsif ($label eq 'innerToMiddle') { $penalty_guard_innerToMiddle_neg++; }
        elsif ($label eq 'middleToOuter') { $penalty_guard_middleToOuter_neg++; }
        return 100;
    }
    if (!defined $table_r->[$distance]) {
        if ($label eq 'loopToMiddle') { $penalty_guard_loopToMiddle_oob++; }
        elsif ($label eq 'innerToLoop') { $penalty_guard_innerToLoop_oob++; }
        elsif ($label eq 'innerToMiddle') { $penalty_guard_innerToMiddle_oob++; }
        elsif ($label eq 'middleToOuter') { $penalty_guard_middleToOuter_oob++; }
        return 100;
    }
    return $table_r->[$distance];
}

sub clamp_tm_target {
    my ($target_ref, $min_val, $max_val, $param_prefix) = @_;
    if ($$target_ref < $min_val) {
        printf STDERR "AVERTISSEMENT : --%s_target_tm (%.1f) hors de la plage demandée [%.1f, %.1f]. Le Tm cible est ajusté à %.1f.\n", $param_prefix, $$target_ref, $min_val, $max_val, $min_val;
        $$target_ref = $min_val;
    } elsif ($$target_ref > $max_val) {
        printf STDERR "AVERTISSEMENT : --%s_target_tm (%.1f) hors de la plage demandée [%.1f, %.1f]. Le Tm cible est ajusté à %.1f.\n", $param_prefix, $$target_ref, $min_val, $max_val, $max_val;
        $$target_ref = $max_val;
    }
}


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
use LLNL::LAVA::PipelineUtils qw(getOligosWithMismatchTolerance set_pipeline_threads buildNativeReversePool analyzeAll enumeratePairs buildMetricsArray reducePairInfosByPenalty reducePrimersByOverlap reduceSignaturesByOverlap flattenInfoData buildBigMerge calculateSignatureIntersection createPerSignatureFiles createAmplificationFiles analyzeSignatureCombinations generateCombinations calculateDynamicPairLengths injectFixedPrimers findPrimerPositionInAlignment computeFixedPrimerWindows); # buildReversePrimers retire (DEPRECATED, remplace par buildNativeReversePool)
use LLNL::LAVA::ForkManager;

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


{ # Fake main() to enforce scope
  my %options;
  my %optionMap =
    (
      "alignment_fasta=s" => \$options{"alignment_fasta"},
      "output_file=s" => \$options{"output_file"}, 
      "threads|cpu=s" => \$options{"threads"},
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

      "signature_common_target_min_percent=f" => \$options{"signature_common_target_min_percent"},
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

      # --- NOUVEAUX PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES (AVEC ALIAS HARMONISÉS) ---
      "primer_min_match_percent=f" => \$options{"primer_min_match_percent"},
      "primer_min_iupac_percent|primer_iupac_min_percent=f" => \$options{"primer_min_iupac_percent"},
      "primer_min_coverage_percent|min_primer_coverage=f" => \$options{"primer_min_coverage_percent"},
      # ---------------------------------------------------------------------------------

      # TODO: Not sure if the pair target lengths should be exposed to the 
      # user, or adjusted based on other parameters
      #"outer_pair_target_length=i" => \$options{"outer_pair_target_length"}, 
      #"middle_pair_target_length=i" => \$options{"middle_pair_target_length"}, 
      #"inner_pair_target_length=i" => \$options{"inner_pair_target_length"}, 

      "option_file|options_file=s" => \$options{"option_file"},
      # --- AMORCES FIXEES (peut etre repete plusieurs fois) ---
      # --- FIXED PRIMERS (can be repeated multiple times) ---
      "fixed_primer=s" => \@{$options{"fixed_primer"}},
      "fixed_primer_optimize=i" => \$options{"fixed_primer_optimize"},
    );

  my %optionDefaults =
    (
      "threads" => "auto",
      "fixed_primer" => [],  # Tableau d'amorces fixees / Array of fixed primers
      "fixed_primer_optimize" => 1,
      "signature_max_length" => 400,
      "total_signature_length" => 250,
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

      "signature_common_target_min_percent" => 70,
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
      "primer_min_iupac_percent" => 98,
      "primer_iupac_min_percent" => 98,
      "primer_min_coverage_percent" => 80,
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
      "    [--signature_max_length <length, default=" .
        $optionDefaults{"signature_max_length"} .
	">\n" .
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

      "    [--signature_common_target_min_percent <float, default=" .
        $optionDefaults{"signature_common_target_min_percent"} .
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
      "    [--primer_min_iupac_percent <percent, default=" .
        $optionDefaults{"primer_min_iupac_percent"} .
	"> (alias: --primer_iupac_min_percent)]\n" .
      "    [--primer_min_coverage_percent <percent, default=" .
        $optionDefaults{"primer_min_coverage_percent"} .
	"> (alias: --min_primer_coverage)]\n" .
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
  GetOptions(%optionMap) or die "FATAL: Options de ligne de commande invalides (Unknown option).\n";
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
      $optionDefaults{"total_signature_length"});

  if ($totalSignatureLength > $signatureMaxLength) {
      print "[ATTENTION] total_signature_length ($totalSignatureLength) est superieur a signature_max_length ($signatureMaxLength).\n";
      print "            Recadrage de la cible sur le plafond ($signatureMaxLength).\n";
      $totalSignatureLength = $signatureMaxLength;
  }

  my $maxTotalDegen = optionWithDefault($options_r, "max_total_degenerate_bases", 2);
  my $maxConsecDegen = optionWithDefault($options_r, "max_consecutive_degenerate_bases", 2);
  my $max3PrimeDegen = optionWithDefault($options_r, "max_3prime_degenerate_bases", 2);
  my $maxToleratedMismatches = optionWithDefault($options_r, "max_tolerated_mismatches", 0);
  my $threePrimeZoneSize = optionWithDefault($options_r, "three_prime_zone_size", 6);
  my $minBaseFrequency = optionWithDefault($options_r, "min_base_frequency", 0.05);
  print "Config: Min Base Frequency = $minBaseFrequency\n";
  my $entropyThreshold = optionWithDefault($options_r, "entropy_threshold", 1.5);
  my $maxTmDiff = optionWithDefault($options_r, "max_tm_diff", 5.0);

  # --- Parsing des amorces fixees / Parsing of fixed primers ---
  # Format accepte : TYPE:SEQUENCE  ou  TYPE:SEQUENCE:POSITION
  # Accepted format: TYPE:SEQUENCE  or  TYPE:SEQUENCE:POSITION
  # Types valides LOOP : F3 B3 F2 B2 F1C B1C FLOOP BLOOP
  # Valid types LOOP:   F3 B3 F2 B2 F1C B1C FLOOP BLOOP
  my @fixedPrimerSpecs = ();
  my $fixedPrimerWindows_r = {};  # Fenetre geometrique calculee / Computed geometric window
  my @raw_fixed = @{ $options{"fixed_primer"} // [] };
  for my $raw (@raw_fixed) {
    my @parts = split(/:/, $raw, 3);
    if (@parts < 2 || !$parts[0] || !$parts[1]) {
      print "[FIXED PRIMER] AVERTISSEMENT: Format invalide '$raw'. Attendu TYPE:SEQUENCE ou TYPE:SEQUENCE:POSITION. Ignore.\n";
      print "[FIXED PRIMER] WARNING: Invalid format '$raw'. Expected TYPE:SEQUENCE or TYPE:SEQUENCE:POSITION. Skipped.\n";
      next;
    }
    my $type = uc($parts[0]);
    $type = "F1C" if $type eq "F1";
    $type = "B1C" if $type eq "B1";
    
    my $spec = {
      type => $type,
      seq  => uc($parts[1]),
      pos  => (defined $parts[2] && $parts[2] =~ /^\d+$/) ? int($parts[2]) : undef,
    };
    # Validation du type pour LOOP
    my %valid_loop_types = map { $_ => 1 } qw(F3 B3 F2 B2 F1C B1C FLOOP BLOOP);
    if (!$valid_loop_types{$spec->{type}}) {
      print STDERR "ERROR: [FIXED PRIMER] Type '$spec->{type}' non reconnu pour LOOP. Types valides: F3 B3 F2 B2 F1C B1C FLOOP BLOOP (F1 et B1 sont acceptes comme alias de F1C et B1C).\n";
      exit(2);
    }
    push @fixedPrimerSpecs, $spec;
    printf("[FIXED PRIMER] Spec enregistree: TYPE=%s SEQ=%s POS=%s\n",
           $spec->{type}, $spec->{seq}, defined $spec->{pos} ? $spec->{pos} : "auto");
  }
  my %isFixedType = ();
  for my $spec (@fixedPrimerSpecs) {
    $isFixedType{$spec->{type}} = 1;
  }



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
  my $outerPrimerMaxTM =
    optionWithDefault($options_r, "outer_primer_max_tm", 
      ($outerPrimerTargetTM + 1.0));
  clamp_tm_target(\$outerPrimerTargetTM, $outerPrimerMinTM, $outerPrimerMaxTM, "outer_primer");

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
  my $loopPrimerMaxTM =
    optionWithDefault($options_r, "loop_primer_max_tm", 
      ($loopPrimerTargetTM + 1.0));
  clamp_tm_target(\$loopPrimerTargetTM, $loopPrimerMinTM, $loopPrimerMaxTM, "loop_primer");


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
  my $middlePrimerMaxTM =
    optionWithDefault($options_r, "middle_primer_max_tm", 
      ($middlePrimerTargetTM + 1.0));
  clamp_tm_target(\$middlePrimerTargetTM, $middlePrimerMinTM, $middlePrimerMaxTM, "middle_primer");

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
  my $innerPrimerMaxTM =
    optionWithDefault($options_r, "inner_primer_max_tm", 
      ($innerPrimerTargetTM + 1.0));
  clamp_tm_target(\$innerPrimerTargetTM, $innerPrimerMinTM, $innerPrimerMaxTM, "inner_primer");

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
    optionWithDefault($options_r, "signature_common_target_min_percent",
      $optionDefaults{"signature_common_target_min_percent"});
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

  # --- RÉCUPÉRATION DES PARAMÈTRES DE TOLÉRANCE AUX MISMATCHES (AVEC ALIAS HARMONISÉS) ---
  $options_r->{"primer_min_iupac_percent"} //= $options_r->{"primer_iupac_min_percent"};
  $options_r->{"primer_iupac_min_percent"} //= $options_r->{"primer_min_iupac_percent"};
  $options_r->{"primer_min_coverage_percent"} //= $options_r->{"min_primer_coverage"};
  $options_r->{"min_primer_coverage"} //= $options_r->{"primer_min_coverage_percent"};

  my $primerMinMatchPercent = optionWithDefault($options_r, "primer_min_match_percent",
    $optionDefaults{"primer_min_match_percent"});
  my $primerIupacMinPercent = optionWithDefault($options_r, "primer_min_iupac_percent", 
    $optionDefaults{"primer_min_iupac_percent"});
  my $minPrimerCoverage = optionWithDefault($options_r, "primer_min_coverage_percent", 
    $optionDefaults{"primer_min_coverage_percent"});
  
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

  set_pipeline_threads($options_r->{"threads"});

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

  # --- Resolution anticipee des positions des amorces fixees ---
  # --- Early position resolution for fixed primers ---
  # On resout les positions maintenant (avant Primer3) pour pouvoir calculer
  # la fenetre geometrique et filtrer les candidats AVANT leur generation.
  # We resolve positions now (before Primer3) to compute the geometric window
  # and filter candidates BEFORE their generation.
  if (@fixedPrimerSpecs) {
    print "\n[FIXED PRIMER WINDOW] Resolution anticipee des positions / Early position resolution...\n";
    for my $spec (@fixedPrimerSpecs) {
      next if defined $spec->{pos};  # position deja fournie / already provided
      my ($found_pos, $found_strand) = findPrimerPositionInAlignment($inputMSA, $spec->{seq}, undef);
      if (defined $found_pos) {
        $spec->{pos}    = $found_pos;
        $spec->{strand} = $found_strand;
        printf("[FIXED PRIMER WINDOW] %s '%s' -> position resolue : %d (brin %s)\n",
               $spec->{type}, $spec->{seq}, $found_pos, $found_strand);
      } else {
        print "[FIXED PRIMER WINDOW] AVERTISSEMENT: position introuvable pour $spec->{type} '$spec->{seq}'. Fenetre non contrainte pour cette amorce.\n";
        print "[FIXED PRIMER WINDOW] WARNING: position not found for $spec->{type} '$spec->{seq}'. No window constraint for this primer.\n";
      }
    }

    my $fixed_window_margin = optionWithDefault($options_r, "fixed_primer_margin", 200);
    my $current_sig_max = optionWithDefault($options_r, "signature_max_length", 320);
    $fixedPrimerWindows_r = computeFixedPrimerWindows(
      \@fixedPrimerSpecs, $current_sig_max, $fixed_window_margin, $sequenceLength
    );
  }


  my $global_included_region = undef;
  if (%{$fixedPrimerWindows_r}) {
    my $win = $fixedPrimerWindows_r->{"F3"};
    if (defined $win) {
      my $len = $win->[1] - $win->[0] + 1;
      if ($len > 0) {
        $global_included_region = $win->[0] . "," . $len;
        print "[FIXED PRIMER WINDOW] Passing SEQUENCE_INCLUDED_REGION=$global_included_region to Primer3\n";
      }
    }
  }

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
      (defined $global_included_region ? ("included_region" => $global_included_region) : ()),
    });

  print "Enumerating outer forward primers\n";
  my @outerForwardPrimers = ();
  if ($isFixedType{"F3"}) {
    print "Skipping outer forward (F3) enumeration because it is fixed.\n";
  } else {
    @outerForwardPrimers = getOligosWithMismatchTolerance($outerEnumerator, $inputMSA,
                                                          $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Outer Forward (F3)");

    print "  Generated \"" . scalar(@outerForwardPrimers) . "\" outer forward primers\n";
  }
  
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
  my @outerReversePrimers = ();
  if ($isFixedType{"B3"}) {
    print "Skipping outer reverse (B3) enumeration because it is fixed.\n";
  } else {
    @outerReversePrimers = buildNativeReversePool(
      $outerEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Outer Reverse (B3)"
    );
    print "  Generated \"" . scalar(@outerReversePrimers) . "\" outer native reverse primers\n";
  }


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
      (defined $global_included_region ? ("included_region" => $global_included_region) : ()),
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
  if ($isFixedType{"BLOOP"}) {
    print "Skipping loop BACK (BLOOP) enumeration because it is fixed.\n";
  } else {
    @loopBackPrimers = getOligosWithMismatchTolerance($loopEnumerator, $inputMSA,
                                                        $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                        $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Loop Back (BLOOP)");

    print "  Generated \"" . scalar(@loopBackPrimers) . "\" loop BACK (BLOOP) primers\n";
  }

  # FLOOP : Option B - genere nativement sur RC(MSA) pour garantir la protection 3'
  # FLOOP: Option B - natively generated on RC(MSA) to guarantee 3-prime protection
  # (Forward Loop = antisens, 3' pointe vers F1c - correspondait avant au 5' du BLOOP source = bug)
  print "Enumerating loop FORWARD (FLOOP) NATIVE reverse primers (Option B)\n";
  if ($isFixedType{"FLOOP"}) {
    print "Skipping loop FORWARD (FLOOP) enumeration because it is fixed.\n";
  } else {
    @loopForwardPrimers = buildNativeReversePool(
      $loopEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Loop Forward (FLOOP)"
    );
    print "  Generated \"" . scalar(@loopForwardPrimers) . "\" loop FORWARD (FLOOP) native primers\n";
  }
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
      (defined $global_included_region ? ("included_region" => $global_included_region) : ()),
    });

  print "Enumerating middle forward primers\n";
  my @middleForwardPrimers = ();
  if ($isFixedType{"F2"}) {
    print "Skipping middle forward (F2) enumeration because it is fixed.\n";
  } else {
    @middleForwardPrimers = getOligosWithMismatchTolerance($middleEnumerator, $inputMSA,
                                                           $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                           $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Middle Forward (F2)");

    print "  Generated \"" . scalar(@middleForwardPrimers) . "\" middle primers\n";
  }

  # Option B : Generation NATIVE des Reverse Middle via Primer3 sur RC(MSA)
  print "Enumerating middle NATIVE reverse primers (Option B)\n";
  my @middleReversePrimers = ();
  if ($isFixedType{"B2"}) {
    print "Skipping middle reverse (B2) enumeration because it is fixed.\n";
  } else {
    @middleReversePrimers = buildNativeReversePool(
      $middleEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Middle Reverse (B2)"
    );
    print "  Generated \"" . scalar(@middleReversePrimers) . "\" middle native reverse primers\n";
  }

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
      (defined $global_included_region ? ("included_region" => $global_included_region) : ()),
    });

  print "Enumerating inner forward primers\n";
  my @innerForwardPrimers = ();
  if ($isFixedType{"F1C"}) {
    print "Skipping inner forward (F1) enumeration because it is fixed.\n";
  } else {
    @innerForwardPrimers = getOligosWithMismatchTolerance($innerEnumerator, $inputMSA,
                                                          $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
                                                          $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency, "Inner Forward (F1)");

    print "  Generated \"" . scalar(@innerForwardPrimers) . "\" inner primers\n";
  }
  
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
  my @innerReversePrimers = ();
  if ($isFixedType{"B1C"}) {
    print "Skipping inner reverse (B1) enumeration because it is fixed.\n";
  } else {
    @innerReversePrimers = buildNativeReversePool(
      $innerEnumerator, $inputMSA,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen, $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      \&checkPrimerMismatchTolerance, \&isIUPACCompatible, \&rev_comp, "Inner Reverse (B1)"
    );
    print "  Generated \"" . scalar(@innerReversePrimers) . "\" inner native reverse primers\n";
  }

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
  # --- FILTRAGE GEOMETRIQUE PAR FENETRE (si amorces fixees) ---
  # --- GEOMETRIC WINDOW FILTERING (if fixed primers are defined) ---
  # Elimination des candidats geometriquement impossibles AVANT l'injection.
  # Eliminating geometrically impossible candidates BEFORE injection.
  # Note: la penalite sigmoide reste inchangee pour les candidats dans la marge.
  # Note: sigmoid penalty remains unchanged for candidates within the margin.
  if (%{$fixedPrimerWindows_r}) {
    my $win_F3    = $fixedPrimerWindows_r->{"F3"}    // [0, 999_999];
    my $win_B3    = $fixedPrimerWindows_r->{"B3"}    // [0, 999_999];
    my $win_F2    = $fixedPrimerWindows_r->{"F2"}    // [0, 999_999];
    my $win_B2    = $fixedPrimerWindows_r->{"B2"}    // [0, 999_999];
    my $win_F1C   = $fixedPrimerWindows_r->{"F1C"}   // [0, 999_999];
    my $win_B1C   = $fixedPrimerWindows_r->{"B1C"}   // [0, 999_999];
    my $win_FLOOP = $fixedPrimerWindows_r->{"FLOOP"} // [0, 999_999];
    my $win_BLOOP = $fixedPrimerWindows_r->{"BLOOP"} // [0, 999_999];

    my $before_fwd_outer  = scalar(@outerForwardPrimers);
    my $before_rev_outer  = scalar(@outerReversePrimers);
    my $before_fwd_middle = scalar(@middleForwardPrimers);
    my $before_rev_middle = scalar(@middleReversePrimers);
    my $before_fwd_inner  = scalar(@innerForwardPrimers);
    my $before_rev_inner  = scalar(@innerReversePrimers);
    my $before_floop      = scalar(@loopForwardPrimers);
    my $before_bloop      = scalar(@loopBackPrimers);

    @outerForwardPrimers  = grep { $_->location() >= $win_F3->[0]    && $_->location() <= $win_F3->[1]    } @outerForwardPrimers;
    @outerReversePrimers  = grep { $_->location() >= $win_B3->[0]    && $_->location() <= $win_B3->[1]    } @outerReversePrimers;
    @middleForwardPrimers = grep { $_->location() >= $win_F2->[0]    && $_->location() <= $win_F2->[1]    } @middleForwardPrimers;
    @middleReversePrimers = grep { $_->location() >= $win_B2->[0]    && $_->location() <= $win_B2->[1]    } @middleReversePrimers;
    @innerForwardPrimers  = grep { $_->location() >= $win_F1C->[0]   && $_->location() <= $win_F1C->[1]   } @innerForwardPrimers;
    @innerReversePrimers  = grep { $_->location() >= $win_B1C->[0]   && $_->location() <= $win_B1C->[1]   } @innerReversePrimers;
    @loopForwardPrimers   = grep { $_->location() >= $win_FLOOP->[0] && $_->location() <= $win_FLOOP->[1] } @loopForwardPrimers;
    @loopBackPrimers      = grep { $_->location() >= $win_BLOOP->[0] && $_->location() <= $win_BLOOP->[1] } @loopBackPrimers;

    printf("[FIXED PRIMER WINDOW] Filtrage LOOP terminé / LOOP filtering done:\n");
    printf("  F3:    %d -> %d | B3:    %d -> %d\n", $before_fwd_outer,  scalar(@outerForwardPrimers),
                                                     $before_rev_outer,  scalar(@outerReversePrimers));
    printf("  F2:    %d -> %d | B2:    %d -> %d\n", $before_fwd_middle, scalar(@middleForwardPrimers),
                                                     $before_rev_middle, scalar(@middleReversePrimers));
    printf("  F1:   %d -> %d | B1:   %d -> %d\n", $before_fwd_inner,  scalar(@innerForwardPrimers),
                                                     $before_rev_inner,  scalar(@innerReversePrimers));
    printf("  FLOOP: %d -> %d | BLOOP: %d -> %d\n", $before_floop, scalar(@loopForwardPrimers),
                                                     $before_bloop, scalar(@loopBackPrimers));
  }

  # --- INJECTION DES AMORCES FIXEES dans les pools correspondants ---

  # --- INJECT FIXED PRIMERS into the corresponding pools ---
  if (@fixedPrimerSpecs) {
    print "\n=== Injection des amorces fixees / Fixed Primer Injection ===\n";
    my %target_tms = (
      "F3" => $outerPrimerTargetTM, "B3" => $outerPrimerTargetTM,
      "F2" => $middlePrimerTargetTM, "B2" => $middlePrimerTargetTM,
      "F1C" => $innerPrimerTargetTM, "B1C" => $innerPrimerTargetTM,
      "FLOOP" => $loopPrimerTargetTM, "BLOOP" => $loopPrimerTargetTM
    );
    
    my $fixed_results_r = injectFixedPrimers(
      $inputMSA, \@fixedPrimerSpecs,
      $primerMinMatchPercent, $primerIupacMinPercent, $minPrimerCoverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen,
      $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      $options{"fixed_primer_optimize"},
      \%target_tms,
      $options_r
    );
    # Fusionner les amorces fixees dans chaque pool / Merge fixed primers into each pool
    unshift @outerForwardPrimers,  @{ $fixed_results_r->{"F3"}    // [] };
    unshift @outerReversePrimers,  @{ $fixed_results_r->{"B3"}    // [] };
    unshift @middleForwardPrimers, @{ $fixed_results_r->{"F2"}    // [] };
    unshift @middleReversePrimers, @{ $fixed_results_r->{"B2"}    // [] };
    unshift @innerForwardPrimers,  @{ $fixed_results_r->{"F1C"}   // [] };
    unshift @innerReversePrimers,  @{ $fixed_results_r->{"B1C"}   // [] };
    unshift @loopForwardPrimers,   @{ $fixed_results_r->{"FLOOP"} // [] };
    unshift @loopBackPrimers,      @{ $fixed_results_r->{"BLOOP"} // [] };
    print "=== Injection terminee / Fixed Primer Injection done ===\n\n";
  }


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
    print "Analyse / Analysis of loop primers ignorée\n";
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
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$outerReversePrimerMeasurements_r};

  my @outerForwardInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
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
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$loopForwardPrimerMeasurements_r};
    @loopBackInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$loopBackPrimerMeasurements_r};

    @loopForwardInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$loopForwardPrimerMeasurements_r};
    @loopBackInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$loopBackPrimerMeasurements_r};
  }

  # Middle primers sorted 2 ways
  my @middleForwardInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$middleReversePrimerMeasurements_r};

  my @middleForwardInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$middleReversePrimerMeasurements_r};

  # Inner primers sorted 2 ways
  my @innerForwardInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByLocation =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getLocation()] } 
    @{$innerReversePrimerMeasurements_r};

  my @innerForwardInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
    map {[$_, $_->getPenalty()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByPenalty =
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getLocation() <=> $b->[0]->getLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getSequence() cmp $b->[0]->getSequence() }
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
  
  print "INFO: Cibles Géométriques Proportionnelles (Cible = $totalSignatureLength pb) :\n";
  print "  -> F3-F2 (12%) : " . $geometry->{'f3_f2_target'} . " nt\n";
  print "  -> F2-F1 (18%) : " . $geometry->{'f2_f1_target'} . " nt\n";
  print "  -> Empan interne F1-B1 (40%) : " . $geometry->{'inner_target'} . " nt\n";
  
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
  
  # --- B&B Initialization Forward ---
  my $min_val_f = sub { my $m = $_[0]; for(@_) { $m = $_ if $_ < $m } return $m; };
  my $minS_innerToLoop_F = @$innerToLoopPenalties_r ? $min_val_f->(@$innerToLoopPenalties_r) * $innerToLoopPenaltyWeight : 0;
  my $minS_loopToMiddle_F = @$loopToMiddlePenalties_r ? $min_val_f->(@$loopToMiddlePenalties_r) * $loopToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_F = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $minS_innerToMiddle_F = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $rmq_middle_f = build_rmq($masterMiddleF_data_r, 2);
  my $rmq_outer_f  = build_rmq($masterOuterF_data_r, 2);
  my $rmq_loop_f   = build_rmq($masterLoopF_data_r, 2) if $includeLoopPrimers;
  my $min_P_outer_F = @$masterOuterF_data_r ? query_rmq($rmq_outer_f, 0, scalar(@$masterOuterF_data_r)-1) * $outerPenaltyWeight : 0;
  
  my $_sig_fwd_pruned = 0;
  my $_sig_fwd_evaluated = 0;
  # --- Counters for zero-signature diagnostic ---
  my $_fwd_rej_geometry = 0;
  my $_fwd_rej_spacing = 0;
  my $_fwd_rej_loopgap = 0;
  my $_fwd_rej_tm_inner_loop = 0;
  my $_fwd_rej_tm_loop_middle = 0;
  my $_fwd_rej_tm_inner_middle = 0;
  my $_fwd_rej_tm_middle_outer = 0;
  my $_fwd_min_delta_tm_inner_loop = 999;
  my $_fwd_min_delta_tm_loop_middle = 999;
  my $_fwd_min_delta_tm_inner_middle = 999;
  my $_fwd_min_delta_tm_middle_outer = 999;
  my $_fwd_min_span_needed = 999999;
  
  my %_fwd_pen_guards;
  my $fwd_prog_dir = "$options_r->{'output_file'}_fwd_prog_$$";
  $fwd_prog_dir = "$options_r->{'output_file'}_fwd_prog_$$" if ref($options_r);
  use File::Path qw(make_path remove_tree);
  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;
  make_path($fwd_prog_dir);
  my $_sig_fwd_t0   = time();
  my $_sig_fwd_done = 0;
  my $_sig_fwd_hits = 0;  # Nombre de signatures Forward trouvees / Forward signatures found
  my $pm_fwd = LLNL::LAVA::ForkManager->new($options_r->{"threads"});
  my $num_fwd_chunks = $pm_fwd->{max_processes} * 12;
  $num_fwd_chunks = 30 if $num_fwd_chunks < 30;
  $num_fwd_chunks = $innerForwardCount if $num_fwd_chunks > $innerForwardCount;
  $num_fwd_chunks = 1 if $num_fwd_chunks < 1;
  my $fwd_chunk_size = int(($innerForwardCount + $num_fwd_chunks - 1) / $num_fwd_chunks);
  $fwd_chunk_size = 1 if $fwd_chunk_size < 1;

  $pm_fwd->run_on_finish(sub {
      my ($pid, $exit_code, $id, $exit_signal, $core_dump, $data_ref) = @_;
      if (defined $data_ref && ref($data_ref) eq 'HASH') {
          foreach my $idx (sort { $a <=> $b } keys %{$data_ref->{infos}}) {
              if (!defined $bestForwardInfos[$idx]) {
                  $forwardSetCount++;
              }
              $bestForwardInfos[$idx] = $data_ref->{infos}->{$idx};
              $bestForwardPenalties[$idx] = $data_ref->{penalties}->{$idx};
          }
          $_sig_fwd_hits += $data_ref->{hits} || 0;
          $_sig_fwd_done += $data_ref->{done} || 0;
          $_sig_fwd_pruned += $data_ref->{pruned} || 0;
          $_sig_fwd_evaluated += $data_ref->{evaluated} || 0;
          
          $_fwd_rej_geometry += $data_ref->{rej_geometry} || 0;
          $_fwd_rej_spacing += $data_ref->{rej_spacing} || 0;
          $_fwd_rej_loopgap += $data_ref->{rej_loopgap} || 0;
          $_fwd_rej_tm_inner_loop += $data_ref->{rej_tm_inner_loop} || 0;
          $_fwd_rej_tm_loop_middle += $data_ref->{rej_tm_loop_middle} || 0;
          $_fwd_rej_tm_inner_middle += $data_ref->{rej_tm_inner_middle} || 0;
          $_fwd_rej_tm_middle_outer += $data_ref->{rej_tm_middle_outer} || 0;
          
          foreach my $k (qw(min_tm_inner_loop min_tm_loop_middle min_tm_inner_middle min_tm_middle_outer min_span_needed)) {
              next unless defined $data_ref->{$k};
              if ($k eq 'min_tm_inner_loop') {
                  $_fwd_min_delta_tm_inner_loop = $data_ref->{$k} if $data_ref->{$k} < $_fwd_min_delta_tm_inner_loop;
              } elsif ($k eq 'min_tm_loop_middle') {
                  $_fwd_min_delta_tm_loop_middle = $data_ref->{$k} if $data_ref->{$k} < $_fwd_min_delta_tm_loop_middle;
              } elsif ($k eq 'min_tm_inner_middle') {
                  $_fwd_min_delta_tm_inner_middle = $data_ref->{$k} if $data_ref->{$k} < $_fwd_min_delta_tm_inner_middle;
              } elsif ($k eq 'min_tm_middle_outer') {
                  $_fwd_min_delta_tm_middle_outer = $data_ref->{$k} if $data_ref->{$k} < $_fwd_min_delta_tm_middle_outer;
              } elsif ($k eq 'min_span_needed') {
                  $_fwd_min_span_needed = $data_ref->{$k} if $data_ref->{$k} < $_fwd_min_span_needed;
              }
          }
          foreach my $k (qw(innerToLoop loopToMiddle innerToMiddle middleToOuter)) {
              $_fwd_pen_guards{"${k}_neg"} += $data_ref->{pen_guards}->{"${k}_neg"} || 0;
              $_fwd_pen_guards{"${k}_oob"} += $data_ref->{pen_guards}->{"${k}_oob"} || 0;
          }
      }
  });

  for (my $chunk_id = 0; $chunk_id < $num_fwd_chunks; $chunk_id++) {
      $pm_fwd->start($chunk_id) and next;
      
      my %chunk_infos = ();
      my %chunk_penalties = ();
      my $chunk_hits = 0;
      my $chunk_done = 0;
      my $chunk_pruned = 0;
      my $chunk_evaluated = 0;
      my $chunk_rej_geometry = 0;
      my $chunk_rej_spacing = 0;
      my $chunk_rej_loopgap = 0;
      my $chunk_rej_tm_inner_loop = 0;
      my $chunk_rej_tm_loop_middle = 0;
      my $chunk_rej_tm_inner_middle = 0;
      my $chunk_rej_tm_middle_outer = 0;
      my $chunk_min_delta_tm_inner_loop = 999;
      my $chunk_min_delta_tm_loop_middle = 999;
      my $chunk_min_delta_tm_inner_middle = 999;
      my $chunk_min_delta_tm_middle_outer = 999;
      my $chunk_min_span_needed = 999999;

      
      for(my $innerIndex = $chunk_id; $innerIndex < $innerForwardCount; $innerIndex += $num_fwd_chunks)
      {
          $chunk_done++;
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
                  if ($loopLocation < $searchStartAt) { $chunk_rej_geometry++; next; }
                  last if ($loopLocation > $loopEndAt);

                  # --- DYNAMIC THERMAL FILTER (Inner vs Loop) ---
                  my $diff = abs($innerTm - $loopTm);
                  if ($diff > $maxTmDiff && !($innerInfo->hasTag("is_fixed") || $loopInfo->hasTag("is_fixed"))) {
                      $chunk_min_delta_tm_inner_loop = $diff if $diff < $chunk_min_delta_tm_inner_loop;
                      $chunk_rej_tm_inner_loop++;
                      next;
                  }
              }
              
              # 3.2 Calculate Search Bounds for Middle Primer (F2) (Structure: F3-F2-LF-F1c)
              my $middleStartAt = $searchStartAt;
              my $middleEndAt = $loopLocation - $loopLength - $minPrimerSpacing + $middlePrimerMaxLength - 1;
              
              # Ensure Middle isn't too close to Inner (respecting loop gap or standard gap)
              # [restauré depuis a6098f5 : le refactor 55328b2 avait supprime cette contrainte]
              if ($includeLoopPrimers) {
                  my $altMiddleEndAt = $innerLocation - ($loopMinGap + 1);
                  $middleEndAt = $altMiddleEndAt if ($altMiddleEndAt < $middleEndAt);
              } else {
                  my $altMiddleEndAt = $innerLocation - $minPrimerSpacing;
                  $middleEndAt = $altMiddleEndAt if ($altMiddleEndAt < $middleEndAt);
              }
              $middleEndAt = 0 if $middleEndAt < 0;
              
              my $innerToLoopDistance = $innerLocation - ($loopLocation + 1);
              
              
              my $m_start = binary_search_first_ge($masterMiddleF_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleF_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_f, $m_start, $m_end) * $middlePenaltyWeight;
                  my $innerToLoopDistance = $includeLoopPrimers ? $innerLocation - ($loopLocation + 1) : 0;
                  
                  my $pen_innerToLoop = penaltyAt($innerToLoopPenalties_r, $innerToLoopDistance, 'innerToLoop');
                  if ($pen_innerToLoop < 0) { $penalty_guard_innerToLoop_neg++; }
                  if ($pen_innerToLoop == 9999) { $penalty_guard_innerToLoop_oob++; }
                  
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + 
                                     ($includeLoopPrimers ? $pen_innerToLoop * $innerToLoopPenaltyWeight : 0);
                  my $min_S_to_mid = $includeLoopPrimers ? $minS_loopToMiddle_F : $minS_innerToMiddle_F;
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_F + $min_S_to_mid + $minS_middleToOuter_F >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleF_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleF_data_r->[$j]};
                          
                          if ($includeLoopPrimers) {
                              my $loopToMiddleDistance = ($loopLocation - $loopLength + 1) - ($middleLocation + $middleLength);
                              if ($loopToMiddleDistance < 0) {
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $needed = $minPrimerSpacing - $loopToMiddleDistance;
                              if ($needed > 0) {
                                  my $span = ($innerLocation + $innerLength) - $middleLocation + $needed;
                                  $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              if ($middleLocation + $middleLength + $loopMinGap > $innerLocation) { $chunk_rej_loopgap++; next; }
                              my $diffLM = abs($loopTm - $midTm);
                              if ($diffLM > $maxTmDiff && !($loopInfo->hasTag("is_fixed") || $middleInfo->hasTag("is_fixed"))) {
                                  $chunk_min_delta_tm_loop_middle = $diffLM if $diffLM < $chunk_min_delta_tm_loop_middle;
                                  $chunk_rej_tm_loop_middle++;
                                  next;
                              }
                          } else {
                              my $innerToMiddleDistance = $innerLocation - ($middleLocation + $middleLength);
                              if ($innerToMiddleDistance < 0) {
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $needed = $minPrimerSpacing - $innerToMiddleDistance;
                              if ($needed > 0) {
                                  my $span = ($innerLocation + $innerLength) - $middleLocation + $needed;
                                  $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $diffIM = abs($innerTm - $midTm);
                              if ($diffIM > $maxTmDiff && !($innerInfo->hasTag("is_fixed") || $middleInfo->hasTag("is_fixed"))) {
                                  $chunk_min_delta_tm_inner_middle = $diffIM if $diffIM < $chunk_min_delta_tm_inner_middle;
                                  $chunk_rej_tm_inner_middle++;
                                  next;
                              }
                          }
                          
                          my $outerStartAt = $searchStartAt;
                          my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing + $outerPrimerMaxLength;
                          my $loopToMiddleDistance = $includeLoopPrimers ? ($loopLocation - $loopLength + 1) - ($middleLocation + $middleLength) : 0;
                          my $innerToMiddleDistance = $innerLocation - ($middleLocation + $middleLength);
                          
                          my $o_start = binary_search_first_ge($masterOuterF_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterF_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_f, $o_start, $o_end) * $outerPenaltyWeight;
                              
                              my $pen_mid = $includeLoopPrimers ? penaltyAt($loopToMiddlePenalties_r, $loopToMiddleDistance, 'loopToMiddle') : penaltyAt($innerToMiddlePenalties_r, $innerToMiddleDistance, 'innerToMiddle');
                              if ($pen_mid < 0) { if($includeLoopPrimers) { $penalty_guard_loopToMiddle_neg++; } else { $penalty_guard_innerToMiddle_neg++; } }
                              if ($pen_mid == 9999) { if($includeLoopPrimers) { $penalty_guard_loopToMiddle_oob++; } else { $penalty_guard_innerToMiddle_oob++; } }

                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($pen_mid * ($includeLoopPrimers ? $loopToMiddlePenaltyWeight : $innerToMiddlePenaltyWeight));
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_F >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterF_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterF_data_r->[$k]};
                                      
                                      my $middleToOuterDistance = $middleLocation - ($outerLocation + $outerLength);
                                      if ($middleToOuterDistance < 0) {
                                          $chunk_rej_spacing++;
                                          next;
                                      }
                                      my $needed = $minPrimerSpacing - $middleToOuterDistance;
                                      if ($needed > 0) {
                                          my $span = ($innerLocation + $innerLength) - $outerLocation + $needed;
                                          $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                          $chunk_rej_spacing++;
                                          next;
                                      }
                                      my $diffMO = abs($midTm - $outTm);
                                      if ($diffMO > $maxTmDiff && !($middleInfo->hasTag("is_fixed") || $outerInfo->hasTag("is_fixed"))) {
                                          $chunk_min_delta_tm_middle_outer = $diffMO if $diffMO < $chunk_min_delta_tm_middle_outer;
                                          $chunk_rej_tm_middle_outer++;
                                          next;
                                      }
                                      
                                      my $pen_out = penaltyAt($middleToOuterPenalties_r, $middleToOuterDistance, 'middleToOuter');
                                      if ($pen_out < 0) { $penalty_guard_middleToOuter_neg++; }
                                      if ($pen_out == 9999) { $penalty_guard_middleToOuter_oob++; }
                                      
                                      my $spacingPenalty = ($includeLoopPrimers ? $pen_innerToLoop * $innerToLoopPenaltyWeight : 0) +
                                                           ($pen_mid * ($includeLoopPrimers ? $loopToMiddlePenaltyWeight : $innerToMiddlePenaltyWeight)) +
                                                           ($pen_out * $middleToOuterPenaltyWeight);
                                      
                                      my $primer3Penalty = $innerPenalty * $innerPenaltyWeight + ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                      
                                      my $detailStr = "";
                                      if ($includeLoopPrimers) {
                                          $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                                                ($pen_innerToLoop * $innerToLoopPenaltyWeight), ($pen_mid * $loopToMiddlePenaltyWeight), ($pen_out * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($loopPenalty * $loopPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($pen_mid * $innerToMiddlePenaltyWeight), ($pen_out * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      }
                                      
                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $bestSetPenalty) {
                                          $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                                          $chunk_infos{$innerIndex} = $includeLoopPrimers ? [$loopInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo];
                                          $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                                          $bestSetPenalty = $currentSetPenalty;
                                      }
                                  }
                              }
                          }
                      }
                  }
              }

          } # End Loop
          
          # Intra-chunk progress reporting
          if ($chunk_done % 5 == 0) {
              my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
              if (open(my $fh, '>', $prog_file_me)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n";
                  close($fh);
              }
                  
              my $total_done = 0;
              my $total_hits = 0;
              foreach my $f (glob("$fwd_prog_dir/chunk_*.prog")) {
                  if (open(my $r, '<', $f)) {
                      my $line = <$r>; close($r);
                      next unless defined $line;
                      chomp $line;
                      my ($d, $h) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                  }
              }
                      
              if ($_LAVA_IS_TTY || 1) {
                  my $elapsed = time() - $_sig_fwd_t0 + 0.001;
                  my $eta = ($total_done < $innerForwardCount) ? int(($innerForwardCount - $total_done) / ($total_done / $elapsed)) : 0;
                  my $rate = $total_done / $elapsed;
                  printf("[LAVA-PROGRESS] Signatures Forward|%d|%d|Sig: %d|%.1f it/s|%d\r", $total_done, $innerForwardCount, $total_hits, $rate, $eta);
                  my $old_h = select(STDOUT); $| = 1; select($old_h);
              }
          }
} # End Inner chunk loop
      
      my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
      if (open(my $fh, '>', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n"; close($fh); }
      
      $pm_fwd->finish(0, {
          infos => \%chunk_infos,
          penalties => \%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
          rej_geometry => $chunk_rej_geometry,
          rej_spacing => $chunk_rej_spacing,
          rej_loopgap => $chunk_rej_loopgap,
          rej_tm_inner_loop => $chunk_rej_tm_inner_loop,
          rej_tm_loop_middle => $chunk_rej_tm_loop_middle,
          rej_tm_inner_middle => $chunk_rej_tm_inner_middle,
          rej_tm_middle_outer => $chunk_rej_tm_middle_outer,
          min_tm_inner_loop => $chunk_min_delta_tm_inner_loop,
          min_tm_loop_middle => $chunk_min_delta_tm_loop_middle,
          min_tm_inner_middle => $chunk_min_delta_tm_inner_middle,
          min_tm_middle_outer => $chunk_min_delta_tm_middle_outer,
          min_span_needed => $chunk_min_span_needed,
          pen_guards => {
              innerToLoop_neg => $penalty_guard_innerToLoop_neg,
              innerToLoop_oob => $penalty_guard_innerToLoop_oob,
              loopToMiddle_neg => $penalty_guard_loopToMiddle_neg,
              loopToMiddle_oob => $penalty_guard_loopToMiddle_oob,
              innerToMiddle_neg => $penalty_guard_innerToMiddle_neg,
              innerToMiddle_oob => $penalty_guard_innerToMiddle_oob,
              middleToOuter_neg => $penalty_guard_middleToOuter_neg,
              middleToOuter_oob => $penalty_guard_middleToOuter_oob,
          }
      });
  } # End chunks
  $pm_fwd->wait_all_children();
  use File::Path qw(remove_tree);
  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;
  
  if ($_sig_fwd_pruned + $_sig_fwd_evaluated > 0) {
      my $pct = ($_sig_fwd_pruned / ($_sig_fwd_pruned + $_sig_fwd_evaluated)) * 100;
      print sprintf("    [Forward B&B] Elagage: %.2f%%\n", $pct);
  }

  my @guard_msgs = ();
  foreach my $k (qw(innerToLoop loopToMiddle innerToMiddle middleToOuter)) {
      if ($_fwd_pen_guards{"${k}_neg"}) {
          push @guard_msgs, "$k(neg: $_fwd_pen_guards{\"${k}_neg\"})";
      }
      if ($_fwd_pen_guards{"${k}_oob"}) {
          push @guard_msgs, "$k(oob: $_fwd_pen_guards{\"${k}_oob\"})";
      }
  }
  if (@guard_msgs) {
      print "    [PENALTY GUARD] Declenchements Forward : " . join(", ", @guard_msgs) . "\n";
  }

  
  # Finaliser la barre Forward / Finalize Forward bar
  # print "  [Forward] $forwardSetCount combinaisons Forward trouvees sur $innerForwardCount amorces F1.\n";

  # Check if anything found
  if($forwardSetCount == 0) {
      print_zero_signature_diagnostic(1, $innerForwardCount, scalar(@$masterOuterF_r),
        $_fwd_rej_geometry, $_fwd_rej_spacing, $_fwd_rej_loopgap,
        $_fwd_rej_tm_inner_loop, $_fwd_rej_tm_loop_middle, $_fwd_rej_tm_inner_middle, $_fwd_rej_tm_middle_outer,
        $_fwd_min_delta_tm_inner_loop, $_fwd_min_delta_tm_loop_middle, $_fwd_min_delta_tm_inner_middle, $_fwd_min_delta_tm_middle_outer,
        $_fwd_min_span_needed, $signatureMaxLength, $maxTmDiff, $options_r);
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
  
  # --- B&B Initialization Reverse ---
  my $minS_innerToLoop_R = @$innerToLoopPenalties_r ? $min_val_f->(@$innerToLoopPenalties_r) * $innerToLoopPenaltyWeight : 0;
  my $minS_loopToMiddle_R = @$loopToMiddlePenalties_r ? $min_val_f->(@$loopToMiddlePenalties_r) * $loopToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_R = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $minS_innerToMiddle_R = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $rmq_middle_r = build_rmq($masterMiddleR_data_r, 2);
  my $rmq_outer_r  = build_rmq($masterOuterR_data_r, 2);
  my $rmq_loop_r   = build_rmq($masterLoopR_data_r, 2) if $includeLoopPrimers;
  my $min_P_outer_R = @$masterOuterR_data_r ? query_rmq($rmq_outer_r, 0, scalar(@$masterOuterR_data_r)-1) * $outerPenaltyWeight : 0;
  
  my $_sig_rev_pruned = 0;
  my $_sig_rev_evaluated = 0;
  my $_rev_rej_geometry = 0;
  my $_rev_rej_spacing = 0;
  my $_rev_rej_loopgap = 0;
  my $_rev_rej_tm_inner_loop = 0;
  my $_rev_rej_tm_loop_middle = 0;
  my $_rev_rej_tm_inner_middle = 0;
  my $_rev_rej_tm_middle_outer = 0;
  my $_rev_min_delta_tm_inner_loop = 999;
  my $_rev_min_delta_tm_loop_middle = 999;
  my $_rev_min_delta_tm_inner_middle = 999;
  my $_rev_min_delta_tm_middle_outer = 999;
  my $_rev_min_span_needed = 999999;
  
  my %_rev_pen_guards;
  my $rev_prog_dir = "$options_r->{'output_file'}_rev_prog_$$";
  $rev_prog_dir = "$options_r->{'output_file'}_rev_prog_$$" if ref($options_r);
  use File::Path qw(make_path remove_tree);
  remove_tree($rev_prog_dir) if -d $rev_prog_dir;
  make_path($rev_prog_dir);
  my $_sig_rev_t0   = time();
  my $_sig_rev_done = 0;
  my $_sig_rev_hits = 0;  # Nombre de signatures Reverse trouvees / Reverse signatures found
  my $pm_rev = LLNL::LAVA::ForkManager->new($options_r->{"threads"});
  my $num_rev_chunks = $pm_rev->{max_processes} * 12;
  $num_rev_chunks = 30 if $num_rev_chunks < 30;
  $num_rev_chunks = $innerReverseCount if $num_rev_chunks > $innerReverseCount;
  $num_rev_chunks = 1 if $num_rev_chunks < 1;
  my $rev_chunk_size = int(($innerReverseCount + $num_rev_chunks - 1) / $num_rev_chunks);
  $rev_chunk_size = 1 if $rev_chunk_size < 1;

  $pm_rev->run_on_finish(sub {
      my ($pid, $exit_code, $id, $exit_signal, $core_dump, $data_ref) = @_;
      if (defined $data_ref && ref($data_ref) eq 'HASH') {
          foreach my $idx (sort { $a <=> $b } keys %{$data_ref->{infos}}) {
              if (!defined $bestReverseInfos[$idx]) {
                  $reverseSetCount++;
              }
              $bestReverseInfos[$idx] = $data_ref->{infos}->{$idx};
              $bestReversePenalties[$idx] = $data_ref->{penalties}->{$idx};
          }
          $_sig_rev_hits += $data_ref->{hits} || 0;
          $_sig_rev_done += $data_ref->{done} || 0;
          $_sig_rev_pruned += $data_ref->{pruned} || 0;
          $_sig_rev_evaluated += $data_ref->{evaluated} || 0;
          
          $_rev_rej_geometry += $data_ref->{rej_geometry} || 0;
          $_rev_rej_spacing += $data_ref->{rej_spacing} || 0;
          $_rev_rej_loopgap += $data_ref->{rej_loopgap} || 0;
          $_rev_rej_tm_inner_loop += $data_ref->{rej_tm_inner_loop} || 0;
          $_rev_rej_tm_loop_middle += $data_ref->{rej_tm_loop_middle} || 0;
          $_rev_rej_tm_inner_middle += $data_ref->{rej_tm_inner_middle} || 0;
          $_rev_rej_tm_middle_outer += $data_ref->{rej_tm_middle_outer} || 0;
          
          foreach my $k (qw(min_tm_inner_loop min_tm_loop_middle min_tm_inner_middle min_tm_middle_outer min_span_needed)) {
              next unless defined $data_ref->{$k};
              if ($k eq 'min_tm_inner_loop') {
                  $_rev_min_delta_tm_inner_loop = $data_ref->{$k} if $data_ref->{$k} < $_rev_min_delta_tm_inner_loop;
              } elsif ($k eq 'min_tm_loop_middle') {
                  $_rev_min_delta_tm_loop_middle = $data_ref->{$k} if $data_ref->{$k} < $_rev_min_delta_tm_loop_middle;
              } elsif ($k eq 'min_tm_inner_middle') {
                  $_rev_min_delta_tm_inner_middle = $data_ref->{$k} if $data_ref->{$k} < $_rev_min_delta_tm_inner_middle;
              } elsif ($k eq 'min_tm_middle_outer') {
                  $_rev_min_delta_tm_middle_outer = $data_ref->{$k} if $data_ref->{$k} < $_rev_min_delta_tm_middle_outer;
              } elsif ($k eq 'min_span_needed') {
                  $_rev_min_span_needed = $data_ref->{$k} if $data_ref->{$k} < $_rev_min_span_needed;
              }
          }
          foreach my $k (qw(innerToLoop loopToMiddle innerToMiddle middleToOuter)) {
              $_rev_pen_guards{"${k}_neg"} += $data_ref->{pen_guards}->{"${k}_neg"} || 0;
              $_rev_pen_guards{"${k}_oob"} += $data_ref->{pen_guards}->{"${k}_oob"} || 0;
          }
      }
  });

  for (my $chunk_id = 0; $chunk_id < $num_rev_chunks; $chunk_id++) {
      $pm_rev->start($chunk_id) and next;
      
      my %chunk_infos = ();
      my %chunk_penalties = ();
      my $chunk_hits = 0;
      my $chunk_done = 0;
      my $chunk_pruned = 0;
      my $chunk_evaluated = 0;
      my $chunk_rej_geometry = 0;
      my $chunk_rej_spacing = 0;
      my $chunk_rej_loopgap = 0;
      my $chunk_rej_tm_inner_loop = 0;
      my $chunk_rej_tm_loop_middle = 0;
      my $chunk_rej_tm_inner_middle = 0;
      my $chunk_rej_tm_middle_outer = 0;
      my $chunk_min_delta_tm_inner_loop = 999;
      my $chunk_min_delta_tm_loop_middle = 999;
      my $chunk_min_delta_tm_inner_middle = 999;
      my $chunk_min_delta_tm_middle_outer = 999;
      my $chunk_min_span_needed = 999999;
      my $penalty_guard_innerToLoop_neg = 0;
      my $penalty_guard_innerToLoop_oob = 0;
      my $penalty_guard_loopToMiddle_neg = 0;
      my $penalty_guard_loopToMiddle_oob = 0;
      my $penalty_guard_innerToMiddle_neg = 0;
      my $penalty_guard_innerToMiddle_oob = 0;
      my $penalty_guard_middleToOuter_neg = 0;
      my $penalty_guard_middleToOuter_oob = 0;
      
      for(my $innerIndex = $chunk_id; $innerIndex < $innerReverseCount; $innerIndex += $num_rev_chunks)
      {
          $chunk_done++;
          my $innerInfo = $masterInnerR_r->[$innerIndex];
          my ($innerLocation, $innerLength, $innerPenalty, $innerTm) = @{$masterInnerR_data_r->[$innerIndex]};
          
          my $bestSetPenalty = 1000000;
          
          # 4.1 Calculate Search Bounds for Loop Primer (Reverse)
          my $searchStartAt = $innerLocation + 1 + $minPrimerSpacing;
          my $searchEndAt = $innerLocation + $signatureMaxLength - $innerLength - 20; 
          
          # Determine Loop List to iterate
          my $currentLoopList_r = $includeLoopPrimers ? $masterLoopR_r : $masterLoopR_r; 
          my $currentLoopData_r = $includeLoopPrimers ? $masterLoopR_data_r : $masterLoopR_data_r;
          my $loopCount = scalar(@{$currentLoopList_r});
          
          if (!$includeLoopPrimers) {
              $currentLoopData_r->[0]->[0] = $searchStartAt; 
          }

          for(my $i = 0; $i < $loopCount; $i++)
          {
              my $loopInfo = $currentLoopList_r->[$i];
              my ($loopLocation, $loopLength, $loopPenalty, $loopTm) = @{$currentLoopData_r->[$i]};
              
              if ($includeLoopPrimers) {
                  next if ($loopLocation < $searchStartAt);
                  last if ($loopLocation > $searchEndAt);
                  my $diff = abs($innerTm - $loopTm);
                  if ($diff > $maxTmDiff) {
                      $chunk_min_delta_tm_inner_loop = $diff if $diff < $chunk_min_delta_tm_inner_loop;
                      $chunk_rej_tm_inner_loop++;
                      next;
                  }
              }
              
              # 4.2 Calculate Search Bounds for Middle Primer (Reverse)
              my $middleStartAt = $loopLocation + $loopLength + $minPrimerSpacing - $middlePrimerMaxLength + 1;
              # [restauré depuis a6098f5 : contrainte loopMinGap/inner supprimee par 55328b2]
              if ($includeLoopPrimers) {
                  my $altMiddleStartAt = $innerLocation + ($loopMinGap + 1);
                  $middleStartAt = $altMiddleStartAt if ($altMiddleStartAt > $middleStartAt);
              } else {
                  my $altMiddleStartAt = $innerLocation + $minPrimerSpacing;
                  $middleStartAt = $altMiddleStartAt if ($altMiddleStartAt > $middleStartAt);
              }
              my $middleEndAt = $searchEndAt;
              
              my $innerToLoopDistance = $loopLocation - ($innerLocation + 1);
              
              
              my $m_start = binary_search_first_ge($masterMiddleR_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleR_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_r, $m_start, $m_end) * $middlePenaltyWeight;
                  my $innerToLoopDistance = $includeLoopPrimers ? $loopLocation - ($innerLocation + 1) : 0;
                  
                  my $pen_innerToLoop = penaltyAt($innerToLoopPenalties_r, $innerToLoopDistance, 'innerToLoop');
                  if ($pen_innerToLoop < 0) { $penalty_guard_innerToLoop_neg++; }
                  if ($pen_innerToLoop == 9999) { $penalty_guard_innerToLoop_oob++; }
                  
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + 
                                     ($includeLoopPrimers ? $pen_innerToLoop * $innerToLoopPenaltyWeight : 0);
                  my $min_S_to_mid = $includeLoopPrimers ? $minS_loopToMiddle_R : $minS_innerToMiddle_R;
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_R + $min_S_to_mid + $minS_middleToOuter_R >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleR_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleR_data_r->[$j]};
                          
                          if ($includeLoopPrimers) {
                              my $loopToMiddleDistance = ($middleLocation - $middleLength + 1) - ($loopLocation + $loopLength);
                              if ($loopToMiddleDistance < 0) {
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $needed = $minPrimerSpacing - $loopToMiddleDistance;
                              if ($needed > 0) {
                                  my $span = $middleLocation - ($innerLocation - $innerLength) + $needed;
                                  $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              if ($middleLocation - $loopMinGap < $innerLocation + $innerLength) { $chunk_rej_loopgap++; next; }
                              my $diffLM = abs($loopTm - $midTm);
                              if ($diffLM > $maxTmDiff && !($loopInfo->hasTag("is_fixed") || $middleInfo->hasTag("is_fixed"))) {
                                  $chunk_min_delta_tm_loop_middle = $diffLM if $diffLM < $chunk_min_delta_tm_loop_middle;
                                  $chunk_rej_tm_loop_middle++;
                                  next;
                              }
                          } else {
                              my $innerToMiddleDistance = ($middleLocation - $middleLength) - $innerLocation;
                              if ($innerToMiddleDistance < 0) {
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $needed = $minPrimerSpacing - $innerToMiddleDistance;
                              if ($needed > 0) {
                                  my $span = $middleLocation - ($innerLocation - $innerLength) + $needed;
                                  $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                  $chunk_rej_spacing++;
                                  next;
                              }
                              my $diffIM = abs($innerTm - $midTm);
                              if ($diffIM > $maxTmDiff && !($innerInfo->hasTag("is_fixed") || $middleInfo->hasTag("is_fixed"))) {
                                  $chunk_min_delta_tm_inner_middle = $diffIM if $diffIM < $chunk_min_delta_tm_inner_middle;
                                  $chunk_rej_tm_inner_middle++;
                                  next;
                              }
                          }
                          
                          my $outerStartAt = $middleLocation + $minPrimerSpacing - $outerPrimerMaxLength + 1;
                          my $outerEndAt = $searchEndAt;
                          
                          my $loopToMiddleDistance = $includeLoopPrimers ? ($middleLocation - $middleLength + 1) - ($loopLocation + $loopLength) : 0;
                          my $innerToMiddleDistance = ($middleLocation - $middleLength) - $innerLocation;
                          
                          my $o_start = binary_search_first_ge($masterOuterR_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterR_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_r, $o_start, $o_end) * $outerPenaltyWeight;
                              
                              my $pen_mid = $includeLoopPrimers ? penaltyAt($loopToMiddlePenalties_r, $loopToMiddleDistance, 'loopToMiddle') : penaltyAt($innerToMiddlePenalties_r, $innerToMiddleDistance, 'innerToMiddle');
                              if ($pen_mid < 0) { if($includeLoopPrimers) { $penalty_guard_loopToMiddle_neg++; } else { $penalty_guard_innerToMiddle_neg++; } }
                              if ($pen_mid == 9999) { if($includeLoopPrimers) { $penalty_guard_loopToMiddle_oob++; } else { $penalty_guard_innerToMiddle_oob++; } }

                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($pen_mid * ($includeLoopPrimers ? $loopToMiddlePenaltyWeight : $innerToMiddlePenaltyWeight));
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_R >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterR_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterR_data_r->[$k]};
                                      
                                      my $middleToOuterDistance = ($outerLocation - $outerLength) - $middleLocation;
                                      if ($middleToOuterDistance < 0) {
                                          $chunk_rej_spacing++;
                                          next;
                                      }
                                      my $needed = $minPrimerSpacing - $middleToOuterDistance;
                                      if ($needed > 0) {
                                          my $span = $outerLocation - ($innerLocation - $innerLength) + $needed;
                                          $chunk_min_span_needed = $span if $span < $chunk_min_span_needed;
                                          $chunk_rej_spacing++;
                                          next;
                                      }
                                      my $diffMO = abs($midTm - $outTm);
                                      if ($diffMO > $maxTmDiff) {
                                          $chunk_min_delta_tm_middle_outer = $diffMO if $diffMO < $chunk_min_delta_tm_middle_outer;
                                          $chunk_rej_tm_middle_outer++;
                                          next;
                                      }
                                      
                                      my $pen_out = penaltyAt($middleToOuterPenalties_r, $middleToOuterDistance, 'middleToOuter');
                                      if ($pen_out < 0) { $penalty_guard_middleToOuter_neg++; }
                                      if ($pen_out == 9999) { $penalty_guard_middleToOuter_oob++; }
                                      
                                      my $spacingPenalty = ($includeLoopPrimers ? $pen_innerToLoop * $innerToLoopPenaltyWeight : 0) +
                                                           ($pen_mid * ($includeLoopPrimers ? $loopToMiddlePenaltyWeight : $innerToMiddlePenaltyWeight)) +
                                                           ($pen_out * $middleToOuterPenaltyWeight);
                                      
                                      my $primer3Penalty = $innerPenalty * $innerPenaltyWeight + ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                      
                                      my $detailStr = "";
                                      if ($includeLoopPrimers) {
                                          $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                                                ($pen_innerToLoop * $innerToLoopPenaltyWeight), ($pen_mid * $loopToMiddlePenaltyWeight), ($pen_out * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($loopPenalty * $loopPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($pen_mid * $innerToMiddlePenaltyWeight), ($pen_out * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      }
                                      
                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $bestSetPenalty) {
                                          $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                                          $chunk_infos{$innerIndex} = $includeLoopPrimers ? [$loopInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo];
                                          $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                                          $bestSetPenalty = $currentSetPenalty;
                                      }
                                  }
                              }
                          }
                      }
                  }
              }

          } # End Loop
          
          # Intra-chunk progress reporting
          if ($chunk_done % 5 == 0) {
              my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";
              if (open(my $fh, '>', $prog_file_me)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n";
                  close($fh);
              }
                  
              my $total_done = 0;
              my $total_hits = 0;
              foreach my $f (glob("$rev_prog_dir/chunk_*.prog")) {
                  if (open(my $r, '<', $f)) {
                      my $line = <$r>; close($r);
                      next unless defined $line;
                      chomp $line;
                      my ($d, $h) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                  }
              }
                      
              if ($_LAVA_IS_TTY || 1) {
                  my $elapsed = time() - $_sig_rev_t0 + 0.001;
                  my $eta = ($total_done < $innerReverseCount) ? int(($innerReverseCount - $total_done) / ($total_done / $elapsed)) : 0;
                  my $rate = $total_done / $elapsed;
                  printf("[LAVA-PROGRESS] Signatures Reverse|%d|%d|Sig: %d|%.1f it/s|%d\r", $total_done, $innerReverseCount, $total_hits, $rate, $eta);
                  my $old_h = select(STDOUT); $| = 1; select($old_h);
              }
          }
} # End Inner chunk loop
      
      my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";
      if (open(my $fh, '>', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n"; close($fh); }
      
      $pm_rev->finish(0, {
          infos => \%chunk_infos,
          penalties => \%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
          rej_geometry => $chunk_rej_geometry,
          rej_spacing => $chunk_rej_spacing,
          rej_loopgap => $chunk_rej_loopgap,
          rej_tm_inner_loop => $chunk_rej_tm_inner_loop,
          rej_tm_loop_middle => $chunk_rej_tm_loop_middle,
          rej_tm_inner_middle => $chunk_rej_tm_inner_middle,
          rej_tm_middle_outer => $chunk_rej_tm_middle_outer,
          min_tm_inner_loop => $chunk_min_delta_tm_inner_loop,
          min_tm_loop_middle => $chunk_min_delta_tm_loop_middle,
          min_tm_inner_middle => $chunk_min_delta_tm_inner_middle,
          min_tm_middle_outer => $chunk_min_delta_tm_middle_outer,
          min_span_needed => $chunk_min_span_needed,
          pen_guards => {
              innerToLoop_neg => $penalty_guard_innerToLoop_neg,
              innerToLoop_oob => $penalty_guard_innerToLoop_oob,
              loopToMiddle_neg => $penalty_guard_loopToMiddle_neg,
              loopToMiddle_oob => $penalty_guard_loopToMiddle_oob,
              innerToMiddle_neg => $penalty_guard_innerToMiddle_neg,
              innerToMiddle_oob => $penalty_guard_innerToMiddle_oob,
              middleToOuter_neg => $penalty_guard_middleToOuter_neg,
              middleToOuter_oob => $penalty_guard_middleToOuter_oob,
          }
      });
  } # End chunks
  $pm_rev->wait_all_children();
  use File::Path qw(remove_tree);
  remove_tree($rev_prog_dir) if -d $rev_prog_dir;
  
  if ($_sig_rev_pruned + $_sig_rev_evaluated > 0) {
      my $pct = ($_sig_rev_pruned / ($_sig_rev_pruned + $_sig_rev_evaluated)) * 100;
      print sprintf("    [Reverse B&B] Elagage: %.2f%%\n", $pct);
  }

  @guard_msgs = ();
  foreach my $k (qw(innerToLoop loopToMiddle innerToMiddle middleToOuter)) {
      if ($_rev_pen_guards{"${k}_neg"}) {
          push @guard_msgs, "$k(neg: $_rev_pen_guards{\"${k}_neg\"})";
      }
      if ($_rev_pen_guards{"${k}_oob"}) {
          push @guard_msgs, "$k(oob: $_rev_pen_guards{\"${k}_oob\"})";
      }
  }
  if (@guard_msgs) {
      print "    [PENALTY GUARD] Declenchements Reverse : " . join(", ", @guard_msgs) . "\n";
  }

  
  # Finaliser la barre Reverse / Finalize Reverse bar
  # print "  [Reverse] $reverseSetCount combinaisons Reverse trouvees sur $innerReverseCount amorces B1.\n";

  if($reverseSetCount == 0) {
      print_zero_signature_diagnostic(0, $innerReverseCount, scalar(@$masterOuterR_r),
        $_rev_rej_geometry, $_rev_rej_spacing, $_rev_rej_loopgap,
        $_rev_rej_tm_inner_loop, $_rev_rej_tm_loop_middle, $_rev_rej_tm_inner_middle, $_rev_rej_tm_middle_outer,
        $_rev_min_delta_tm_inner_loop, $_rev_min_delta_tm_loop_middle, $_rev_min_delta_tm_inner_middle, $_rev_min_delta_tm_middle_outer,
        $_rev_min_span_needed, $signatureMaxLength, $maxTmDiff, $options_r);
      exit 0;
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
          # next if ($inner_gap > 100); # Too far apart (Inner Gap Limit)

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
  my $total_sigs_to_validate = scalar(@{$allFoundSignatures_r});
  print "Validating and calculating coverage for $total_sigs_to_validate signatures...\n";
  my $validated_count = 0;

  my $val_pm = LLNL::LAVA::ForkManager->new($options_r->{"threads"});
  my $actual_threads = $val_pm->{max_processes};
  my %validation_results;

  $val_pm->run_on_finish(sub {
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_r) = @_;
      if (defined($data_r) && ref($data_r) eq 'ARRAY') {
          foreach my $res (@$data_r) {
              my ($idx, $cov, $status, $final_ids_r, $primer_cov_r) = @$res;
              $validation_results{$idx} = {
                  coverage => $cov,
                  status   => $status,
                  final_ids => $final_ids_r,
                  primer_cov => $primer_cov_r
              };
          }
      }
  });

  # Chunking
  my $val_chunk_size = POSIX::ceil($total_sigs_to_validate / ($actual_threads * 4)); # Entrelacement
  $val_chunk_size = 100 if $val_chunk_size < 100;
  
  my @val_chunks;
  for(my $i = 0; $i < $total_sigs_to_validate; $i += $val_chunk_size) {
      my $end = $i + $val_chunk_size - 1;
      $end = $total_sigs_to_validate - 1 if $end >= $total_sigs_to_validate;
      push @val_chunks, [$i, $end];
  }
  
  # Distribuer en round-robin
  my @val_worker_batches;
  for(my $i = 0; $i < $actual_threads; $i++) {
      push @val_worker_batches, [];
  }
  for(my $i = 0; $i < scalar(@val_chunks); $i++) {
      my $worker_idx = $i % $actual_threads;
      push @{$val_worker_batches[$worker_idx]}, $val_chunks[$i];
  }

  for(my $w = 0; $w < $actual_threads; $w++) {
      my $batch = $val_worker_batches[$w];
      next if scalar(@$batch) == 0;
      
      $val_pm->start and next;
      
      my @results_for_worker;
      foreach my $chunk (@$batch) {
          my ($start, $end) = @$chunk;
          for(my $idx = $start; $idx <= $end; $idx++) {
              my $signature = $allFoundSignatures_r->[$idx];
              
              my ($final_ids_r, $coverage, $status) = calculateSignatureIntersection(
                  $signature, 
                  $inputMSA->num_sequences(), 
                  $signatureCommonTargetMinPercent,
                  $includeLoopPrimers,
                  "loop"
              );
              
              my $primer_cov_r = [];
              eval { $primer_cov_r = $signature->getTag("primer_coverage_details"); };
              
              push @results_for_worker, [$idx, $coverage, $status, $final_ids_r, $primer_cov_r];
          }
      }
      $val_pm->finish(0, \@results_for_worker);
  }
  
  $val_pm->wait_all_children;

  # Ré-appliquer les résultats dans le parent
  for(my $idx = 0; $idx < $total_sigs_to_validate; $idx++) {
      my $signature = $allFoundSignatures_r->[$idx];
      if (exists $validation_results{$idx}) {
          my $res = $validation_results{$idx};
          $signature->setTag("signature_intersection_ids", $res->{final_ids});
          $signature->setTag("signature_coverage_percent", sprintf("%.2f", $res->{coverage}));
          $signature->setTag("signature_target_count", scalar(@{$res->{final_ids}}));
          $signature->setTag("validation_status", $res->{status});
          $signature->setTag("primer_coverage_details", $res->{primer_cov});
          $validated_count++;
          
          if ($validated_count % 1000 == 0 || $validated_count == $total_sigs_to_validate) {
              print "[LAVA-PROGRESS] Validated $validated_count / $total_sigs_to_validate signatures...\n";
          }
      }
  }

  print "Validation complete.\n";


  my @allSignatures = 
    map {$_->[0]}
    sort { $a->[1] <=> $b->[1] || $a->[0]->getStartLocation() <=> $b->[0]->getStartLocation() || $a->[0]->getLength() <=> $b->[0]->getLength() || $a->[0]->getLocationSummary() cmp $b->[0]->getLocationSummary() }
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
    print "\n🔍 Analyse / Analysis of combinaisons sur les $num_signatures signatures finales après réduction...\n";
    
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
      print "Analyse / Analysis of combinaisons sauvegardée dans: $combinations_file\n\n";
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

# --- Branch & Bound Helpers ---
sub build_rmq {
    my ($data_r, $col_idx) = @_;
    my $n = scalar(@$data_r);
    return [] if $n == 0;
    my $log_n = int(log($n) / log(2)) + 1;
    my @st;
    for my $i (0 .. $n - 1) {
        $st[$i][0] = $data_r->[$i]->[$col_idx];
    }
    for my $j (1 .. $log_n) {
        my $len = 1 << ($j - 1);
        for my $i (0 .. $n - 1) {
            if ($i + $len < $n) {
                my $a = $st[$i][$j - 1];
                my $b = $st[$i + $len][$j - 1];
                $st[$i][$j] = ($a < $b) ? $a : $b;
            } else {
                $st[$i][$j] = $st[$i][$j - 1];
            }
        }
    }
    return \@st;
}

sub query_rmq {
    my ($st_r, $L, $R) = @_;
    return 1000000 if (!defined $st_r || scalar(@$st_r) == 0 || $L > $R || $L < 0);
    $R = scalar(@$st_r) - 1 if $R >= scalar(@$st_r);
    my $j = int(log($R - $L + 1) / log(2));
    my $len = 1 << $j;
    my $a = $st_r->[$L][$j];
    my $b = $st_r->[$R - $len + 1][$j];
    return ($a < $b) ? $a : $b;
}

sub binary_search_first_ge {
    my ($data_r, $target_loc) = @_;
    my $L = 0;
    my $R = scalar(@$data_r) - 1;
    my $ans = -1;
    while ($L <= $R) {
        my $mid = $L + (($R - $L) >> 1);
        if ($data_r->[$mid]->[0] >= $target_loc) {
            $ans = $mid;
            $R = $mid - 1;
        } else {
            $L = $mid + 1;
        }
    }
    return $ans;
}

sub binary_search_last_le {
    my ($data_r, $target_loc) = @_;
    my $L = 0;
    my $R = scalar(@$data_r) - 1;
    my $ans = -1;
    while ($L <= $R) {
        my $mid = $L + (($R - $L) >> 1);
        if ($data_r->[$mid]->[0] <= $target_loc) {
            $ans = $mid;
            $L = $mid + 1;
        } else {
            $R = $mid - 1;
        }
    }
    return $ans;
}


sub print_zero_signature_diagnostic {
    my ($is_forward, $innerCount, $outerCount, 
        $rej_geometry, $rej_spacing, $rej_loopgap, 
        $rej_tm_inner_loop, $rej_tm_loop_middle, $rej_tm_inner_middle, $rej_tm_middle_outer,
        $min_tm_inner_loop, $min_tm_loop_middle, $min_tm_inner_middle, $min_tm_middle_outer,
        $min_span_needed, $signatureMaxLength, $maxTmDiff, $opts_r) = @_;

    my $outBase = defined $opts_r->{'output_file'} ? $opts_r->{'output_file'} : "lava_run";
    my $diagFile = "${outBase}_diagnostic.txt";
    open(my $dfh, '>', $diagFile) or warn "Cannot open $diagFile\n";
    
    my $side = $is_forward ? "Forward" : "Reverse";
    
    my $msg = "================================================================================\n";
    $msg .= "DIAGNOSTIC LAVA : AUCUNE DEMI-SIGNATURE $side TROUVÉE\n";
    $msg .= "================================================================================\n\n";

    if ($innerCount == 0 || $outerCount == 0) {
        $msg .= "[CAS 1] PROBLÈME EN AMONT : Pools d'amorces vides.\n";
        $msg .= "Aucune amorce n'a passé la validation initiale (0 candidat généré).\n";
        $msg .= "Le problème se situe avant l'assemblage : vos séquences sont probablement trop\n";
        $msg .= "divergentes dans cette région, ou les contraintes sont trop strictes.\n\n";
        $msg .= "Pistes de résolution :\n";
        $msg .= "  - Abaisser --min_primer_coverage (actuel : " . (defined $opts_r->{'min_primer_coverage'} ? $opts_r->{'min_primer_coverage'} : 100) . "%)\n";
        $msg .= "  - Augmenter --max_total_degenerate_bases (actuel : " . (defined $opts_r->{'max_total_degenerate_bases'} ? $opts_r->{'max_total_degenerate_bases'} : 0) . ")\n";
        $msg .= "  - Augmenter --max_tolerated_mismatches (actuel : " . (defined $opts_r->{'max_tolerated_mismatches'} ? $opts_r->{'max_tolerated_mismatches'} : 0) . ")\n";
        $msg .= "  - Élargir la plage de Tm (tm_min, tm_max) des amorces individuelles.\n";
    } else {
        $msg .= "[CAS 2] PROBLÈME D'ASSEMBLAGE : Incompatibilité des contraintes.\n";
        $msg .= "Les pools sont non vides (ex: $innerCount amorces F1c/B1c), mais aucune combinaison\n";
        $msg .= "ne respecte toutes les contraintes de température et d'espacement.\n\n";
        
        my @causes;
        push @causes, { name => "Thermodynamique (Inner/Loop)", rej => $rej_tm_inner_loop, min => $min_tm_inner_loop, param => "max_tm_diff", type => "tm" } if $rej_tm_inner_loop > 0;
        push @causes, { name => "Thermodynamique (Loop/Middle)", rej => $rej_tm_loop_middle, min => $min_tm_loop_middle, param => "max_tm_diff", type => "tm" } if $rej_tm_loop_middle > 0;
        push @causes, { name => "Thermodynamique (Inner/Middle)", rej => $rej_tm_inner_middle, min => $min_tm_inner_middle, param => "max_tm_diff", type => "tm" } if $rej_tm_inner_middle > 0;
        push @causes, { name => "Thermodynamique (Middle/Outer)", rej => $rej_tm_middle_outer, min => $min_tm_middle_outer, param => "max_tm_diff", type => "tm" } if $rej_tm_middle_outer > 0;
        
        push @causes, { name => "Espacement / Géométrie (Trop court)", rej => $rej_spacing, min => $min_span_needed, param => "signature_max_length", type => "span" } if $rej_spacing > 0;
        push @causes, { name => "Espacement (Loop Gap)", rej => $rej_loopgap, min => 0, param => "loop_min_gap", type => "gap" } if $rej_loopgap > 0;
        push @causes, { name => "Contraintes Géométriques Générales", rej => $rej_geometry, min => 0, param => "search_range", type => "geom" } if $rej_geometry > 0;
        
        @causes = sort { $b->{rej} <=> $a->{rej} } @causes;
        
        my $total_rej = $rej_geometry + $rej_spacing + $rej_loopgap + $rej_tm_inner_loop + $rej_tm_loop_middle + $rej_tm_inner_middle + $rej_tm_middle_outer;
        
        if (@causes && $total_rej > 0) {
            my $primary = $causes[0];
            my $pct = sprintf("%.1f", ($primary->{rej} / $total_rej) * 100);
            $msg .= "-> Cause principale : $primary->{name} ($pct% des rejets)\n";
            
            if ($primary->{type} eq 'tm') {
                $msg .= "   - Meilleur écart atteignable : " . sprintf("%.2f", $primary->{min}) . " °C.\n";
                $msg .= "   - Votre limite --max_tm_diff est : " . sprintf("%.2f", $maxTmDiff) . " °C.\n";
                my $suggest = int($primary->{min} + 1.5);
                if ($suggest > 10) {
                    $msg .= "   => SUGGESTION : Les séquences sont très divergentes. Une valeur > 10 °C n'est pas recommandée.\n";
                    $msg .= "      Envisagez de segmenter la cible ou de relâcher les paramètres d'amont.\n";
                } else {
                    $msg .= "   => SUGGESTION : Essayez --max_tm_diff $suggest\n";
                }
            } elsif ($primary->{type} eq 'span') {
                $msg .= "   - Empan minimal nécessaire : $primary->{min} nt.\n";
                $msg .= "   - Votre --signature_max_length est : $signatureMaxLength nt.\n";
                my $suggest = $primary->{min} + 20;
                $msg .= "   => SUGGESTION : Essayez --signature_max_length $suggest\n";
            } elsif ($primary->{type} eq 'gap') {
                $msg .= "   - Collision avec le loop_min_gap détectée.\n";
                $msg .= "   => SUGGESTION : Essayez de réduire --loop_min_gap ou d'augmenter l'espace total.\n";
            }
            
            if (scalar(@causes) > 1) {
                my $sec = $causes[1];
                my $pct2 = sprintf("%.1f", ($sec->{rej} / $total_rej) * 100);
                $msg .= "\n-> Cause secondaire : $sec->{name} ($pct2% des rejets)\n";
                if ($sec->{type} eq 'tm') {
                    $msg .= "   - Meilleur écart atteignable : " . sprintf("%.2f", $sec->{min}) . " °C.\n";
                } elsif ($sec->{type} eq 'span') {
                    $msg .= "   - Empan minimal nécessaire : $sec->{min} nt.\n";
                }
            }
        }
    }
    
    if ($is_forward) {
        $msg .= "\nLe scan Reverse n'a pas été exécuté car aucune demi-signature Forward n'a été trouvée.\n";
    }
    $msg .= "================================================================================\n";
    
    print $msg;
    print $dfh $msg if $dfh;
    close($dfh) if $dfh;
}
