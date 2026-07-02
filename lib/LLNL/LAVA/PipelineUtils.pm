################################################################################
#
# LLNL::LAVA::PipelineUtils - Fonctions utilitaires partagées du pipeline LAVA / Shared utility functions of LAVA pipeline
# Shared utility functions for the LAVA LAMP primer design pipeline
#
# Ce module factorise les fonctions communes entre lava_stem_primer.pl
# et lava_loop_primer.pl, éliminant la duplication de code. / and lava_loop_primer.pl, eliminating code duplication.
# This module factors out common functions from lava_stem_primer.pl
# and lava_loop_primer.pl, eliminating code duplication.
#
# Phase 34-36 code audit (2026).
#
# Copyright (c) 2026, Cheikh Talibouya <cheikhtalibouya.toure04@gmail.com | cheikhtalibouya.toure@pasteur.sn>. All rights reserved.
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the conditions of the BSD 3-Clause License are met.
# See the LICENSE file at the root of this project for full terms.
# Voir le fichier LICENSE à la racine du projet pour les termes complets.
#
################################################################################

package LLNL::LAVA::PipelineUtils;

use strict;
use warnings;
use Carp;

use Exporter 'import';
our @EXPORT_OK = qw(
  buildReversePrimers
  analyzeAll
  enumeratePairs
  buildMetricsArray
  reducePairInfosByPenalty
  reducePrimersByOverlap
  reduceSignaturesByOverlap
  flattenInfoData
  buildBigMerge
  calculateSignatureIntersection
  createPerSignatureFiles
  createAmplificationFiles
  analyzeSignatureCombinations
  generateCombinations
  calculateDynamicPairLengths
  reducePrimersByWindow
  buildNativeReversePool
);

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";
use LLNL::LAVA::PrimerSet::PCRPair;
use Bio::SeqIO;
use POSIX qw(floor);
use Time::HiRes qw(time);
use IO::Handle;
STDERR->autoflush(1);
# Detection TTY pour barre en place (comme tqdm) / TTY detection for in-place bar (like tqdm)
our $_LAVA_IS_TTY = -t STDERR ? 1 : 0;

# Auto-detection de Term::ProgressBar (equivalent tqdm) / Auto-detect Term::ProgressBar (tqdm equivalent)
# Si le module n'est pas installe, on utilise une barre ASCII maison sans dependance externe.
# If the module is not installed, a built-in ASCII bar is used with no external dependency.
our $HAS_TERM_PROGRESSBAR = eval { require Term::ProgressBar; Term::ProgressBar->import(); 1 } || 0;

#-------------------------------------------------------------------------------
# _make_progress_bar($total, $label)
# Cree une barre de progression Term::ProgressBar ou un objet fallback ASCII.
# Creates a Term::ProgressBar or a fallback ASCII progress bar object.
#-------------------------------------------------------------------------------
sub _make_progress_bar {
  my ($total, $label) = @_;
  $label //= "Traitement";

  if ($HAS_TERM_PROGRESSBAR && -t STDOUT) {
    # Mode riche : Term::ProgressBar avec ETA / Rich mode: Term::ProgressBar with ETA
    my $bar = Term::ProgressBar->new({
      name   => $label,
      count  => $total,
      ETA    => 'linear',
      remove => 0,
      fh     => \*STDERR,
    });
    $bar->minor(0);
    return { type => 'term', bar => $bar, total => $total, done => 0 };
  } else {
    # Mode fallback : barre ASCII maison / Fallback mode: built-in ASCII bar
    my $t0 = time();
    return {
      type   => 'ascii',
      total  => $total,
      done   => 0,
      label  => $label,
      t0     => $t0,
      last_print => 0,
    };
  }
}

#-------------------------------------------------------------------------------
# _update_progress($bar_r, $done, $extra_info)
# Met a jour la barre. Appeler a chaque iteration (freq limitee internalement).
# Updates the bar. Call at every iteration (frequency limited internally).
# $extra_info = hashref optionnel : { strict=>N, degen=>N, rejected=>N }
#-------------------------------------------------------------------------------
sub _update_progress {
  my ($bar_r, $done, $extra_r) = @_;
  $bar_r->{done} = $done;
  my $total = $bar_r->{total};

  if ($bar_r->{type} eq 'term') {
    $bar_r->{bar}->update($done);
    return;
  }

  # Fallback ASCII : afficher toutes les 200 iterations ou a 100%
  # Fallback ASCII: print every 200 iterations or at 100%
  return if ($done % 200 != 0 && $done != $total);

  my $now     = time();
  my $elapsed = $now - $bar_r->{t0} + 0.001;
  my $pct     = $total > 0 ? int($done / $total * 100) : 0;
  my $filled  = int($pct / 5);   # 20 segments
  my $empty   = 20 - $filled;
  my $bar_str = "#" x $filled . "-" x $empty;

  my $eta_str = "";
  if ($done > 0 && $done < $total) {
    my $rate    = $done / $elapsed;
    my $remain  = ($total - $done) / $rate;
    $eta_str = sprintf(" ETA:%ds", int($remain));
  }

  my $extra_str = "";
  if ($extra_r) {
    my @parts;
    push @parts, "OK:"   . ($extra_r->{strict}   // 0) if exists $extra_r->{strict};
    push @parts, "DEG:"  . ($extra_r->{degen}    // 0) if exists $extra_r->{degen};
    push @parts, "REJ:"  . ($extra_r->{rejected} // 0) if exists $extra_r->{rejected};
    $extra_str = " | " . join(" ", @parts) if @parts;
  }

  # Emission vers STDOUT pour Flask (toujours, pas seulement en TTY)
  # Emit to STDOUT for Flask (always, not only in TTY)
  my $rate_str = ($done > 0 && $elapsed > 0) ? sprintf('%.0f it/s', $done / $elapsed) : '? it/s';
  my $eta_val  = ($done > 0 && $done < $total && $elapsed > 0)
                 ? int(($total - $done) / ($done / $elapsed)) : 0;
  printf("[LAVA-PROGRESS] %s|%d|%d|%s|%s|%d\n",
         $bar_r->{label} // 'Reverse Validation', $done, $total,
         $extra_str, $rate_str, $eta_val);
}

sub _finish_progress {
  my ($bar_r) = @_;
  if ($bar_r->{type} eq 'term') {
    $bar_r->{bar}->update($bar_r->{total});
  } else {
    # Finaliser la ligne / Finalize the line
    _update_progress($bar_r, $bar_r->{total});
    # Progression finale envoyee via LAVA-PROGRESS stdout / Final progress sent via LAVA-PROGRESS stdout
  }
}


################################################################################
# FONCTIONS EXPORTÉES / EXPORTED FUNCTIONS
################################################################################

=head1 NAME

LLNL::LAVA::PipelineUtils - Fonctions utilitaires partagées entre STEM et LOOP

=head1 DESCRIPTION

Ce module centralise les fonctions de pipeline qui sont identiques entre
les scripts lava_stem_primer.pl et lava_loop_primer.pl. Cela garantit
qu'une correction de bug est automatiquement appliquée aux deux pipelines.

This module centralizes pipeline functions that are identical between
lava_stem_primer.pl and lava_loop_primer.pl. This ensures that a bug fix
is automatically applied to both pipelines.

=cut

#-------------------------------------------------------------------------------


=head2 buildReversePrimers (DEPRECATED)

B<DEPRECATED> — Cette fonction est conservée pour compatibilité ascendante uniquement.
Elle ne doit plus être appelée dans les nouveaux scripts.

B<DEPRECATED> — This function is retained for backward compatibility only.
It must not be called in new scripts.

Utiliser buildNativeReversePool à la place, qui génère nativement les amorces
reverse via Primer3 sur le reverse-complément du MSA, avec validation IUPAC complète.

Use buildNativeReversePool instead, which natively generates reverse primers
via Primer3 on the reverse-complement of the MSA, with full IUPAC validation.

=cut

sub buildReversePrimers
{
  # DEPRECATED: Remplacer par buildNativeReversePool / Replace with buildNativeReversePool
  # Conservé pour compatibilité ascendante uniquement / Retained for backward compatibility only
  my ($primers_r) = @_;

  my @reversePrimers = ();
  foreach my $primer(@{$primers_r})
  {
    my $minusStrandPrimer = $primer->clone();
    $minusStrandPrimer->reverseComplement();
    push(@reversePrimers, $minusStrandPrimer);
  }

  return @reversePrimers;
}

#-------------------------------------------------------------------------------

=head2 buildNativeReversePool

  Option B - Generation native des amorces Reverse.
  Genere les amorces du brin moins directement via Primer3 sur le RC de l alignement,
  puis les valide contre les sequences RC. Cela garantit que la zone 3 des amorces
  Reverse est correctement protegee (pas de base degeneree au 3 prime).

  Native Reverse primer generation (Option B).
  Generates minus-strand primers directly via Primer3 on the RC of the alignment,
  then validates them against RC sequences. This ensures the 3-prime zone of
  Reverse primers is correctly protected (no degenerate base at 3-prime).

  Parametres / Parameters:
    enumerator        - Objet Primer3Conserved configure
    alignment         - MSA BioPerl (Bio::SimpleAlign)
    min_match_percent - Seuil de concordance stricte
    min_iupac_percent - Seuil de couverture IUPAC
    min_primer_coverage - Couverture minimale par amorce
    maxTotalDegen, maxConsecDegen, max3PrimeDegen - Contraintes de degenerescence
    maxToleratedMismatches - Mismatches toleres hors zone 3prime
    threePrimeZoneSize - Taille de la zone 3prime stricte
    minBaseFrequency  - Frequence minimale pour inclure une base

=cut

sub buildNativeReversePool {
  my ($enumerator, $alignment,
      $min_match_percent, $min_iupac_percent, $min_primer_coverage,
      $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen,
      $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency,
      $checkPrimerMismatchTolerance_ref, $isIUPACCompatible_ref, $rev_comp_ref) = @_;

  # --- 1. Construire le RC de l alignement complet ---
  # Build the RC of the full alignment
  my $alignmentLength = $alignment->length();
  print "  [NativeReverse] Longueur alignement / Alignment length: $alignmentLength nt\n";

  # Creer un nouvel alignement avec les sequences RC / Create new alignment with RC sequences
  my $rcAlignment = Bio::SimpleAlign->new();
  my $seqIndex = 0;
  foreach my $seq ($alignment->each_seq()) {
    my $seqStr = $seq->seq();
    $seqStr = uc($seqStr);
    # Complement des bases (avec gaps preserves pour l entropie) / Complement bases (gaps preserved for entropy)
    (my $comp = $seqStr) =~ tr/ATCGatcgRYMKSWHBVDNrymkswhbvdn-/TAGCtagcYRKMSWDVBHNyrkmswdvbhn-/;
    my $rcStr = reverse($comp);
    my $rcSeq = Bio::LocatableSeq->new(
      -id  => $seq->id() . "_RC",
      -seq => $rcStr,
    );
    $rcAlignment->add_seq($rcSeq);
    $seqIndex++;
  }
  print "  [NativeReverse] RC MSA construit avec $seqIndex sequences / RC MSA built with $seqIndex sequences\n";

  # --- 2. Preparer les sequences RC pour la validation ---
  # Prepare RC sequences for validation (gaps/N remplace par N)
  my @rcSequences = ();
  foreach my $seq ($rcAlignment->each_seq()) {
    my $s = $seq->seq();
    $s = uc($s);
    $s =~ s/[^ATCG]/N/g;
    push @rcSequences, $s;
  }

  # --- 3. Lancer Primer3 sur le RC de l alignement ---
  # Run Primer3 on the RC alignment (generates native minus-strand primers)
  my @candidatePrimers = $enumerator->getOligos($rcAlignment);
  print "  [NativeReverse] Primer3 a genere / generated " . scalar(@candidatePrimers) . " candidats sur RC MSA\n";

  # --- 4. Valider chaque candidat contre les sequences RC ---
  # Validate each candidate against RC sequences (SENSE orientation vs RC targets)
  my @validatedPrimers = ();
  my $strict_count   = 0;
  my $degen_count    = 0;
  my $rejected_count = 0;

  my $nb_rev_candidates = scalar(@candidatePrimers);
  print "INFO: [NativeReverse] Analyse de $nb_rev_candidates amorces candidates Reverse...\n";
  print "  - Strict match: ${min_match_percent}% | IUPAC: ${min_iupac_percent}% | Coverage: ${min_primer_coverage}%\n";
  print "  - MaxDegen: $maxTotalDegen total / $max3PrimeDegen au 3prime\n";

  # Barre de progression / Progress bar (Term::ProgressBar ou fallback ASCII)
  my $rev_progress = _make_progress_bar($nb_rev_candidates, "Reverse Validation");

  foreach my $primer (@candidatePrimers) {
    my $posInRC  = $primer->location();  # Position 0-indexee dans RC(Seq1)
    my $length   = $primer->length();
    my $origSeq  = $primer->sequence();

    # Validation contre les sequences RC avec tous les parametres identiques aux Forward
    # Validation against RC sequences with identical parameters as Forward primers
    my ($finalSeq, $coveragePct, $isDegenerate, $compatibleIds) =
      $checkPrimerMismatchTolerance_ref->(
        \@rcSequences, $posInRC, $length, $origSeq,
        $min_match_percent, $min_iupac_percent, $min_primer_coverage,
        $maxTotalDegen, $maxConsecDegen, $max3PrimeDegen,
        $maxToleratedMismatches, $threePrimeZoneSize, $minBaseFrequency
      );

    # --- 5. Convertir la position RC -> position genome original ---
    # RC position p -> original genome position: alignmentLength - 1 - p (5' du brin -)
    my $genomicLocation = $alignmentLength - 1 - $posInRC;

    # Meme logique de rejet/acceptation que getOligosWithMismatchTolerance
    # Same rejection/acceptance logic as getOligosWithMismatchTolerance
    if ($coveragePct >= $min_primer_coverage) {
      # Creer l objet Oligo natif en orientation minus / Create native Oligo in minus orientation
      my $validatedPrimer = $primer->clone();
      $validatedPrimer->sequence($finalSeq);
      $validatedPrimer->location($genomicLocation);   # 5' du brin -, coordonnee originale
      $validatedPrimer->setTag("strand", "minus");    # Brin moins natif

      if ($isDegenerate) {
        $validatedPrimer->setTag("is_degenerate", 1);
        $validatedPrimer->setTag("original_sequence", $origSeq);
        $validatedPrimer->setTag("iupac_coverage", sprintf("%.1f", $coveragePct));
        $validatedPrimer->setTag("compatible_sequence_ids", $compatibleIds);
        $degen_count++;
        print "REVERSE DEGENERATE acceptee - PosRC: $posInRC -> GenomPos: $genomicLocation, Couv: " .
              sprintf("%.1f", $coveragePct) . "%, Seq: $finalSeq\n";
        _update_progress($rev_progress, $strict_count+$degen_count+$rejected_count,
          { strict=>$strict_count, degen=>$degen_count, rejected=>$rejected_count });
      } else {
        $validatedPrimer->setTag("is_degenerate", 0);
        $validatedPrimer->setTag("iupac_coverage", "100.0");
        $validatedPrimer->setTag("compatible_sequence_ids", $compatibleIds);
        $strict_count++;
        print "REVERSE STRICT acceptee   - PosRC: $posInRC -> GenomPos: $genomicLocation, Couv: 100.0%, Seq: $finalSeq\n";
        _update_progress($rev_progress, $strict_count+$degen_count+$rejected_count,
          { strict=>$strict_count, degen=>$degen_count, rejected=>$rejected_count });
      }

      push @validatedPrimers, $validatedPrimer;
    } else {
      $rejected_count++;
      print "REVERSE REJECTED           - PosRC: $posInRC -> GenomPos: $genomicLocation, Couv: " .
            sprintf("%.1f", $coveragePct) . "% < ${min_primer_coverage}%\n";
      _update_progress($rev_progress, $strict_count+$degen_count+$rejected_count,
        { strict=>$strict_count, degen=>$degen_count, rejected=>$rejected_count });
    }
  }

  _finish_progress($rev_progress);
  print "  [NativeReverse] Resultats / Results:\n";
  print "    - Strictes / Strict: $strict_count\n";
  print "    - Degenerees / Degenerate: $degen_count\n";
  print "    - Rejetees / Rejected: $rejected_count\n";
  print "    - Total valide / Total validated: " . scalar(@validatedPrimers) . "/" . scalar(@candidatePrimers) . "\n\n";

  return @validatedPrimers;
}

#-------------------------------------------------------------------------------

=head2 analyzeAll

  Analyse tous les primers avec l'analyseur donné (PrimerAnalyzer).
  Analyzes all primers with the given analyzer (PrimerAnalyzer).

=cut

sub analyzeAll
{
  my ($primers_r, $analyzer) = @_;

  my $measurements_r = [];
  foreach my $primer(@{$primers_r})
  {
    my $primerInfo = $analyzer->analyze($primer);
    push(@{$measurements_r}, $primerInfo);
  }

  return $measurements_r;
}

#-------------------------------------------------------------------------------

=head2 enumeratePairs

  Énumère toutes les paires Forward/Reverse compatibles dans les limites de longueur.
  Enumerates all compatible Forward/Reverse pairs within length constraints.

  Paramètres requis / Required parameters:
    sorted_forward_infos - Référence vers les infos forward triées par position
    sorted_reverse_infos - Référence vers les infos reverse triées par position
    max_length           - Longueur maximale de la paire

  Paramètre optionnel / Optional parameter:
    min_length           - Longueur minimale de la paire (défaut: 1)

=cut

sub enumeratePairs
{
  my ($paramHash_r) = @_;

  my $forwardInfos_r = optionRequired($paramHash_r, "sorted_forward_infos");
  my $reverseInfos_r = optionRequired($paramHash_r, "sorted_reverse_infos");
  my $maxLength = optionRequired($paramHash_r, "max_length");
  my $minLength = optionWithDefault($paramHash_r, "min_length", 1);

  my $pcrPairs_r = [];
  my $forwardInfoCount = scalar(@{$forwardInfos_r});
  my $reverseInfoCount = scalar(@{$reverseInfos_r});
  my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner loop
  for(my $i = 0; $i < $forwardInfoCount; $i++)
  {
    my $forwardInfo = $forwardInfos_r->[$i];
    my $forwardStart = $forwardInfo->getLocation();
    my $forwardEnd = $forwardStart + $forwardInfo->getLength() - 1;

    # Used to bound the upper end of the inner loop search
    my $maxReverseLocation = $forwardStart + $maxLength - 1;
    
    # Used to help bound the lower end of the inner loop search
    my $previousCompatibleIndexFound = $FALSE;

    for(my $j = $previousFirstCompatibleIndex; $j < $reverseInfoCount; $j++)
    {
      # Careful - these reverse primer locations are expressed as their
      # positive strand starts and ends.  This means the reverse primer's
      # "start" location IS NOT the 5' start of the primer.
      my $reverseInfo = $reverseInfos_r->[$j];
      my $reverseEnd = $reverseInfo->getLocation();
      my $reverseStart = $reverseEnd - $reverseInfo->getLength() + 1;

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

      # Enforce lower length boundary
      if(($reverseEnd - $forwardStart + 1) < $minLength)
      {
        next;
      }

      # Enforce upper length boundary by using the calculated max location
      if($reverseEnd > $maxReverseLocation)
      {
        last;
      }
         
      my $newPair = LLNL::LAVA::PrimerSet::PCRPair->new(
        {
          "forward_info" => $forwardInfo,
          "reverse_info" => $reverseInfo,
        });
      push(@{$pcrPairs_r}, $newPair);
    }
  }

  return $pcrPairs_r;
}

#-------------------------------------------------------------------------------

=head2 buildMetricsArray

  Construit un tableau de métriques (pénalité, position, longueur) pour chaque paire.
  Builds a metrics array (penalty, position, length) for each pair.

=cut

sub buildMetricsArray
{
  my ($pairInfos_r) = @_;

  my $pairCount = scalar(@{$pairInfos_r});

  my $metricsArray_r = [];
  for(my $i = 0; $i < $pairCount; $i++)
  { 
    my $pairInfo = $pairInfos_r->[$i];

    my $penalty = $pairInfo->getPenalty();
    my $pair = $pairInfo->getAnalyzedPair();

    my $start = $pair->getStartLocation();
    my $end = $pair->getEndLocation();
    my $length = $pair->getLength();

    my $forwardInfo = $pair->getForwardInfo();
    my $reverseInfo = $pair->getReverseInfo();

    my $forwardLength = $forwardInfo->getLength();
    my $reverseLength = $reverseInfo->getLength();

    my $clearAt = $start + $forwardLength;
    my $clearThrough = $end - $reverseLength;

    # Embed the pair in the metric array for independent handling
    $metricsArray_r->[$i] = [$penalty, $start, $end, $length, 
      $forwardLength, $reverseLength, $clearAt, $clearThrough];
  }

  return $metricsArray_r;
}

#-------------------------------------------------------------------------------

=head2 reducePairInfosByPenalty

  Réduit les paires en ne gardant que les N meilleures par pénalité.
  Reduces pairs by keeping only the N best by penalty score.

  Paramètres / Parameters:
    pair_infos - Référence vers les infos de paires
    max_pairs  - Nombre maximum de paires à conserver

=cut

sub reducePairInfosByPenalty
{
  my ($paramHash_r) = @_;

  my $pairInfos_r = optionRequired($paramHash_r, "pair_infos");
  my $maxPairs = optionRequired($paramHash_r, "max_pairs");

  # Sort the pairs by penalty / Tri par pénalité / Sort by penalty
  my @sortedInfos =
    map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    map {[$_, $_->getPenalty()] }
    @{$pairInfos_r}; 

  # Extract the top N pairs / Extraire les N meilleures paires
  my @bestInfos = ();
  my $infoCount = scalar(@sortedInfos);
  for(my $i = 0; ($i < $infoCount) && ($i < $maxPairs); $i++)
  {
    $bestInfos[$i] = $sortedInfos[$i];
  }

  # Re-sort by start location / Re-tri par position
  @bestInfos =
    map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    map {[$_, $_->getAnalyzedPair()->getStartLocation()]}
    @bestInfos;

  return \@bestInfos;
}

#-------------------------------------------------------------------------------

=head2 reducePrimersByOverlap

  Filtre les primers par chevauchement maximal autorisé.
  Filters primers by maximum allowed overlap percentage.

  Priorité aux primers avec la plus faible pénalité.
  Priority to primers with lowest penalty score.

  Paramètres / Parameters:
    max_overlap_percent     - Pourcentage de chevauchement maximal (0-100)
    info_sorted_by_location - Infos triées par position génomique
    info_sorted_by_penalty  - Infos triées par score de pénalité

=cut

sub reducePrimersByOverlap
{
  my ($paramHash_r) = @_;

  my $maxOverlapPercent = optionRequired($paramHash_r, "max_overlap_percent");
  my $sortedByLocation_r = optionRequired($paramHash_r, "info_sorted_by_location");
  my $sortedByPenalty_r = optionRequired($paramHash_r, "info_sorted_by_penalty");

  # Vérification de cohérence / Consistency check
  if(scalar(@{$sortedByLocation_r}) != scalar(@{$sortedByPenalty_r}))
  {
    confess("data error - lists don't have identical length! " .
      scalar(@{$sortedByLocation_r}) .
      " for location sorting, and " .
      scalar(@{$sortedByPenalty_r}) .
      " for penalty sorting.");
  }

  my $primerCount = scalar(@{$sortedByLocation_r});
  print " $primerCount->";
   
  # Short-cut if we're at 100% overlap
  if($maxOverlapPercent == 100)
  {
    my @primerList = ();
    foreach my $primerInfo(@{$sortedByPenalty_r})
    {
      push(@primerList, $primerInfo);
    }

    # Sort by location for return
    @primerList = 
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getLocation()] } 
      @primerList;
  
    print "" .
      scalar(@primerList) .
      ")\n";
  
    return \@primerList;
  } # End shortcut for 100% overlap

  # The combination of location/length should be unique for each primer 
  my %unavailablePrimers = ();
  my @primerList = ();

  # Pré-calcul des lookups pour performance / Pre-compute lookups for performance
  my @byPenaltyInfoLookup = ();
  for(my $infoIndex = 0; $infoIndex < $primerCount; $infoIndex++)
  {
    my $primerInfo = $sortedByPenalty_r->[$infoIndex];
    my $location = $primerInfo->getLocation();
    my $length = $primerInfo->getLength();
    my $primerTitle = "$location:$length"; 
    $byPenaltyInfoLookup[$infoIndex] = [$location, $length, $primerTitle];
  }
  my @byLocationInfoLookup = ();
  for(my $infoIndex = 0; $infoIndex < $primerCount; $infoIndex++)
  {
    my $primerInfo = $sortedByLocation_r->[$infoIndex];
    my $location = $primerInfo->getLocation();
    my $length = $primerInfo->getLength();
    my $primerTitle = "$location:$length"; 
    $byLocationInfoLookup[$infoIndex] = [$location, $length, $primerTitle];
  }

  # Priority to lowest penalty primers, knock out overlaps for each accepted primer. 
  for(my $outerIndex = 0; $outerIndex < $primerCount; $outerIndex++)
  {
    my $primerInfo = $sortedByPenalty_r->[$outerIndex];
    my ($location, $length, $primerTitle) = @{$byPenaltyInfoLookup[$outerIndex]};

    # Skip this primer if its unavailable
    if(exists $unavailablePrimers{$primerTitle})
    {
      next;
    }

    # Mark this new primer as unavailable now that it's chosen
    $unavailablePrimers{$primerTitle} = $TRUE;
    push(@primerList, $primerInfo);

    # Set upstream and downstream bounds for overlap checks
    my $upstreamStart = $location - (2 * $length);
    if($upstreamStart < 0)
    {
      $upstreamStart = 0;
    }
    my $downstreamEnd = $location + $length - 1;

    # Iterate over primers within bounds and check overlap percent
    # Mark all primers with more than the cutoff percent as unavailable
    for(my $innerIndex = 0; $innerIndex < $primerCount; $innerIndex++)
    {
      my $innerInfo = $sortedByLocation_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPrimerTitle) = @{$byLocationInfoLookup[$innerIndex]};

      if($innerLocation < $upstreamStart)
      {
	next;
      }

      if($innerLocation > $downstreamEnd)
      {
	last;
      }

      # Assume inner is longer, swap if not true 
      my $shorterStart = $location;
      my $shorterLength = $length;
      my $shorterEnd = $shorterStart + $shorterLength - 1;
      
      my $longerStart = $innerLocation;
      my $longerLength = $innerLength;
      my $longerEnd = $longerStart + $longerLength - 1;
     
      if($innerLength > $length)
      {
        $shorterStart = $innerLocation;
	$shorterLength = $innerLength;
	$shorterEnd = $shorterStart + $shorterLength - 1;

	$longerStart = $location;
	$longerLength = $length;
	$longerEnd = $longerStart + $longerLength - 1;
      }
    
      my $earlierStart = $shorterStart;
      if($earlierStart > $longerStart)
      {
	$earlierStart = $longerStart;
      }
      my $laterEnd = $shorterEnd;
      if($laterEnd < $longerEnd)
      {
	$laterEnd = $longerEnd;
      }

      
      my $overlapCount = ($shorterLength + $longerLength) - ($laterEnd - $earlierStart + 1);
      # Negative value is distance between oligos
      if($overlapCount < 0)
      {
	$overlapCount = 0;
      }
      my $overlapPercent = ($overlapCount / $shorterLength) * 100;

      if($overlapPercent > $maxOverlapPercent)
      {
        $unavailablePrimers{$innerPrimerTitle} = $TRUE;
      }
    } # End foreach primer by location
  } # End foreach primer by penalty

  # Sort primers by location for return
  @primerList = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @primerList;

  print "" .
    scalar(@primerList) .
    ")\n";

  return \@primerList;
}

#-------------------------------------------------------------------------------

=head2 reduceSignaturesByOverlap

  Réduit les signatures LAMP en éliminant celles qui se chevauchent trop.
  Reduces LAMP signatures by eliminating those with excessive overlap.

  Paramètres / Parameters:
    signatures          - Référence vers le tableau de signatures
    max_overlap_percent - Pourcentage de chevauchement maximal (défaut: 99)
    resolve_overlap_by  - Stratégie de résolution: "penalty" ou "coverage" (défaut: "penalty")
    sort_by_score       - Trier par score (défaut: $TRUE)
    sort_by_location    - Trier par position (défaut: $FALSE)

=cut

sub reduceSignaturesByOverlap
{
  my ($paramHash_r) = @_;

  my $signatures_r = optionRequired($paramHash_r, "signatures");
  my $maxOverlapPercent = optionWithDefault($paramHash_r, "max_overlap_percent", 99);
  my $resolveOverlapBy = optionWithDefault($paramHash_r, "resolve_overlap_by", "penalty");
  my $sortByScore = optionWithDefault($paramHash_r, "sort_by_score", $FALSE);
  my $sortByLocation = optionWithDefault($paramHash_r, "sort_by_location", $FALSE);
  if($sortByScore == $FALSE && $sortByLocation == $FALSE)
  {
    $sortByScore = $TRUE;
  }
  if($sortByScore == $TRUE && $sortByLocation == $TRUE)
  {
    $sortByLocation = $FALSE;
  }

  my $signatureCount= scalar(@{$signatures_r});
  print " Reducing signatures ha $signatureCount->";
  
  # Short-cut if we're at 100% overlap
  if($maxOverlapPercent == 100)
  {
    my @signatureList = ();
    foreach my $signature(@{$signatures_r})
    {
      push(@signatureList, $signature);
  }

    # Sort by location for return
    @signatureList = 
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getStartLocation()] } 
      @signatureList;
  
    print "" .
      scalar(@signatureList) .
      ")\n";
  
    return \@signatureList;
  } # End shortcut for 100% overlap

  my @byPenaltyLookup = ();
  my @byLocationLookup = ();
  for(my $sigIndex = 0; $sigIndex < $signatureCount; $sigIndex++)
  {
    my $signature = $signatures_r->[$sigIndex];

    my $penalty = $signature->getTag("lamp_penalty");
    my $coverage = $signature->getTag("signature_coverage_percent") || 0;
    my $startLocation = $signature->getStartLocation();
    my $length = $signature->getLength();

    $byPenaltyLookup[$sigIndex] = [$sigIndex, $penalty, $startLocation,
      $length, $coverage];
    $byLocationLookup[$sigIndex] = [$sigIndex, $penalty, $startLocation,
      $length, $coverage];
  }

  # Sort by resolve strategy / Tri selon la stratégie de résolution / Sort according to resolution strategy
  if ($resolveOverlapBy eq "coverage") {
    @byPenaltyLookup = 
      map {$_->[0]}
      sort {$b->[2] <=> $a->[2] || $a->[1] <=> $b->[1]}
      map {[$_, $_->[1], $_->[4]]} 
      @byPenaltyLookup;
  } else {
    @byPenaltyLookup = 
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->[1]]} 
      @byPenaltyLookup;
  }

  # Sort the other data lookup by location
  @byLocationLookup = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->[2]]}
    @byLocationLookup;

  my %unavailableSignatures = ();
  my @signatureInfoList = ();

  # Priority to lowest penalty, knock out overlaps for each accepted signature
  for(my $outerIndex = 0; $outerIndex < $signatureCount; $outerIndex++)
  {
    my ($originalIndex, $penalty, $location, $length) =
      @{$byPenaltyLookup[$outerIndex]};

    # Skip if unavailable
    if(exists $unavailableSignatures{$originalIndex})
    {
      next;
    }

    # Mark as chosen
    $unavailableSignatures{$originalIndex} = $TRUE;
    push(@signatureInfoList, $byPenaltyLookup[$outerIndex]);

    # Set bounds for overlap checks
    my $upstreamStart = $location - (5 * $length);
    if($upstreamStart < 0)
    {
      $upstreamStart = 0;
    }
    my $downstreamEnd = $location + $length - 1;

    # Iterate and check overlap percent
    for(my $innerIndex = 0; $innerIndex < $signatureCount; $innerIndex++)
    {
      my ($innerOriginalIndex, $innerPenalty, $innerLocation, $innerLength) =
        @{$byLocationLookup[$innerIndex]};

      if($innerLocation < $upstreamStart)
      {
	next;
      }

      if($innerLocation > $downstreamEnd)
      {
	last;
      }

      # Overlap calculation / Calcul du chevauchement
      my $shorterStart = $location;
      my $shorterLength = $length;
      my $shorterEnd = $shorterStart + $shorterLength - 1;
      
      my $longerStart = $innerLocation;
      my $longerLength = $innerLength;
      my $longerEnd = $longerStart + $longerLength - 1;
     
      if($innerLength > $length)
      {
        $shorterStart = $innerLocation;
	$shorterLength = $innerLength;
	$shorterEnd = $shorterStart + $shorterLength - 1;

	$longerStart = $location;
	$longerLength = $length;
	$longerEnd = $longerStart + $longerLength - 1;
      }
    
      my $earlierStart = $shorterStart;
      if($earlierStart > $longerStart)
      {
	$earlierStart = $longerStart;
      }
      my $laterEnd = $shorterEnd;
      if($laterEnd < $longerEnd)
      {
	$laterEnd = $longerEnd;
      }

      my $overlapCount = ($shorterLength + $longerLength) - ($laterEnd - $earlierStart + 1);
      if($overlapCount < 0)
      {
	$overlapCount = 0;
      }
      my $overlapPercent = ($overlapCount / $shorterLength) * 100;

      if($overlapPercent > $maxOverlapPercent)
      {
        $unavailableSignatures{$innerOriginalIndex} = $TRUE;
      }
    } # End foreach signature by location
  } # End foreach signature by penalty

  # Sort for return / Tri pour le retour
  my $sortingIndex = 2; # Assume sorting by location
  if($sortByScore == $TRUE)
  {
    $sortingIndex = 1;
  }
  @signatureInfoList = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->[$sortingIndex]]} 
    @signatureInfoList;

  # Build the answer array from the original signature set indexes
  my @signatureList = ();
  foreach my $info_r(@signatureInfoList)
  {
    my ($originalIndex, $penalty, $location, $length) = @{$info_r};
    push(@signatureList, $signatures_r->[$originalIndex]);
  }

  print "" .
    scalar(@signatureList) .
    ")\n";

  return \@signatureList;
}

#-------------------------------------------------------------------------------

=head2 flattenInfoData

  Aplatit les infos de primers en tableau de données [position, longueur, pénalité, Tm].
  Flattens primer info objects into data arrays [location, length, penalty, Tm].

  Paramètres / Parameters:
    info_set_ref - Référence vers le tableau d'infos de primers

=cut

sub flattenInfoData
{
  my ($paramHash_r) = @_;

  my $infoSet_r = optionRequired($paramHash_r, "info_set_ref");

  my $flattenedData_r = [];
  my $infoSetLength = scalar(@{$infoSet_r});
  for(my $infoIndex = 0; $infoIndex < $infoSetLength; $infoIndex++)
  {
    my $info = $infoSet_r->[$infoIndex];
    
    # Extraction robuste du Tm / Robust Tm extraction
    my $tm = 0;
    eval {
        my $oligo = $info->getAnalyzedPrimer();
        $tm = $oligo->getTag("primer3_tm");
    };
    if($@) { $tm = 0; } # Default if missing

    $flattenedData_r->[$infoIndex] = 
      [$info->getLocation(), $info->getLength(), $info->getPenalty(), $tm];
  }

  return $flattenedData_r;
}

#-------------------------------------------------------------------------------

=head2 buildBigMerge

  Architecture Single-Pass (LAVA 2026) commune aux scripts LOOP et STEM.
  Construit les listes maîtresses d'amorces (Master Lists) en une seule passe
  avec max_overlap_percent = 100, préservant toute la diversité combinatoire
  nécessaire à l'assemblage des signatures LAMP.
  
  Single-Pass Architecture (LAVA 2026) shared by LOOP and STEM scripts.
  Builds primer master lists in a single pass with max_overlap_percent = 100,
  preserving the full combinatorial diversity needed for LAMP signature assembly.
  The overlap reduction is intentionally deferred to the final signature reduction.

  Arguments (hash ref) :
    inner_f_loc   - \\@innerForwardInfoByLocation
    inner_f_pen   - \\@innerForwardInfoByPenalty
    inner_r_loc   - \\@innerReverseInfoByLocation
    inner_r_pen   - \\@innerReverseInfoByPenalty
    special_f_loc - \\@loopForwardInfoByLocation  (ou stemForward, ou [])
    special_f_pen - \\@loopForwardInfoByPenalty
    special_r_loc - \\@loopBackInfoByLocation     (ou stemBack, ou [])
    special_r_pen - \\@loopBackInfoByPenalty
    include_special - 1 ou 0 (activer les amorces Loop/STEM / enable Loop/STEM primers)
    middle_f_loc  - \\@middleForwardInfoByLocation
    middle_f_pen  - \\@middleForwardInfoByPenalty
    middle_r_loc  - \\@middleReverseInfoByLocation
    middle_r_pen  - \\@middleReverseInfoByPenalty
    outer_f_loc   - \\@outerForwardInfoByLocation
    outer_f_pen   - \\@outerForwardInfoByPenalty
    outer_r_loc   - \\@outerReverseInfoByLocation
    outer_r_pen   - \\@outerReverseInfoByPenalty

  Retourne un hash ref avec les clés / Returns a hash ref with keys:
    inner_f, inner_f_data, inner_r, inner_r_data,
    special_f, special_f_data, special_r, special_r_data,
    middle_f, middle_f_data, middle_r, middle_r_data,
    outer_f,  outer_f_data,  outer_r,  outer_r_data

=cut

#-------------------------------------------------------------------------------

=head2 reducePrimersByWindow

  Reduit le nombre de candidats amorces en ne gardant que les K meilleurs par
  fenetre genomique de W nucleotides. Preserves la diversite spatiale tout en
  eliminant les candidats redondants proches.

  Reduces primer candidates by keeping only the K best (by penalty) per
  genomic window of W nucleotides. Preserves spatial diversity while
  eliminating redundant nearby candidates.

  Params:
    primers_by_location : arrayref de PrimerInfo tries par position
    window_size         : largeur de fenetre en pb (ex: 10)
    max_per_window      : nb max de primers gardes par fenetre (ex: 3)

  Returns: arrayref de PrimerInfo tries par position

=cut

sub reducePrimersByWindow
{
  my ($primers_r, $window_size, $max_per_window) = @_;

  return $primers_r unless scalar(@{$primers_r}) > 0;

  # 1. Trier globalement les candidats d'amorces par pénalité croissante (les meilleures en premier)
  # 1. Sort all primer candidates globally by ascending penalty (best thermodynamic candidates first)
  my @sorted_by_penalty = sort {
    $a->getPenalty() <=> $b->getPenalty()
  } @{$primers_r};

  my @selected = ();

  # 2. Zones d'exclusion actives : tableau de structures {start => X, end => Y, count => C}
  # 2. Active exclusion zones: array of structures {start => X, end => Y, count => C}
  my @exclusion_zones = ();

  foreach my $primerInfo (@sorted_by_penalty) {
    my $location = $primerInfo->getLocation();

    # Vérifier si l'amorce candidate tombe dans une zone d'exclusion existante
    # Check if the candidate primer falls into an existing exclusion zone
    my $matched_zone_r = undef;
    foreach my $zone_r (@exclusion_zones) {
      if ($location >= $zone_r->{start} && $location <= $zone_r->{end}) {
        $matched_zone_r = $zone_r;
        last;
      }
    }

    if (defined $matched_zone_r) {
      # Si l'amorce est dans une zone d'exclusion, on ne la garde que si le quota max par zone n'est pas atteint
      # If the primer is in an exclusion zone, only keep it if the maximum quota per zone is not reached
      if ($matched_zone_r->{count} < $max_per_window) {
        push @selected, $primerInfo;
        $matched_zone_r->{count}++;
      }
    } else {
      # Si elle est en dehors de toute zone, on l'accepte et on crée une nouvelle zone d'exclusion centrée sur elle
      # If it's outside any zone, accept it and create a new exclusion zone centered around it
      push @selected, $primerInfo;

      my $half_window = int($window_size / 2);
      my $start = $location - $half_window;
      my $end   = $location + $half_window;

      push @exclusion_zones, {
        start => $start,
        end   => $end,
        count => 1
      };
    }
  }

  # 3. Réordonner la liste finale par position génomique croissante pour préserver la cohérence des scans géométriques
  # 3. Re-sort the final selected list by ascending genomic location to preserve consistency for geometric scans
  my @reduced = sort {
    $a->getLocation() <=> $b->getLocation()
  } @selected;

  return \@reduced;
}

#-------------------------------------------------------------------------------

sub buildBigMerge
{
  my ($p) = @_;

  # Sous-routine interne : reduire + aplatir une liste / Internal helper: reduce + flatten a list
  my $merge = sub {
    my ($loc_r, $pen_r) = @_;
    return [] unless scalar(@{$loc_r});
    my $result = reducePrimersByOverlap({
      "max_overlap_percent"     => 100,
      "info_sorted_by_location" => $loc_r,
      "info_sorted_by_penalty"  => $pen_r,
    });
    # Reduction spatiale si parametres fournis / Spatial reduction if parameters provided
    if ($p->{window_size} && $p->{max_per_window}) {
      my $before = scalar(@{$result});
      $result = reducePrimersByWindow($result, $p->{window_size}, $p->{max_per_window});
      my $after = scalar(@{$result});
      print "    [Fenetre W=$p->{window_size}nt, max=$p->{max_per_window}/fenetre] $before -> $after candidats\n" if $before != $after;
    }
    return $result;
  };

  my $masterInnerF_r  = $merge->($p->{inner_f_loc},  $p->{inner_f_pen});
  my $masterInnerR_r  = $merge->($p->{inner_r_loc},  $p->{inner_r_pen});
  my $masterMiddleF_r = $merge->($p->{middle_f_loc}, $p->{middle_f_pen});
  my $masterMiddleR_r = $merge->($p->{middle_r_loc}, $p->{middle_r_pen});
  my $masterOuterF_r  = $merge->($p->{outer_f_loc},  $p->{outer_f_pen});
  my $masterOuterR_r  = $merge->($p->{outer_r_loc},  $p->{outer_r_pen});

  # Amorces spéciales (Loop ou STEM) : seulement si activées / Special primers (Loop or STEM): only if enabled
  my $masterSpecialF_r = [];
  my $masterSpecialR_r = [];
  if ($p->{include_special}) {
    $masterSpecialF_r = $merge->($p->{special_f_loc}, $p->{special_f_pen});
    $masterSpecialR_r = $merge->($p->{special_r_loc}, $p->{special_r_pen});
  }

  return {
    inner_f        => $masterInnerF_r,
    inner_f_data   => flattenInfoData({"info_set_ref" => $masterInnerF_r}),
    inner_r        => $masterInnerR_r,
    inner_r_data   => flattenInfoData({"info_set_ref" => $masterInnerR_r}),
    special_f      => $masterSpecialF_r,
    special_f_data => flattenInfoData({"info_set_ref" => $masterSpecialF_r}),
    special_r      => $masterSpecialR_r,
    special_r_data => flattenInfoData({"info_set_ref" => $masterSpecialR_r}),
    middle_f       => $masterMiddleF_r,
    middle_f_data  => flattenInfoData({"info_set_ref" => $masterMiddleF_r}),
    middle_r       => $masterMiddleR_r,
    middle_r_data  => flattenInfoData({"info_set_ref" => $masterMiddleR_r}),
    outer_f        => $masterOuterF_r,
    outer_f_data   => flattenInfoData({"info_set_ref" => $masterOuterF_r}),
    outer_r        => $masterOuterR_r,
    outer_r_data   => flattenInfoData({"info_set_ref" => $masterOuterR_r}),
  };
}

#-------------------------------------------------------------------------------

=head2 generateCombinations

  Génère toutes les combinaisons de taille N à partir d'un tableau.
  Generates all combinations of size N from an array.

=cut

sub generateCombinations {
  my ($array_ref, $size) = @_;
  my @array = @{$array_ref};
  my @combinations = ();
  
  return ([]) if $size == 0;
  return map { [$_] } @array if $size == 1;
  
  for my $i (0 .. $#array - $size + 1) {
    my $element = $array[$i];
    my @remaining = @array[($i + 1) .. $#array];
    my @sub_combinations = generateCombinations(\@remaining, $size - 1);
    for my $sub_combination (@sub_combinations) {
      push @combinations, [$element, @{$sub_combination}];
    }
  }
  return @combinations;
}

#-------------------------------------------------------------------------------

=head2 analyzeSignatureCombinations

  Analyse toutes les combinaisons de signatures pour maximiser la couverture.
  Analyzes all signature combinations to maximize coverage.

=cut

sub analyzeSignatureCombinations {
  my ($signatures_ref, $total_sequences) = @_;
  
  my @signatures = @{$signatures_ref};
  my $num_signatures = scalar(@signatures);
  
  print "\n  ANALYSE DES COMBINAISONS DE SIGNATURES\n";
  print "   Nombre de signatures: $num_signatures\n";
  print "   Nombre total de sequences: $total_sequences\n\n";
  
  print "  REFERENCE DES SIGNATURES:\n";
  for my $i (0 .. $#signatures) {
    my $signature = $signatures[$i];
    my $coverage = $signature->getTag("signature_coverage_percent") || 0;
    my $penalty = $signature->getTag("lamp_penalty") || 0;
    my $target_count = $signature->getTag("signature_target_count") || 0;
    printf "   Sig%d: %d sequences (%.1f%%) - Penalite: %.1f\n", 
           $i + 1, $target_count, $coverage, $penalty;
  }
  print "\n";
  
  my %results_by_size = ();
  my $max_combination_size = ($num_signatures > 8) ? 4 : $num_signatures;
  
  if ($max_combination_size < $num_signatures) {
    print "   Limitation des combinaisons a $max_combination_size par $max_combination_size (sur $num_signatures signatures)\n";
  }
  
  for my $combination_size (1 .. $max_combination_size) {
    print "  COMBINAISONS $combination_size par $combination_size:\n";
    
    my @combinations = generateCombinations(\@signatures, $combination_size);
    my @combination_results = ();
    
    for my $combination_ref (@combinations) {
      my @combination = @{$combination_ref};
      my %union_sequences = ();
      my @signature_names = ();
      
      for my $i (0 .. $#combination) {
        my $signature = $combination[$i];
        my $global_index = -1;
        for my $j (0 .. $#signatures) {
          if ($signatures[$j] == $signature) {
            $global_index = $j + 1;
            last;
          }
        }
        my $coverage = $signature->getTag("signature_coverage_percent") || 0;
        my $penalty = $signature->getTag("lamp_penalty") || 0;
        my $sig_name = sprintf("Sig%d(%.1f%%/P%.1f)", $global_index, $coverage, $penalty);
        push @signature_names, $sig_name;
        
        my $amplified_sequences_ref;
        eval { $amplified_sequences_ref = $signature->getTag("signature_intersection_ids"); };
        if (!$@ && defined $amplified_sequences_ref) {
          for my $seq_id (@{$amplified_sequences_ref}) {
            $union_sequences{$seq_id} = 1;
          }
        }
      }
      
      my $union_count = scalar(keys %union_sequences);
      my $union_coverage = ($total_sequences > 0) ? ($union_count / $total_sequences) * 100 : 0.0;
      
      push @combination_results, {
        signatures => \@combination,
        signature_names => \@signature_names,
        union_sequences => [keys %union_sequences],
        union_count => $union_count,
        union_coverage => $union_coverage
      };
    }
    
    @combination_results = sort { $b->{union_coverage} <=> $a->{union_coverage} } @combination_results;
    $results_by_size{$combination_size} = \@combination_results;
    
    my $max_display = ($combination_size == 1) ? 10 : 5;
    my $displayed = 0;
    for my $result (@combination_results) {
      last if $displayed >= $max_display;
      my $names_str = join(" + ", @{$result->{signature_names}});
      printf "   %s: %d sequences (%.2f%%)\n", 
             $names_str, $result->{union_count}, $result->{union_coverage};
      $displayed++;
    }
    if (scalar(@combination_results) > $max_display) {
      printf "   ... et %d autres combinaisons\n", scalar(@combination_results) - $max_display;
    }
    print "\n";
  }
  
  return \%results_by_size;
}

#-------------------------------------------------------------------------------

=head2 calculateDynamicPairLengths

  Calcule dynamiquement les longueurs cibles des paires Middle et Inner
  a partir des contraintes de distance maximale entre niveaux de primers.
  Dynamically computes target pair lengths from max distance constraints.

  Parametres / Parameters:
    outer_pair_target  - Longueur cible de la paire Outer
    max_dist_outer_middle - Distance max F3-F2
    max_dist_middle_inner - Distance max F2-F1
    min_inner_pair_spacing - Espacement minimum F1-B1

  Retourne / Returns: ($middlePairTarget, $innerPairTarget)

=cut

sub calculateDynamicPairLengths {
  my ($outer_pair_target, $max_dist_outer_middle, $max_dist_middle_inner, $min_inner_pair_spacing) = @_;
  
  print "\nINFO: Calcul dynamique des longueurs cibles active.\n";
  print "INFO: Dynamic target length calculation activated.\n";
  
  my $middle_pair_target = $outer_pair_target - (2 * $max_dist_outer_middle);
  print "  -> Cible Middle Pair calculee : $middle_pair_target nt (distance max: $max_dist_outer_middle nt)\n";
  
  my $inner_pair_target = $middle_pair_target - (2 * $max_dist_middle_inner);
  print "  -> Cible Inner Pair calculee : $inner_pair_target nt (distance max: $max_dist_middle_inner nt)\n";
  
  if ($inner_pair_target < $min_inner_pair_spacing) {
    die "ERREUR: La cible Inner Pair ($inner_pair_target nt) < distance minimale F1-B1 ($min_inner_pair_spacing nt).\n" .
        "ERROR: Inner Pair target ($inner_pair_target nt) < minimum F1-B1 spacing ($min_inner_pair_spacing nt).\n" .
        "Veuillez ajuster vos contraintes de distance.\n";
  }
  print "--------------------------------------------------\n\n";
  
  return ($middle_pair_target, $inner_pair_target);
}

#-------------------------------------------------------------------------------

=head2 calculateSignatureIntersection

  Calcule l'intersection des sequences compatibles de tous les primers d'une signature.
  Computes the intersection of compatible sequences across all primers in a signature.
  
  Version unifiee STEM+LOOP : prend le meilleur des deux scripts.
  Unified STEM+LOOP version: takes the best from both scripts.

  Parametres / Parameters:
    $signature          - Objet signature LAMP
    $total_sequences    - Nombre total de sequences dans l'alignement
    $min_signature_coverage - Seuil minimal de couverture (defaut: 70%)
    $include_extra_primers  - Booleen: inclure les primers enrichis
    $extra_primer_type      - "stem" ou "loop"

=cut

sub calculateSignatureIntersection {
  my ($signature, $total_sequences, $min_signature_coverage, $include_extra_primers, $extra_primer_type) = @_;
  
  # Valeurs par defaut / Default values
  $min_signature_coverage = 70 unless defined $min_signature_coverage;
  $extra_primer_type = "" unless defined $extra_primer_type;
  
  my $type_label = uc($extra_primer_type) || "CLASSIC";
  my $mode_text = $include_extra_primers ? "LAMP+$type_label" : "LAMP classique (6 primers)";
  print "\n  VALIDATION DE SIGNATURE $mode_text (seuil: ${min_signature_coverage}%)\n";
  
  # Primers principaux : F3/B3 (outer), F2/B2 (middle), F1/B1 (inner)
  my @all_primers = ();
  my @primer_names = ("F3", "B3", "F2", "B2", "F1", "B1");
  
  my $outerPair = $signature->getOuterInfo()->getAnalyzedPair();
  my $middlePair = $signature->getMiddleInfo()->getAnalyzedPair();
  my $innerPair = $signature->getInnerInfo()->getAnalyzedPair();
  
  push @all_primers, $outerPair->getForwardInfo()->getAnalyzedPrimer();
  push @all_primers, $outerPair->getReverseInfo()->getAnalyzedPrimer();
  push @all_primers, $middlePair->getForwardInfo()->getAnalyzedPrimer();
  push @all_primers, $middlePair->getReverseInfo()->getAnalyzedPrimer();
  push @all_primers, $innerPair->getForwardInfo()->getAnalyzedPrimer();
  push @all_primers, $innerPair->getReverseInfo()->getAnalyzedPrimer();
  
  # Ajouter les primers enrichis (STEM ou LOOP) si demande
  # Add enriched primers (STEM or LOOP) if requested
  if ($include_extra_primers) {
    my $f_tag = "f" . lc($extra_primer_type) . "_info";
    my $b_tag = "b" . lc($extra_primer_type) . "_info";
    my $F_label = "F" . uc($extra_primer_type);
    my $B_label = "B" . uc($extra_primer_type);
    
    # Verification defensive par primer individuel (pattern STEM, plus robuste)
    # Defensive per-primer check (STEM pattern, more robust)
    my $f_info;
    eval { $f_info = $signature->getTag($f_tag); };
    if (!$@ && defined $f_info) {
      my $f_primer = $f_info->getAnalyzedPrimer();
      # Creer le tag par defaut si absent / Create default tag if missing
      eval { $f_primer->getTag("compatible_sequence_ids"); };
      if ($@) {
        print "   [FIX] $F_label sans tag compatible_sequence_ids - creation par defaut\n";
        my @all_seq_ids = (0 .. $total_sequences - 1);
        $f_primer->setTag("compatible_sequence_ids", \@all_seq_ids);
      }
      push @all_primers, $f_primer;
      push @primer_names, $F_label;
    }
    
    my $b_info;
    eval { $b_info = $signature->getTag($b_tag); };
    if (!$@ && defined $b_info) {
      my $b_primer = $b_info->getAnalyzedPrimer();
      eval { $b_primer->getTag("compatible_sequence_ids"); };
      if ($@) {
        print "   [FIX] $B_label sans tag compatible_sequence_ids - creation par defaut\n";
        my @all_seq_ids = (0 .. $total_sequences - 1);
        $b_primer->setTag("compatible_sequence_ids", \@all_seq_ids);
      }
      push @all_primers, $b_primer;
      push @primer_names, $B_label;
    }
    
    if (!defined $f_info && !defined $b_info) {
      print "   [WARN] $type_label primers demandes mais non trouves dans la signature\n";
    } else {
      print "   $type_label primers ajoutes a l'analyse\n";
    }
  }
  
  print "   Primers a analyser: " . scalar(@all_primers) . " (" . join(", ", @primer_names) . ")\n";
  return ([], 0.0, "Aucun primer disponible") if @all_primers == 0;
  
  # Phase 1: Verifier les sequences compatibles de chaque primer
  # Phase 1: Check compatible sequences for each primer
  my @primer_coverage_data = ();
  for my $i (0 .. $#all_primers) {
    my $primer = $all_primers[$i];
    my $primer_name = $primer_names[$i];
    
    my $compatible_sequences_ref;
    eval { $compatible_sequences_ref = $primer->getTag("compatible_sequence_ids"); };
    
    if ($@ || !defined $compatible_sequences_ref) {
      print "   [FAIL] $primer_name: PAS de tag 'compatible_sequence_ids'\n";
      return ([], 0.0, "Primer $primer_name sans sequences compatibles");
    }
    
    my $count = scalar(@{$compatible_sequences_ref});
    my $pct = ($total_sequences > 0) ? ($count / $total_sequences) * 100 : 0.0;
    print "   [OK] $primer_name: $count sequences compatibles (${pct}%)\n";
    
    push @primer_coverage_data, {
      name => $primer_name, sequences => $compatible_sequences_ref,
      count => $count, percent => $pct
    };
  }
  
  # Phase 2: Intersection successive / Successive intersection
  print "\n  CALCUL DE L'INTERSECTION:\n";
  my %intersection_set = map { $_ => 1 } @{$primer_coverage_data[0]->{sequences}};
  
  for my $i (1 .. $#primer_coverage_data) {
    my %primer_set = map { $_ => 1 } @{$primer_coverage_data[$i]->{sequences}};
    my %new_intersection = ();
    for my $seq_id (keys %intersection_set) {
      $new_intersection{$seq_id} = 1 if exists $primer_set{$seq_id};
    }
    %intersection_set = %new_intersection;
    last if scalar(keys %intersection_set) == 0;
  }
  
  # Phase 3: Validation finale / Final validation
  my @final_ids = keys %intersection_set;
  my $coverage_count = scalar(@final_ids);
  my $coverage_pct = ($total_sequences > 0) ? ($coverage_count / $total_sequences) * 100 : 0.0;
  
  print "\n  RESULTAT FINAL:\n";
  print "   Sequences amplifiees par TOUTES les amorces: $coverage_count/$total_sequences\n";
  printf "   Pourcentage de couverture: %.2f%%\n", $coverage_pct;
  
  my $validation_status;
  if ($coverage_pct >= $min_signature_coverage) {
    print "   [PASS] SIGNATURE VALIDEE (${coverage_pct}% >= ${min_signature_coverage}%)\n";
    $validation_status = "VALIDEE";
  } else {
    print "   [FAIL] SIGNATURE REJETEE (${coverage_pct}% < ${min_signature_coverage}%)\n";
    $validation_status = "REJETEE - Couverture insuffisante";
  }
  
  # Stocker les tags de validation / Store validation tags
  $signature->setTag("signature_intersection_ids", \@final_ids);
  $signature->setTag("signature_coverage_percent", sprintf("%.2f", $coverage_pct));
  $signature->setTag("signature_target_count", $coverage_count);
  $signature->setTag("validation_status", $validation_status);
  $signature->setTag("primer_coverage_details", \@primer_coverage_data);
  
  return (\@final_ids, $coverage_pct, $validation_status);
}

#-------------------------------------------------------------------------------

=head2 createPerSignatureFiles

  Cree un fichier par signature avec rapport de validation.
  Creates a per-signature file with validation report.
  
  Version unifiee : eval defensif (STEM) + header algorithme (LOOP).
  Unified version: defensive eval (STEM) + algorithm header (LOOP).

=cut

sub createPerSignatureFiles {
  my ($signatures_ref, $sequence_names_ref, $output_base_name, $primer_type) = @_;
  
  $primer_type = "LAMP" unless defined $primer_type;
  return unless @{$signatures_ref} > 0;
  
  print "\n  Creation des fichiers par signature individuelle...\n";
  
  my $signatures_dir = "${output_base_name}_signatures_individuelles";
  unless (-d $signatures_dir) {
    mkdir($signatures_dir) or die "Impossible de creer le dossier $signatures_dir: $!";
  }
  
  my $total_signatures = scalar(@{$signatures_ref});
  print "   Nombre de signatures finales: $total_signatures\n";
  
  for my $sig_index (0 .. $#{$signatures_ref}) {
    eval {
      my $signature = $signatures_ref->[$sig_index];
      
      # Recuperation securisee des tags (pattern STEM, defensif)
      # Safe tag retrieval (STEM pattern, defensive)
      my $amplified_sequences = [];
      my $coverage_percent = 0;
      my $coverage_count = 0;
      my $validation_status = "INCONNUE";
      my $primer_details = [];
      
      eval { $amplified_sequences = $signature->getTag("signature_intersection_ids"); };
      if ($@ || !defined $amplified_sequences) { 
        my $total_seqs = scalar(@{$sequence_names_ref});
        $amplified_sequences = [0 .. $total_seqs - 1];
      }
      eval { $coverage_percent = $signature->getTag("signature_coverage_percent"); };
      if ($@) { $coverage_percent = 100.0; }
      eval { $coverage_count = $signature->getTag("signature_target_count"); };
      if ($@) { $coverage_count = scalar(@{$amplified_sequences}); }
      eval { $validation_status = $signature->getTag("validation_status"); };
      if ($@) { $validation_status = "VALIDEE (par defaut)"; }
      eval { $primer_details = $signature->getTag("primer_coverage_details"); };
      if ($@) { $primer_details = []; }
      
      my $status_short = ($validation_status eq "VALIDEE") ? "VALID" : "REJECT";
      my $sig_filename = sprintf("%s/signature_%02d_%s_%d_seq_%.1fpc.txt", 
                                $signatures_dir, $sig_index + 1, $status_short, $coverage_count, $coverage_percent);
      
      open(my $sig_fh, '>', $sig_filename) or die "Cannot open $sig_filename: $!";
      
      # En-tete unifie / Unified header
      print $sig_fh "# ========================================\n";
      print $sig_fh "# SIGNATURE " . ($sig_index + 1) . " - RAPPORT DE VALIDATION\n";
      print $sig_fh "# ========================================\n";
      print $sig_fh "# Type: $primer_type\n";
      print $sig_fh "# Statut de validation: $validation_status\n";
      print $sig_fh "# Sequences amplifiees: $coverage_count\n";
      print $sig_fh "# Couverture: ${coverage_percent}%\n";
      print $sig_fh "# Date: " . localtime() . "\n";
      print $sig_fh "#\n";
      print $sig_fh "# ALGORITHME: Intersection de toutes les amorces\n";
      print $sig_fh "#\n";
      
      # Details de couverture / Coverage details
      if (@{$primer_details} > 0) {
        print $sig_fh "# COUVERTURE DES PRIMERS INDIVIDUELS:\n";
        for my $pd (@{$primer_details}) {
          printf $sig_fh "#   %s: %d sequences (%.1f%%)\n", $pd->{name}, $pd->{count}, $pd->{percent};
        }
        print $sig_fh "#\n";
      }
      
      print $sig_fh "# DETAILS DE LA SIGNATURE $primer_type:\n";
      
      # Primers de base (avec eval defensif) / Core primers (with defensive eval)
      for my $pname (["F3", "getF3"], ["B3", "getB3"], ["F2", "getF2"], 
                     ["B2", "getB2"], ["F1", "getF1"], ["B1", "getB1"]) {
        eval {
          my $method = $pname->[1];
          my $seq = $signature->$method();
          print $sig_fh "# $pname->[0]: " . ($seq || "N/A") . "\n";
        };
        if ($@) { print $sig_fh "# $pname->[0]: ERROR\n"; }
      }
      
      # Primers enrichis (STEM ou LOOP) / Enriched primers
      my $f_tag = "f" . lc($primer_type) . "_info";
      my $b_tag = "b" . lc($primer_type) . "_info";
      my $F_label = "F" . uc($primer_type);
      my $B_label = "B" . uc($primer_type);
      
      eval {
        my $f_info = $signature->getTag($f_tag);
        print $sig_fh "# $F_label: " . ($f_info ? $f_info->getSequence() : "N/A") . "\n";
      };
      if ($@) { print $sig_fh "# $F_label: N/A\n"; }
      eval {
        my $b_info = $signature->getTag($b_tag);
        print $sig_fh "# $B_label: " . ($b_info ? $b_info->getSequence() : "N/A") . "\n";
      };
      if ($@) { print $sig_fh "# $B_label: N/A\n"; }
      
      print $sig_fh "#\n";
      print $sig_fh "# === SEQUENCES AMPLIFIEES ===\n";
      
      if ($amplified_sequences && @{$amplified_sequences} > 0) {
        for my $seq_id (@{$amplified_sequences}) {
          if ($seq_id < @{$sequence_names_ref}) {
            print $sig_fh $sequence_names_ref->[$seq_id] . "\n";
          }
        }
      } else {
        print $sig_fh "# Aucune sequence amplifiee\n";
      }
      
      close($sig_fh);
      print "   [OK] Signature " . ($sig_index + 1) . ": $sig_filename\n";
    };
    if ($@) {
      print "   [ERR] Signature " . ($sig_index + 1) . ": $@\n";
    }
  }
  
  print "   Dossier: $signatures_dir ($total_signatures fichiers)\n\n";
  return $signatures_dir;
}

#-------------------------------------------------------------------------------

=head2 createAmplificationFiles

  Cree les fichiers FASTA des sequences amplifiees et exclues.
  Creates FASTA files for amplified and excluded sequences.
  
  Version LOOP (index-based avec bounds checking).
  LOOP version (index-based with bounds checking).

=cut

sub createAmplificationFiles {
  my ($signatures_ref, $sequence_objects_ref, $sequence_names_ref, $output_base_name) = @_;
  
  return unless @{$signatures_ref} > 0;
  
  # Collecter les IDs amplifies / Collect amplified IDs
  my %all_amplified_ids = ();
  for my $signature (@{$signatures_ref}) {
    my $intersection_ids = $signature->getTag("signature_intersection_ids");
    if (defined $intersection_ids && @{$intersection_ids} > 0) {
      for my $seq_id (@{$intersection_ids}) {
        $all_amplified_ids{$seq_id} = 1;
      }
    }
  }
  
  my @amplified_ids = keys %all_amplified_ids;
  my @excluded_ids = ();
  for my $i (0 .. $#{$sequence_objects_ref}) {
    push @excluded_ids, $i unless exists $all_amplified_ids{$i};
  }
  
  my $amplified_file = "${output_base_name}_amplified.fasta";
  my $excluded_file = "${output_base_name}_excluded.fasta";
  
  # Fichier FASTA amplifie / Amplified FASTA
  if (@amplified_ids > 0) {
    my $out = Bio::SeqIO->new(-file => ">$amplified_file", -format => 'fasta');
    for my $sid (@amplified_ids) {
      $out->write_seq($sequence_objects_ref->[$sid]) if $sid < @{$sequence_objects_ref};
    }
    print "Fichier cree: $amplified_file (" . scalar(@amplified_ids) . " sequences)\n";
  }
  
  # Fichier des noms / Names file
  my $names_file = "${output_base_name}_amplified_noms.txt";
  if (@amplified_ids > 0) {
    open(my $fh, '>', $names_file) or die "Cannot open $names_file: $!";
    for my $sid (@amplified_ids) {
      print $fh $sequence_objects_ref->[$sid]->display_id() . "\n" if $sid < @{$sequence_objects_ref};
    }
    close($fh);
    print "Fichier noms cree: $names_file (" . scalar(@amplified_ids) . " noms)\n";
  }
  
  # Fichier FASTA exclu / Excluded FASTA
  if (@excluded_ids > 0) {
    my $out = Bio::SeqIO->new(-file => ">$excluded_file", -format => 'fasta');
    for my $sid (@excluded_ids) {
      $out->write_seq($sequence_objects_ref->[$sid]) if $sid < @{$sequence_objects_ref};
    }
    print "Fichier cree: $excluded_file (" . scalar(@excluded_ids) . " sequences)\n";
  }
  
  # Statistiques / Statistics
  my $total = scalar(@{$sequence_objects_ref});
  my $amp_count = scalar(@amplified_ids);
  my $amp_pct = ($total > 0) ? ($amp_count / $total) * 100 : 0;
  
  print "\n=== RESUME D'AMPLIFICATION ===\n";
  printf "Total: %d | Amplifiees: %d (%.2f%%) | Exclues: %d (%.2f%%)\n",
         $total, $amp_count, $amp_pct, scalar(@excluded_ids), 100 - $amp_pct;
  print "==============================\n\n";
  
  return ($amplified_file, $excluded_file);
}

1;

__END__

=head1 AUTHOR

LAVA-DNA Fork (2026) - Code audit Phase 34-36

=head1 SEE ALSO

L<LLNL::LAVA::Core>, L<LLNL::LAVA::Validator>

=cut

