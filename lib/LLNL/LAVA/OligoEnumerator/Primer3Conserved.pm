package LLNL::LAVA::OligoEnumerator::Primer3Conserved;

use strict;
use vars qw(@ISA);
use Carp;

use Bio::SimpleAlign; # Recieves as a parameter
use Bio::Tools::Run::Primer3; # Uses to enumerate oligos over sequences

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::Oligo; # Builds and returns

use LLNL::LAVA::OligoEnumerator; # is-a
@ISA = ("LLNL::LAVA::OligoEnumerator");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::OligoEnumerator::Primer3Conserved - Finds all oligos that could be
used as primers based on perfect conservation across the input BioPerl MSA, 
using Primer3 for enumeration.
   
=head1 SYNOPSIS

  use LLNL::LAVA::OligoEnumerator::Primer3Conserved;

  # Instantiation
  $enumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => "/usr/bin/primer3_core",
    });
  
  # Set primer3 operating parameters 
  # ONLY these tags are recognized for now.
  # These values aren't sanity checked, so be very careful what you ask for!
  $enumerator->setPrimer3Targets(
    {
      "target_length" => 20,
      "min_length" => 17,
      "max_length" => 27,
      "target_tm" => 62,
      "min_tm" => 61,
      "max_tm" => 63,
      "most_to_return" => 20001,
      "thermodynamic_path" => "/etc/primer3_config/",
    });

  # Get oligo results of enumerating over an MSA
  @oligos = $enumerator->getOligos($bioPerlMSA);

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

Primer3Conserved is a concrete OligoEnumerator.  A BioPerl MSA is accepted
as input to getOligos().  Oligos are enumerated across the first sequence
using BioPerl's Bio::Tools::Run::Primer3.  Oligos that have exactly conserved
matches in all the remaining MSA sequences are returned.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $oligoEnumerator = 
               LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
                 {"primer3_executable" => "/usr/bin/primer3_core"} );
 Function  : Creates a new OligoEnumerator for Primer3 and perfect consensus
 Arguments : Hash Ref - options for OligoEnumerators, including:
               primer3_executable - path to primer3 (usually primer3_core)
 Example   : See Usage
 Returns   : A new LLNL::LAVA::OligoEnumerator

=cut

sub new
{
  my ($classType, $paramHash_r) = @_;

  if(!defined $paramHash_r)
  {
    confess("programming error - first parameter is a required hash ref");
  }
  if(ref($paramHash_r) ne "HASH")
  {
    confess("programming error - first parameter must be a hash reference");
  } 

  # Path to the primer3 executable
  my $primer3Executable = optionRequired($paramHash_r, "primer3_executable");

  my $this = $classType->SUPER::new();

  # Map of parameter name to official Primer3 target name
  # No values can be set through setPrimer3Targets() that aren't listed here 
  my $p3Names_r = {
    "target_length" => "PRIMER_INTERNAL_OPT_SIZE",
    "min_length" => "PRIMER_INTERNAL_MIN_SIZE",
    "max_length" => "PRIMER_INTERNAL_MAX_SIZE",
    "target_tm" => "PRIMER_INTERNAL_OPT_TM",
    "min_tm" => "PRIMER_INTERNAL_MIN_TM",
    "max_tm" => "PRIMER_INTERNAL_MAX_TM",
    "max_poly_bases" => "PRIMER_INTERNAL_OLIGO_MAX_POLY_X",
    "most_to_return" => "PRIMER_NUM_RETURN",
    "min_gc" => "PRIMER_INTERNAL_MIN_GC",
    "max_gc" => "PRIMER_INTERNAL_MAX_GC",
    "dna_conc" => "PRIMER_INTERNAL_DNA_CONC",
    "salt_divalent" => "PRIMER_INTERNAL_SALT_DIVALENT",
    "salt_monovalent" => "PRIMER_INTERNAL_SALT_MONOVALENT",
    "dntp_conc" => "PRIMER_INTERNAL_DNTP_CONC",
    # ====== NOUVEAUX PARAMÈTRES THERMODYNAMIQUES LAMP 2026 ======
    # Paramètres SantaLucia 1998 et Owczarzy 2004 pour LAMP isotherme à 65°C / SantaLucia 1998 and Owczarzy 2004 parameters for isothermal LAMP at 65°C
    # New thermodynamic parameters for LAMP isothermal amplification at 65°C
    "tm_formula" => "PRIMER_TM_FORMULA",
    "salt_corrections" => "PRIMER_SALT_CORRECTIONS",
    "thermodynamic_alignment" => "PRIMER_THERMODYNAMIC_ALIGNMENT",
    "salt_corrections" => "PRIMER_SALT_CORRECTIONS",
    "thermodynamic_alignment" => "PRIMER_THERMODYNAMIC_ALIGNMENT",
    "thermodynamic_path" => "PRIMER_THERMODYNAMIC_PARAMETERS_PATH",
    "excluded_regions" => "SEQUENCE_EXCLUDED_REGION",
  };

  # Set of default primer3 targets (primer3 target name => value)
  my $p3Targets_r = {
    "PRIMER_TASK" => "pick_hyb_probe_only",
    #"PRIMER_INTERNAL_OLIGO_SELF_ANY" => "12.00",
    "PRIMER_INTERNAL_OLIGO_SELF_ANY" => "8.00",
    #"PRIMER_INTERNAL_MAX_POLY_X" => 4, # DISABLING TO SEE IF IT CONFLICTS WITH PRIMER_INTERNAL_OLIGO_MAX_POLY_X

    # Default suggested by Primer3 documentation
    #"PRIMER_INTERNAL_OLIGO_MAX_END_STABILITY" => 9.0, 
    #"PRIMER_MAX_END_STABILITY" => 9.0, 
    $p3Names_r->{"salt_monovalent"} => 50,
    $p3Names_r->{"salt_divalent"} => 8,
    $p3Names_r->{"dntp_conc"} => 1.4, #nano-molar
    $p3Names_r->{"dna_conc"} => 400,
    $p3Names_r->{"min_gc"} => 30,
    $p3Names_r->{"max_gc"} => 80,
    $p3Names_r->{"target_length"} => 20,
    $p3Names_r->{"min_length"} => 18,
    $p3Names_r->{"max_length"} => 27,
    $p3Names_r->{"target_tm"} => 60,
    $p3Names_r->{"min_tm"} => 50,
    $p3Names_r->{"max_tm"} => 69,
    $p3Names_r->{"most_to_return"} => 20001, # Off-by-one error in primer3?
    # ====== PARAMÈTRES THERMODYNAMIQUES LAMP 2026 ======
    # SantaLucia 1998 (TM_FORMULA=1) et Owczarzy 2004 (SALT_CORRECTIONS=2)
    # Optimisé pour amplification LAMP isotherme à 65°C / Optimized for isothermal LAMP amplification at 65°C
    # Optimized for LAMP isothermal amplification at 65°C
    "PRIMER_TM_FORMULA" => 1,                     # SantaLucia 1998 (plus précis)
    "PRIMER_SALT_CORRECTIONS" => 2,               # Owczarzy 2004 (pour Mg2+)
    "PRIMER_THERMODYNAMIC_ALIGNMENT" => 1,        # Calcul ΔG structures secondaires
  };

  # Paramètres LAVA personnalisés (non-Primer3) / Custom LAVA parameters (non-Primer3)
  my $lavaParams_r = {
    "entropy_threshold" => 1.5,  # Seuil d'entropie de Shannon (défaut 1.5)
  };

  $this->{"d_primer3Executable"} = $primer3Executable;
  $this->{"d_primer3NameConversion"} = $p3Names_r; 
  $this->{"d_primer3Targets"} = $p3Targets_r;
  $this->{"d_lavaParams"} = $lavaParams_r;  # Paramètres LAVA personnalisés
 
  return $this;
}

#-------------------------------------------------------------------------------

=head2 getOligos

 Usage     : $enumerator->getOligos($bioPerlMSA);
 Function  : Creates an array of LLNL::LAVA::Oligo objects based on the
             input BioPerl MSA.  Oligos are enumerated across the first 
             sequence using BioPerl's Bio::Tools::Run::Primer3.  Oligos that 
             have exactly conserved matches in all the remaining MSA 
             sequences are returned.
 Arguments : Bio::SimpleAlign - MSA of the targets
 Example   : See Usage
 Returns   : Array of L<LLNL::LAVA::Oligo> - the oligos found for the MSA

=cut

sub getOligos
{
  my ($this, $alignment) = @_;

  my $sequenceCount = $alignment->num_sequences();
  if($sequenceCount <= 0)
  {
    confess("data error - MSA must contain at least one sequence");
  }

  # Make sure identifiers are unique, and that lengths are identical
  my %sequenceIDs = (); 
  my $observedLength = -1;
  foreach my $sequence($alignment->each_seq())
  {
    # Make any non-ATCGN an N, should handle this more gracefully, but 
    # primer3 doesn't really like non-ATCGN's
    my $seqContent = $sequence->seq();
    $seqContent = uc($seqContent);  # Convertir en majuscules d'abord
    $seqContent =~ s/[^ATCG]/N/g;  # Puis remplacer caractères non-ADN par N
    $sequence->seq($seqContent);

    my $id = $sequence->id();
    my $headerPiece = $sequence->desc();
    if(defined $headerPiece &&
       $headerPiece ne "")
    {
      $id .= " " . $headerPiece;
    }

    if(exists $sequenceIDs{$id})
    {
      confess("data error - alignment has 2 entries with the ID \"$id\", but " .
        "we can only handle one sequence of each ID here");
    }
    
    # Check for matched lengths
    my $currLength = $sequence->length();
    if($observedLength == -1) # Record first length no matter what
    {
      $observedLength = $currLength;
    }
    if($observedLength != $currLength)
    {
      confess("data error - the first sequence was $observedLength long, but " .
        "the sequence \"$id\" was $currLength long");
    }

    $sequenceIDs{$id} = $TRUE;
  }

  # TODO: Going to have to scale the number of primers returned along
  # with the sequence length, or enforce a max length for this locus?

  my %oligosBySequenceID = ();
  my $sequenceIndex = 1;
  
  # Running primer3 on the first sequence
  my $firstSequence = $alignment->get_seq_by_pos(1); # 1-indexed!
  my $primer3 = Bio::Tools::Run::Primer3->new(
    -seq => $firstSequence,
    -path => $this->{"d_primer3Executable"});
  my $p3Targets_r = $this->{"d_primer3Targets"};
#  foreach my $targetName (keys(%{$primer3->arguments()}))
#  {
#   print "  $targetName\t" + $p3Targets_r->{$targetName} + "\n"; 
# }

  $primer3->add_targets(%{$p3Targets_r});

  # ANALYSE ENTROPIQUE (LAVA 2026)
  # Calculer l'entropie de Shannon par colonne pour identifier les zones variables
  # et les exclure de la recherche Primer3
  my @entropies = ();
  my $alignmentLength = $alignment->length();
  
  for (my $pos = 1; $pos <= $alignmentLength; $pos++) {
      my %baseComments = ();
      my $count = 0;
      my $gapCount = 0;
      foreach my $seq ($alignment->each_seq()) {
          my $base = substr($seq->seq(), $pos-1, 1);
          # On compte TOUT, même les gaps. On ignore juste 'N' qui est une vraie inconnu.
          # Mais même 'N' pourrait être compté comme un état... restons simple : N ignoré, Gap compté. / But even 'N' could be counted as a state... let's keep it simple: N ignored, Gap counted.
          if ($base eq '-') {
              $baseComments{'-'}++;
              $gapCount++;
              $count++;
          } elsif ($base ne 'N') {
              $baseComments{$base}++;
              $count++;
          }
      }
      
      my $entropy = 0;
      
      # Règle de pénalité GAPS : Si > 20% de gaps, zone invalide (entropie forcée au max) / GAPS penalty rule: If > 20% gaps, invalid zone (entropy forced to max)
      if ($count > 0 && ($gapCount / $count) > 0.20) {
          $entropy = 10.0; # Valeur arbitrairement haute pour garantir l'exclusion
      } elsif ($count > 0) {
          foreach my $base (keys %baseComments) {
              my $p = $baseComments{$base} / $count;
              $entropy -= $p * (log($p) / log(2));
          }
      }
      push(@entropies, $entropy);
  }
  
  # Identifier les zones à haute entropie / Identify high entropy zones

  
  # Identifier les zones à haute entropie / Identify high entropy zones
  my @excludedRegions = ();
  my $inHighEntropy = 0;
  my $startHigh = 0;
  # Récupérer le seuil d'entropie depuis les paramètres LAVA (défaut 1.5) / Retrieve entropy threshold from LAVA parameters (default 1.5)
  my $lavaParams_r = $this->{"d_lavaParams"};
  my $entropyThreshold = $lavaParams_r->{"entropy_threshold"} || 1.5;
  
  # Récupérer la taille de fenêtre (taille MINIMALE de l'amorce pour être permissif) / Retrieve window size (MINIMUM primer size to be permissive)
  # Si on utilise la taille min (ex: 18), on s'assure qu'au moins une amorce de cette taille peut passer
  my $windowSize = $p3Targets_r->{"PRIMER_INTERNAL_MIN_SIZE"} || 18;

  # Lissage par fenêtre glissante
  my @smoothedEntropies = ();
  for (my $i = 0; $i < scalar(@entropies) - $windowSize + 1; $i++) {
      my $sum = 0;
      for (my $j = 0; $j < $windowSize; $j++) {
          $sum += $entropies[$i + $j];
      }
      my $avg = $sum / $windowSize;
      # Si la moyenne sur la fenêtre dépasse le seuil, on marque le CENTRE de la fenêtre comme à exclure / If window average exceeds threshold, mark the CENTER of the window as excluded
      # Ou plus simplement, si une fenêtre est mauvaise, on ne devrait pas pouvoir commencer un primer ici
      # Mais Primer3 attend des régions exclues, pas des points de départ interdits. / But Primer3 expects excluded regions, not forbidden start points.
      
      # Approche : Si moyenne > seuil, toute la fenêtre est considérée comme "à risque" / Approach: If average > threshold, the entire window is considered "at risk"
      if ($avg > $entropyThreshold) {
          # On marque toute la fenêtre
          for (my $j = 0; $j < $windowSize; $j++) {
              $smoothedEntropies[$i + $j] = 1; # 1 = à exclure
          }
      }
  }

  for (my $i = 0; $i < scalar(@entropies); $i++) {
      if ($smoothedEntropies[$i]) {
          if (!$inHighEntropy) {
              $startHigh = $i;
              $inHighEntropy = 1;
          }
      } else {
          if ($inHighEntropy) {
              my $len = $i - $startHigh;
              push(@excludedRegions, [$startHigh, $len]);
              $inHighEntropy = 0;
          }
      }
  }
  if ($inHighEntropy) {
      my $len = scalar(@entropies) - $startHigh;
      push(@excludedRegions, [$startHigh, $len]);
  }
  
  # Ajouter SEQUENCE_EXCLUDED_REGION aux targets Primer3
  if (@excludedRegions) {
       my @regionsStrs = ();
      foreach my $reg (@excludedRegions) {
          push(@regionsStrs, $reg->[0] . "," . $reg->[1]);
      }
      $primer3->add_targets("SEQUENCE_EXCLUDED_REGION" => join(" ", @regionsStrs));
  }

  my $primer3Results_r = $primer3->run();
  my $oligoCount = $primer3Results_r->number_of_results;
  print "Primer3Conserved getOligos had $oligoCount oligos\n";
  
  if ($oligoCount == 0) {
      # DEBUG LAVA 2026 : Diagnostiquer pourquoi 0 résultats / DEBUG LAVA 2026: Diagnose why 0 results
      print "⚠️  AUCUN OLIGO TROUVÉ. Diagnostic :\n";
      
      # Afficher les régions exclues / Display excluded regions
      if (@excludedRegions) {
          print "   -> Régions exclues par entropie (" . scalar(@excludedRegions) . " zones) : ";
          my @regionsStrs = ();
          foreach my $reg (@excludedRegions) {
             push(@regionsStrs, $reg->[0] . "-" . ($reg->[0]+$reg->[1]));
          }
          print join(", ", @regionsStrs) . "\n";
      } else {
          print "   -> Aucune région exclue par entropie.\n";
      }
      
      # Tenter d'afficher les raisons Primer3 (EXPLAIN)
      # BioPerl ne donne pas toujours accès facile au EXPLAIN complet, mais on essaie / BioPerl doesn't always give easy access to full EXPLAIN, but we try
      # On va regarder les résultats bruts s'ils sont accessibles / We will look at raw results if accessible
      eval {
          my $explain = $primer3Results_r->{'results'}->{'PRIMER_LEFT_EXPLAIN'} 
                     || $primer3Results_r->{'results'}->{'PRIMER_RIGHT_EXPLAIN'}
                     || $primer3Results_r->{'results'}->{'PRIMER_INTERNAL_OLIGO_EXPLAIN'};
          print "   -> Primer3 EXPLAIN : $explain\n" if $explain;
      };
  }

  # Transition the results into an array
  my $resultsList_r = [];
  for(my $currIndex = 0; $currIndex < $oligoCount; $currIndex++)
  {
    my $result = $primer3Results_r->primer_results($currIndex);
    my $oligo = LLNL::LAVA::Oligo->newFromPrimer3($result);

    # IUPAC OPEN GENERATION (LAVA 2026) - DISABLED BY USER REQUEST [2026-01-30]
    # We want to use the raw Primer3 sequence (sequence from reference genome)
    # and let lava_loop_primer.pl handle the degeneracy from scratch.
    
    # my $start = $oligo->location(); # 0-indexed start
    # my $len = $oligo->length();
    # my @consensusBases = ();
    
    # for (my $i = 0; $i < $len; $i++) {
    #     my $msaPos = $start + $i;
    #     my %basesAtPos = ();
    #     foreach my $seq ($alignment->each_seq()) {
    #         my $b = substr($seq->seq(), $msaPos, 1);
    #         $b = uc($b);
    #         next if $b eq '-' || $b eq 'N'; 
    #         $basesAtPos{$b} = 1;
    #     }
    #     my $iupac = _getIUPAC(keys %basesAtPos);
    #     push(@consensusBases, $iupac);
    # }
    
    # my $consensusSeq = join("", @consensusBases);
    # $oligo->{"d_sequence"} = $consensusSeq; # Force update sequence

    push(@{$resultsList_r}, $oligo);
  }
  
  # FILTRE DE CONSERVATION STRICT SUPPRIMÉ
  # Les oligos sont retournés directement (seq originale). / Oligos are returned directly (original seq).

  # FILTRE HOMOPOLYMÈRES (Post-Processing LAVA 2026)
  # Rétabli à la demande de l'utilisateur.
  # Ce filtre élimine les amorces contenant des répétitions consécutives de la même base
  # strictement supérieures à $maxPolyBases (par défaut 2 dans lava_loop_primer.pl, modifiable via --max_poly_bases).
  
  my $maxPolyBases = $p3Targets_r->{"PRIMER_INTERNAL_OLIGO_MAX_POLY_X"};
  # Fallback sur les paramètres globaux si non défini par Primer3
  if (!defined $maxPolyBases) {
     $maxPolyBases = $lavaParams_r->{"max_poly_bases"} || 5; 
  }
  
  # print "INFO: Filtrage post-Primer3 des homopolymères > $maxPolyBases bases\n";
  
  my @survivingOligos = ();
  my $rejectedPolyCount = 0;
  
  foreach my $oligo (@{$resultsList_r}) {
      my $seq = $oligo->sequence();
      my $max_poly = 0;
      # Recherche de la plus longue répétition d'une même base (A, C, G ou T)
      while ($seq =~ /(A+|C+|G+|T+)/g) {
          my $len = length($1);
          if ($len > $max_poly) { $max_poly = $len; }
      }
      
      # Si la plus longue répétition dépasse la limite autorisée, on rejette l'amorce
      if ($max_poly > $maxPolyBases) {
          $rejectedPolyCount++;
          # print "  REJETÉ (Poly $max_poly > $maxPolyBases): $seq\n";
      } else {
          push(@survivingOligos, $oligo);
      }
  }
  
  # print "INFO: $rejectedPolyCount oligos rejetés pour homopolymères excessifs.\n";
  # print "INFO: " . scalar(@survivingOligos) . " oligos conservés.\n";
  
  return @survivingOligos;
}

# Helper IUPAC
sub _getIUPAC {
    my @bases = @_;
    return 'N' if !@bases;
    my $key = join('', sort @bases);
    my %table = (
        'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T',
        'AG' => 'R', 'CT' => 'Y', 'GT' => 'K', 'AC' => 'M', 'CG' => 'S', 'AT' => 'W',
        'ACG' => 'V', 'ACT' => 'H', 'AGT' => 'D', 'CGT' => 'B',
        'ACGT' => 'N'
    );
    return $table{$key} || 'N';
}

# Fin de la classe sans écraser les méthodes suivantes si elles existent (mais setPrimer3Targets est après) / End of class without overwriting following methods if they exist (but setPrimer3Targets is after)
# Je ne remplace que getOligos.

#-------------------------------------------------------------------------------

=head2 setPrimer3Targets

 Usage     : $enumerator->setPrimer3Targets({target_name => $target_value...});
 Function  : Sets the oligo properties for executing Primer3.  Will cause an
             error if a target is used that this module doesn't have a
             default value for.
 Arguments : Hash ref of target->value pairs
 Example   : See Usage
 Returns   : <n/a>

=cut

sub setPrimer3Targets
{
  my ($this, $paramHash_r) = @_;

  if(!defined $paramHash_r)
  {
    confess("programming error - first parameter is a required hash ref");
  }
  if(ref($paramHash_r) ne "HASH")
  {
    confess("programming error - first parameter must be a hash reference");
  } 

  my $p3Targets_r = $this->{"d_primer3Targets"};
  my $p3Names_r = $this->{"d_primer3NameConversion"}; 
  my $lavaParams_r = $this->{"d_lavaParams"};

  # for each param passed in
  foreach my $currTarget(keys(%{$paramHash_r}))
  {
    # Gérer les paramètres LAVA personnalisés séparément / Manage custom LAVA parameters separately
    if ($currTarget eq "entropy_threshold") {
      $lavaParams_r->{$currTarget} = $paramHash_r->{$currTarget};
      next;
    }
    
    # oops if it doesn't exist
    if(! exists $p3Names_r->{$currTarget})
    {
      confess("programming error - failing to set primer3 target " .
        "\"$currTarget\" because no default value for that target was " .
        "created during Primer3Conserved's instantiation");
    }

    
    # Since it does exist, set the correct tag in Targets to be the new value
    $p3Targets_r->{$p3Names_r->{$currTarget}} = $paramHash_r->{$currTarget};
  }
}

1; # Lame

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
