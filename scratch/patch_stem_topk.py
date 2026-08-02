import re

with open("lava_stem_primer.pl", "r") as f:
    content = f.read()

# 1. Add options
options_replace = """      "signature_max_length=i" => \$options{"signature_max_length"},
      "signature_min_length=i" => \$options{"signature_min_length"},
      "total_signature_length=i" => \$options{"total_signature_length"},"""
options_new = """      "signature_max_length=i" => \$options{"signature_max_length"},
      "signature_min_length=i" => \$options{"signature_min_length"},
      "total_signature_length=i" => \$options{"total_signature_length"},
      
      "half_signature_candidates=i" => \$options{"half_signature_candidates"},
      "max_signature_penalty=i" => \$options{"max_signature_penalty"},"""
content = content.replace(options_replace, options_new, 1)

defaults_replace = """  my $totalSignatureLength = 
    optionWithDefault($options_r, "total_signature_length",
      $signatureMaxLength); # Default to max length if not specified"""
defaults_new = """  my $totalSignatureLength = 
    optionWithDefault($options_r, "total_signature_length",
      $signatureMaxLength); # Default to max length if not specified

  my $halfSignatureCandidates = optionWithDefault($options_r, "half_signature_candidates", 5);
  my $maxSignaturePenalty = optionWithDefault($options_r, "max_signature_penalty", 1000000);"""
content = content.replace(defaults_replace, defaults_new, 1)

# 2. Add msaNumSeqs
content = content.replace(
    'my $_sig_fwd_pruned = 0;',
    'my $msaNumSeqs = $inputMSA->num_sequences();\n  my $_sig_fwd_pruned = 0;'
)

# 3. Replace BestSetPenalty globally to kthBestPenalty
content = content.replace('$bestSetPenalty', '$kthBestPenalty')

fwd_best_set = """          my $kthBestPenalty = 1000000;"""
fwd_top_k = """          my @topCandidates = ();
          my $kthBestPenalty = $maxSignaturePenalty;"""
content = content.replace(fwd_best_set, fwd_top_k, 1)

# Replace forward save logic
fwd_save = """                  my $forwardSetPenalty = $spacingPenalty + $primer3Penalty;
                  if($forwardSetPenalty < $kthBestPenalty)
                  {
                    $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                    $chunk_infos{$innerIndex} = [$stemInfo, $middleInfo, $outerInfo];
                    $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                    $kthBestPenalty = $forwardSetPenalty;
                  }"""
fwd_save_topk = """                  my $forwardSetPenalty = $spacingPenalty + $primer3Penalty;
                  if($forwardSetPenalty < $kthBestPenalty)
                  {
                    my $innerBitVec = $masterInnerF_data_r->[$innerIndex][4];
                    my $middleBitVec = $masterMiddleF_data_r->[$i][4];
                    my $outerBitVec = $masterOuterF_data_r->[$k][4];
                    my $stemBitVec = $masterStemF_data_r->[$j][4];
                    my $intersection = $innerBitVec & $middleBitVec & $outerBitVec & $stemBitVec;
                    my $coverage_count = unpack("%32b*", $intersection);
                    my $coveragePct = ($msaNumSeqs > 0) ? ($coverage_count / $msaNumSeqs) * 100 : 0;
                    
                    my $candidate = {
                        infos => [$stemInfo, $middleInfo, $outerInfo],
                        penalties => [$spacingPenalty, $primer3Penalty, $detailStr],
                        coverage => $coveragePct,
                        total_penalty => $forwardSetPenalty
                    };
                    
                    push @topCandidates, $candidate;
                    @topCandidates = sort { $b->{coverage} <=> $a->{coverage} || $a->{total_penalty} <=> $b->{total_penalty} } @topCandidates;
                    if (scalar(@topCandidates) > $halfSignatureCandidates) {
                        pop @topCandidates;
                    }
                    if (scalar(@topCandidates) == $halfSignatureCandidates) {
                        $kthBestPenalty = $topCandidates[-1]->{total_penalty};
                    }
                  }"""
content = content.replace(fwd_save, fwd_save_topk, 1)

fwd_end = """            } # End forward STEM iteration"""
fwd_end_topk = """            } # End forward STEM iteration
          if (scalar(@topCandidates) > 0) {
              $chunk_hits++;
              $chunk_infos{$innerIndex} = \\@topCandidates;
          }"""
content = content.replace(fwd_end, fwd_end_topk, 1)


# Reverse loop
content = content.replace(fwd_best_set, fwd_top_k, 1)

rev_save = """                  my $reverseSetPenalty = $spacingPenalty + $primer3Penalty;
                  if($reverseSetPenalty < $kthBestPenalty)
                  {
                    $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                    $chunk_infos{$innerIndex} = [$stemInfo, $middleInfo, $outerInfo];
                    $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                    $kthBestPenalty = $reverseSetPenalty;
                  }"""
rev_save_topk = """                  my $reverseSetPenalty = $spacingPenalty + $primer3Penalty;
                  if($reverseSetPenalty < $kthBestPenalty)
                  {
                    my $innerBitVec = $masterInnerR_data_r->[$innerIndex][4];
                    my $middleBitVec = $masterMiddleR_data_r->[$i][4];
                    my $outerBitVec = $masterOuterR_data_r->[$k][4];
                    my $stemBitVec = $masterStemR_data_r->[$j][4];
                    my $intersection = $innerBitVec & $middleBitVec & $outerBitVec & $stemBitVec;
                    my $coverage_count = unpack("%32b*", $intersection);
                    my $coveragePct = ($msaNumSeqs > 0) ? ($coverage_count / $msaNumSeqs) * 100 : 0;
                    
                    my $candidate = {
                        infos => [$stemInfo, $middleInfo, $outerInfo],
                        penalties => [$spacingPenalty, $primer3Penalty, $detailStr],
                        coverage => $coveragePct,
                        total_penalty => $reverseSetPenalty
                    };
                    
                    push @topCandidates, $candidate;
                    @topCandidates = sort { $b->{coverage} <=> $a->{coverage} || $a->{total_penalty} <=> $b->{total_penalty} } @topCandidates;
                    if (scalar(@topCandidates) > $halfSignatureCandidates) {
                        pop @topCandidates;
                    }
                    if (scalar(@topCandidates) == $halfSignatureCandidates) {
                        $kthBestPenalty = $topCandidates[-1]->{total_penalty};
                    }
                  }"""
content = content.replace(rev_save, rev_save_topk, 1)

rev_end = """            } # End reverse STEM iteration"""
rev_end_topk = """            } # End reverse STEM iteration
          if (scalar(@topCandidates) > 0) {
              $chunk_hits++;
              $chunk_infos{$innerIndex} = \\@topCandidates;
          }"""
content = content.replace(rev_end, rev_end_topk, 1)


# Combination phase
comb_start = """  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      next unless defined $bestForwardInfos[$i];
      
      my $innerF = $masterInnerF_r->[$i];
      my $f_set_infos = $bestForwardInfos[$i]; # [StemF, MidF, OutF]
       
      # InnerF (F1c) Location data"""
comb_start_new = """  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      next unless defined $bestForwardInfos[$i];
      my $innerF = $masterInnerF_r->[$i];
      foreach my $f_cand (@{$bestForwardInfos[$i]}) {
          my $f_set_infos = $f_cand->{infos};
          
          # InnerF (F1c) Location data"""
content = content.replace(comb_start, comb_start_new, 1)

comb_rev = """      for(my $j = 0; $j < scalar(@{$masterInnerR_r}); $j++) {
          next unless defined $bestReverseInfos[$j];
          
          my $innerR = $masterInnerR_r->[$j];
          my $r_set_infos = $bestReverseInfos[$j]; # [StemR, MidR, OutR]"""
comb_rev_new = """      for(my $j = 0; $j < scalar(@{$masterInnerR_r}); $j++) {
          next unless defined $bestReverseInfos[$j];
          my $innerR = $masterInnerR_r->[$j];
          foreach my $r_cand (@{$bestReverseInfos[$j]}) {
              my $r_set_infos = $r_cand->{infos};"""
content = content.replace(comb_rev, comb_rev_new, 1)

comb_penalties = """          # Add total penalty tag
          my $f_penalty = $bestForwardPenalties[$i]->[0] + $bestForwardPenalties[$i]->[1];
          my $r_penalty = $bestReversePenalties[$j]->[0] + $bestReversePenalties[$j]->[1];
          $lampSignature->setTag("lamp_penalty", $f_penalty + $r_penalty);
          $lampSignature->setTag("penalty_notes", sprintf("Total F:%.1f R:%.1f | F{%s} | R{%s}", $f_penalty, $r_penalty, $bestForwardPenalties[$i]->[2], $bestReversePenalties[$j]->[2]));
          
          push(@{$allFoundSignatures_r}, $lampSignature);
          $combinedSignatureCount++;
      }
  }"""
comb_penalties_new = """          # Add total penalty tag
          my $f_penalty = $f_cand->{penalties}->[0] + $f_cand->{penalties}->[1];
          my $r_penalty = $r_cand->{penalties}->[0] + $r_cand->{penalties}->[1];
          $lampSignature->setTag("lamp_penalty", $f_penalty + $r_penalty);
          $lampSignature->setTag("penalty_notes", sprintf("Total F:%.1f R:%.1f | F{%s} | R{%s}", $f_penalty, $r_penalty, $f_cand->{penalties}->[2], $r_cand->{penalties}->[2]));
          
          push(@{$allFoundSignatures_r}, $lampSignature);
          $combinedSignatureCount++;
          }
      }
  }"""
content = content.replace(comb_penalties, comb_penalties_new, 1)

with open("lava_stem_primer.pl.new", "w") as f:
    f.write(content)
