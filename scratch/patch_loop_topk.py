import re

with open("lava_loop_primer.pl", "r") as f:
    content = f.read()

# 1. Add $msaNumSeqs before Forward B&B
content = content.replace(
    'my $_sig_fwd_pruned = 0;',
    'my $msaNumSeqs = $inputMSA->num_sequences();\n  my $_sig_fwd_pruned = 0;'
)

# Replace pruning conditions for Forward and Reverse
# We will just replace ALL occurrences of `$bestSetPenalty` with `$kthBestPenalty`
content = content.replace('$bestSetPenalty', '$kthBestPenalty')

# 2. Forward loop: bestSetPenalty to topCandidates
fwd_best_set = """          my $kthBestPenalty = 1000000;"""
fwd_top_k = """          my @topCandidates = ();
          my $kthBestPenalty = $maxSignaturePenalty;"""
content = content.replace(fwd_best_set, fwd_top_k, 1)

# 3. Replace currentSetPenalty logic in Forward
fwd_save = """                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $kthBestPenalty) {
                                          $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                                          $chunk_infos{$innerIndex} = $includeLoopPrimers ? [$loopInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo];
                                          $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                                          $kthBestPenalty = $currentSetPenalty;
                                      }"""
fwd_save_topk = """                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $kthBestPenalty) {
                                          my $innerBitVec = $masterInnerF_data_r->[$innerIndex][4];
                                          my $middleBitVec = $masterMiddleF_data_r->[$i][4];
                                          my $outerBitVec = $masterOuterF_data_r->[$k][4];
                                          my $intersection = $innerBitVec & $middleBitVec & $outerBitVec;
                                          if ($includeLoopPrimers) {
                                              my $loopBitVec = $masterLoopF_data_r->[$j][4];
                                              $intersection = $intersection & $loopBitVec if defined $loopBitVec;
                                          }
                                          my $coverage_count = unpack("%32b*", $intersection);
                                          my $coveragePct = ($msaNumSeqs > 0) ? ($coverage_count / $msaNumSeqs) * 100 : 0;
                                          
                                          my $candidate = {
                                              infos => $includeLoopPrimers ? [$loopInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo],
                                              penalties => [$spacingPenalty, $primer3Penalty, $detailStr],
                                              coverage => $coveragePct,
                                              total_penalty => $currentSetPenalty
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

# Save chunk_infos in Forward (at the end of inner loop)
fwd_end = """          } # End Loop"""
fwd_end_topk = """          if (scalar(@topCandidates) > 0) {
              $chunk_hits++;
              $chunk_infos{$innerIndex} = \\@topCandidates;
          }
          } # End Loop"""
content = content.replace(fwd_end, fwd_end_topk, 1)

# Now for Reverse
rev_best_set = """          my $kthBestPenalty = 1000000;"""
rev_top_k = """          my @topCandidates = ();
          my $kthBestPenalty = $maxSignaturePenalty;"""
content = content.replace(rev_best_set, rev_top_k, 1)

# For reverse save, the exact same string might be found (because they were identical originally!)
rev_save_topk = """                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $kthBestPenalty) {
                                          my $innerBitVec = $masterInnerR_data_r->[$innerIndex][4];
                                          my $middleBitVec = $masterMiddleR_data_r->[$i][4];
                                          my $outerBitVec = $masterOuterR_data_r->[$k][4];
                                          my $intersection = $innerBitVec & $middleBitVec & $outerBitVec;
                                          if ($includeLoopPrimers) {
                                              my $loopBitVec = $masterLoopR_data_r->[$j][4];
                                              $intersection = $intersection & $loopBitVec if defined $loopBitVec;
                                          }
                                          my $coverage_count = unpack("%32b*", $intersection);
                                          my $coveragePct = ($msaNumSeqs > 0) ? ($coverage_count / $msaNumSeqs) * 100 : 0;
                                          
                                          my $candidate = {
                                              infos => $includeLoopPrimers ? [$loopInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo],
                                              penalties => [$spacingPenalty, $primer3Penalty, $detailStr],
                                              coverage => $coveragePct,
                                              total_penalty => $currentSetPenalty
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
content = content.replace(fwd_save, rev_save_topk, 1) # Note: fwd_save is the original target string, replace the second occurrence!

rev_end = """          } # End Loop"""
rev_end_topk = """          if (scalar(@topCandidates) > 0) {
              $chunk_hits++;
              $chunk_infos{$innerIndex} = \\@topCandidates;
          }
          } # End Loop"""
content = content.replace(rev_end, rev_end_topk, 1) # Replaces the second occurrence

# Combination Phase
comb_start = """  for(my $i = 0; $i < $innerForwardCount; $i++)
  {
      my $f_innerInfo = $masterInnerF_r->[$i];
      if(defined($bestForwardInfos[$i]))
      {
          my ($f_loopInfo, $f_middleInfo, $f_outerInfo);
          if ($includeLoopPrimers) {
              ($f_loopInfo, $f_middleInfo, $f_outerInfo) = @{$bestForwardInfos[$i]};
          } else {
              ($f_middleInfo, $f_outerInfo) = @{$bestForwardInfos[$i]};
          }
          
          my ($f_spacingPenalty, $f_primer3Penalty, $f_detailStr) = @{$bestForwardPenalties[$i]};"""

comb_start_topk = """  for(my $i = 0; $i < $innerForwardCount; $i++)
  {
      my $f_innerInfo = $masterInnerF_r->[$i];
      my $f_candidates = $bestForwardInfos[$i];
      if(defined($f_candidates))
      {
        foreach my $f_cand (@$f_candidates) {
          my ($f_loopInfo, $f_middleInfo, $f_outerInfo);
          if ($includeLoopPrimers) {
              ($f_loopInfo, $f_middleInfo, $f_outerInfo) = @{$f_cand->{infos}};
          } else {
              ($f_middleInfo, $f_outerInfo) = @{$f_cand->{infos}};
          }
          
          my ($f_spacingPenalty, $f_primer3Penalty, $f_detailStr) = @{$f_cand->{penalties}};"""

content = content.replace(comb_start, comb_start_topk, 1)

comb_rev = """          for(my $j = 0; $j < $innerReverseCount; $j++)
          {
              my $r_innerInfo = $masterInnerR_r->[$j];
              if(defined($bestReverseInfos[$j]))
              {
                  my ($r_loopInfo, $r_middleInfo, $r_outerInfo);
                  if ($includeLoopPrimers) {
                      ($r_loopInfo, $r_middleInfo, $r_outerInfo) = @{$bestReverseInfos[$j]};
                  } else {
                      ($r_middleInfo, $r_outerInfo) = @{$bestReverseInfos[$j]};
                  }
                  my ($r_spacingPenalty, $r_primer3Penalty, $r_detailStr) = @{$bestReversePenalties[$j]};"""

comb_rev_topk = """          for(my $j = 0; $j < $innerReverseCount; $j++)
          {
              my $r_innerInfo = $masterInnerR_r->[$j];
              my $r_candidates = $bestReverseInfos[$j];
              if(defined($r_candidates))
              {
                foreach my $r_cand (@$r_candidates) {
                  my ($r_loopInfo, $r_middleInfo, $r_outerInfo);
                  if ($includeLoopPrimers) {
                      ($r_loopInfo, $r_middleInfo, $r_outerInfo) = @{$r_cand->{infos}};
                  } else {
                      ($r_middleInfo, $r_outerInfo) = @{$r_cand->{infos}};
                  }
                  my ($r_spacingPenalty, $r_primer3Penalty, $r_detailStr) = @{$r_cand->{penalties}};"""

content = content.replace(comb_rev, comb_rev_topk, 1)

# Add closing braces for foreach
content = re.sub(r'(my \$comboPenalty = \$spacingPenalty \+ \$primer3Penalty;.*?)(\n\s+})\n\s+}(\n\s+})\n\s+}', 
                 r'\1\n                  }\n              }\n          }\n        }\n      }\n  }', content, flags=re.DOTALL)

with open("lava_loop_primer.pl.new", "w") as f:
    f.write(content)
