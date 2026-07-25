import sys
import re

def process_file(filename, is_loop):
    with open(filename, 'r') as f:
        content = f.read()

    helpers = """
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
    return \\@st;
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
"""
    if "sub build_rmq" not in content:
        content += helpers

    # --- FORWARD SCAN MODIFICATIONS ---
    fwd_init_target = r'(my \$_sig_fwd_t0\s*=\s*time\(\);)'
    
    if is_loop:
        fwd_init = """
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
  my $fwd_prog_file = "$options_r->{'output_file'}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {
      for (1..30) { print $fh "0,0,0,0\\n"; }
      close $fh;
  }
  
  \\1"""
    else:
        fwd_init = """
  # --- B&B Initialization Forward ---
  my $min_val_f = sub { my $m = $_[0]; for(@_) { $m = $_ if $_ < $m } return $m; };
  my $minS_stemToMiddle_F = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_F = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $rmq_middle_f = build_rmq($masterMiddleF_data_r, 2);
  my $rmq_outer_f  = build_rmq($masterOuterF_data_r, 2);
  my $min_P_outer_F = @$masterOuterF_data_r ? query_rmq($rmq_outer_f, 0, scalar(@$masterOuterF_data_r)-1) * $outerPenaltyWeight : 0;
  
  my $_sig_fwd_pruned = 0;
  my $_sig_fwd_evaluated = 0;
  my $fwd_prog_file = "$options{'output_file'}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {
      for (1..30) { print $fh "0,0,0,0\\n"; } 
      close $fh;
  }
  
  \\1"""
  
    content = re.sub(fwd_init_target, fwd_init, content, 1)

    content = re.sub(r'my \$fwd_chunk_size = int\(\(\$innerForwardCount \+ \$num_fwd_chunks - 1\) / \$num_fwd_chunks\);\s*\$fwd_chunk_size = 1 if \$fwd_chunk_size < 1;', r'my $fwd_chunk_size = int(($innerForwardCount + $num_fwd_chunks - 1) / $num_fwd_chunks);\n  $fwd_chunk_size = 1 if $fwd_chunk_size < 1;', content, 1)
    
    finish_fwd_target = r'(\$_sig_fwd_hits \+= \$data_ref->\{hits\} \|\| 0;\s*\$_sig_fwd_done \+= \$data_ref->\{done\} \|\| 0;)\s*if \(\$_LAVA_IS_TTY \|\| 1\) \{.*?\n\s*\}'
    finish_fwd_repl = r'\1\n          $_sig_fwd_pruned += $data_ref->{pruned} || 0;\n          $_sig_fwd_evaluated += $data_ref->{evaluated} || 0;'
    content = re.sub(finish_fwd_target, finish_fwd_repl, content, flags=re.DOTALL)

    chunk_fwd_target = r'for \(my \$chunk_start = 0; \$chunk_start < \$innerForwardCount; \$chunk_start \+= \$fwd_chunk_size\) \{.*?\$pm_fwd->start\(\$chunk_start\) and next;\s*my %chunk_infos = \(\);\s*my %chunk_penalties = \(\);\s*my \$chunk_hits = 0;\s*my \$chunk_done = 0;\s*for\(my \$innerIndex = \$chunk_start; \$innerIndex <= \$chunk_end; \$innerIndex\+\+\)'
    chunk_fwd_repl = r'''for (my $chunk_id = 0; $chunk_id < $num_fwd_chunks; $chunk_id++) {
      $pm_fwd->start($chunk_id) and next;
      
      my %chunk_infos = ();
      my %chunk_penalties = ();
      my $chunk_hits = 0;
      my $chunk_done = 0;
      my $chunk_pruned = 0;
      my $chunk_evaluated = 0;
      
      for(my $innerIndex = $chunk_id; $innerIndex < $innerForwardCount; $innerIndex += $num_fwd_chunks)'''
    content = re.sub(chunk_fwd_target, chunk_fwd_repl, content, flags=re.DOTALL)

    if is_loop:
        mid_out_fwd_target = r'my \$middleCount = scalar\(\@\{\$masterMiddleF_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{.*?\} \s*# End Middle'
        bb_mid_out_fwd = """
              my $m_start = binary_search_first_ge($masterMiddleF_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleF_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_f, $m_start, $m_end) * $middlePenaltyWeight;
                  my $innerToLoopDistance = $includeLoopPrimers ? $innerLocation - ($loopLocation + 1) : 0;
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + 
                                     ($includeLoopPrimers ? $innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight : 0);
                  my $min_S_to_mid = $includeLoopPrimers ? $minS_loopToMiddle_F : $minS_innerToMiddle_F;
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_F + $min_S_to_mid + $minS_middleToOuter_F >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleF_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleF_data_r->[$j]};
                          
                          if ($includeLoopPrimers) {
                              next if ($middleLocation + $middleLength + $minPrimerSpacing > $loopLocation - $loopLength + 1);
                              next if ($middleLocation + $middleLength + $loopMinGap > $innerLocation);
                              next if (abs($loopTm - $midTm) > $maxTmDiff);
                          } else {
                              next if (abs($innerTm - $midTm) > $maxTmDiff);
                          }
                          
                          my $outerStartAt = $searchStartAt;
                          my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;
                          my $loopToMiddleDistance = $includeLoopPrimers ? ($loopLocation - $loopLength + 1) - ($middleLocation + $middleLength) : 0;
                          my $innerToMiddleDistance = $innerLocation - ($middleLocation + $middleLength);
                          
                          my $o_start = binary_search_first_ge($masterOuterF_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterF_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_f, $o_start, $o_end) * $outerPenaltyWeight;
                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($includeLoopPrimers ? $loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight 
                                                                      : $innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight);
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_F >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterF_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterF_data_r->[$k]};
                                      
                                      next if ($outerLocation + $outerLength + $minPrimerSpacing > $middleLocation);
                                      next if (abs($midTm - $outTm) > $maxTmDiff);
                                      
                                      my $middleToOuterDistance = $middleLocation - ($outerLocation + $outerLength);
                                      my $spacingPenalty = 0;
                                      my $primer3Penalty = 0;
                                      my $detailStr = "";
                                      
                                      if ($includeLoopPrimers) {
                                          $spacingPenalty = ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight) +
                                                            ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $loopPenalty * $loopPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                                                ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight),
                                                ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($loopPenalty * $loopPenaltyWeight),
                                                ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $spacingPenalty = ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
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
"""
    else:
        mid_out_fwd_target = r'my \$middleCount = scalar\(\@\{\$masterMiddleF_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{.*?\} \s*# End forward middle iteration'
        bb_mid_out_fwd = """
              my $m_start = binary_search_first_ge($masterMiddleF_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleF_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_f, $m_start, $m_end) * $middlePenaltyWeight;
                  my $innerToStemDistance = $includeStemPrimers ? $stemLocation - ($innerLocation + $innerLength) : 0;
                  if($innerToStemDistance < 0) { $innerToStemDistance = 0; }
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeStemPrimers ? $stemPenalty * $stemPenaltyWeight : 0) + 
                                     ($includeStemPrimers ? $innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight : 0);
                  my $min_S_to_mid = $includeStemPrimers ? $minS_stemToMiddle_F : 0; 
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_F + $min_S_to_mid + $minS_middleToOuter_F >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleF_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleF_data_r->[$j]};
                          
                          if ($includeStemPrimers) {
                              next if ($stemLocation + $stemLength + $minPrimerSpacing > $middleLocation);
                              next if (abs($stemTm - $midTm) > $maxTmDiff);
                          } else {
                              next if (abs($innerTm - $midTm) > $maxTmDiff);
                          }
                          
                          my $outerStartAt = $searchStartAt;
                          my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;
                          my $innerToMiddleDistance = $middleLocation - ($innerLocation + $innerLength);
                          if($innerToMiddleDistance < 0) { $innerToMiddleDistance = 0; }
                          
                          my $o_start = binary_search_first_ge($masterOuterF_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterF_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_f, $o_start, $o_end) * $outerPenaltyWeight;
                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight);
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_F >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterF_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterF_data_r->[$k]};
                                      
                                      next if ($outerLocation + $outerLength + $minPrimerSpacing > $middleLocation);
                                      next if (abs($midTm - $outTm) > $maxTmDiff);
                                      
                                      my $middleToOuterDistance = $middleLocation - ($outerLocation + $outerLength);
                                      if($middleToOuterDistance < 0) { $middleToOuterDistance = 0; }
                                      
                                      my $spacingPenalty = 0;
                                      my $primer3Penalty = 0;
                                      my $detailStr = "";
                                      
                                      if ($includeStemPrimers) {
                                          $spacingPenalty = ($innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight) +
                                                            ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $stemPenalty * $stemPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_S:%.1f I_M:%.1f M_O:%.1f] Thm[I:%.1f S:%.1f M:%.1f O:%.1f]", 
                                                ($innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight),
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($stemPenalty * $stemPenaltyWeight),
                                                ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $spacingPenalty = ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      }
                                      
                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $bestSetPenalty) {
                                          $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                                          $chunk_infos{$innerIndex} = $includeStemPrimers ? [$stemInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo];
                                          $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                                          $bestSetPenalty = $currentSetPenalty;
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
"""

    content = re.sub(mid_out_fwd_target, bb_mid_out_fwd, content, flags=re.DOTALL)

    if is_loop:
        finish_fwd_block_target = r'(\} \s*# End Loop\n)\s*(\} \s*# End Inner chunk loop\n)\s*\$pm_fwd->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*\}\);\s*\} \s*# End chunks\s*\$pm_fwd->wait_all_children\(\);'
    else:
        finish_fwd_block_target = r'(\} \s*# End forward STEM iteration\n)\s*(\} \s*# End forward inner chunk loop\n)\s*\$pm_fwd->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*\}\);\s*\} \s*# End chunks\s*\$pm_fwd->wait_all_children\(\);'
        
    finish_fwd_block_repl = r"""\1          
          # Intra-chunk progress reporting
          $chunk_done++;
          if ($chunk_done % 5 == 0 || $chunk_done == $fwd_chunk_size) {
              if (open(my $fh, '>>', $fwd_prog_file)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n";
                  close($fh);
                  
                  if (open(my $fh_read, '<', $fwd_prog_file)) {
                      my $total_done = 0;
                      my $total_hits = 0;
                      while(<$fh_read>) { 
                          chomp; 
                          next unless $_;
                          my ($d, $h) = split /,/, $_;
                          $total_done += $d;
                          $total_hits += $h;
                      }
                      close($fh_read);
                      
                      if ($_LAVA_IS_TTY || 1) {
                          my $elapsed = time() - $_sig_fwd_t0 + 0.001;
                          my $eta = ($total_done < $innerForwardCount) ? int(($innerForwardCount - $total_done) / ($total_done / $elapsed)) : 0;
                          my $rate = $total_done / $elapsed;
                          printf("[LAVA-PROGRESS] Signatures Forward|%d|%d|Sig: %d|%.1f it/s|%d\r", $total_done, $innerForwardCount, $total_hits, $rate, $eta);
                          my $old_h = select(STDOUT); $| = 1; select($old_h);
                      }
                  }
              }
          }
\2      
      $pm_fwd->finish(0, {
          infos => \%chunk_infos,
          penalties => \%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_fwd->wait_all_children();
  unlink $fwd_prog_file if -e $fwd_prog_file;
  
  if ($_sig_fwd_evaluated > 0) {
      my $pct = ($_sig_fwd_pruned / $_sig_fwd_evaluated) * 100;
      printf("  [Forward B&B] Elagage: %.2f%% (%d / %d branches evaluees)\\n", $pct, $_sig_fwd_pruned, $_sig_fwd_evaluated);
  }
"""
    content = re.sub(finish_fwd_block_target, finish_fwd_block_repl, content, flags=re.DOTALL)



    # REVERSE SCAN MODIFICATIONS
    rev_init_target = r'(my \$_sig_rev_t0\s*=\s*time\(\);)'
    
    if is_loop:
        rev_init = """
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
  my $rev_prog_file = "$options_r->{'output_file'}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {
      for (1..30) { print $fh "0,0,0,0\\n"; } 
      close $fh;
  }
  
  \\1"""
    else:
        rev_init = """
  # --- B&B Initialization Reverse ---
  my $minS_stemToMiddle_R = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_R = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $rmq_middle_r = build_rmq($masterMiddleR_data_r, 2);
  my $rmq_outer_r  = build_rmq($masterOuterR_data_r, 2);
  my $min_P_outer_R = @$masterOuterR_data_r ? query_rmq($rmq_outer_r, 0, scalar(@$masterOuterR_data_r)-1) * $outerPenaltyWeight : 0;
  
  my $_sig_rev_pruned = 0;
  my $_sig_rev_evaluated = 0;
  my $rev_prog_file = "$options{'output_file'}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {
      for (1..30) { print $fh "0,0,0,0\\n"; }
      close $fh;
  }
  
  \\1"""
    content = re.sub(rev_init_target, rev_init, content, 1)

    content = re.sub(r'\$num_rev_chunks = \$innerReverseCount if \$num_rev_chunks > \$innerReverseCount;\s*my \$rev_chunk_size', r'$num_rev_chunks = $innerReverseCount if $num_rev_chunks > $innerReverseCount;\n  $num_rev_chunks = 1 if $num_rev_chunks < 1;\n  my $rev_chunk_size', content, 1)

    finish_rev_target = r'(\$_sig_rev_hits \+= \$data_ref->\{hits\} \|\| 0;\s*\$_sig_rev_done \+= \$data_ref->\{done\} \|\| 0;)\s*if \(\$_LAVA_IS_TTY \|\| 1\) \{.*?\n\s*\}'
    finish_rev_repl = r'\1\n          $_sig_rev_pruned += $data_ref->{pruned} || 0;\n          $_sig_rev_evaluated += $data_ref->{evaluated} || 0;'
    content = re.sub(finish_rev_target, finish_rev_repl, content, flags=re.DOTALL)

    chunk_rev_target = r'for \(my \$chunk_start = 0; \$chunk_start < \$innerReverseCount; \$chunk_start \+= \$rev_chunk_size\) \{.*?\$pm_rev->start\(\$chunk_start\) and next;\s*my %chunk_infos = \(\);\s*my %chunk_penalties = \(\);\s*my \$chunk_hits = 0;\s*my \$chunk_done = 0;\s*for\(my \$innerIndex = \$chunk_start; \$innerIndex <= \$chunk_end; \$innerIndex\+\+\)'
    chunk_rev_repl = r'''for (my $chunk_id = 0; $chunk_id < $num_rev_chunks; $chunk_id++) {
      $pm_rev->start($chunk_id) and next;
      
      my %chunk_infos = ();
      my %chunk_penalties = ();
      my $chunk_hits = 0;
      my $chunk_done = 0;
      my $chunk_pruned = 0;
      my $chunk_evaluated = 0;
      
      for(my $innerIndex = $chunk_id; $innerIndex < $innerReverseCount; $innerIndex += $num_rev_chunks)'''
    content = re.sub(chunk_rev_target, chunk_rev_repl, content, flags=re.DOTALL)


    # B&B for REVERSE loop (Correct coordinate logic!)
    if is_loop:
        mid_out_rev_target = r'my \$middleCount = scalar\(\@\{\$masterMiddleR_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{.*?\} \s*# End Middle'
        
        bb_mid_out_rev = """
              my $m_start = binary_search_first_ge($masterMiddleR_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleR_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_r, $m_start, $m_end) * $middlePenaltyWeight;
                  my $innerToLoopDistance = $includeLoopPrimers ? ($loopLocation - $loopLength) - $innerLocation : 0;
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeLoopPrimers ? $loopPenalty * $loopPenaltyWeight : 0) + 
                                     ($includeLoopPrimers ? $innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight : 0);
                  my $min_S_to_mid = $includeLoopPrimers ? $minS_loopToMiddle_R : $minS_innerToMiddle_R;
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_R + $min_S_to_mid + $minS_middleToOuter_R >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleR_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleR_data_r->[$j]};
                          
                          if ($includeLoopPrimers) {
                              next if ($middleLocation - $middleLength + 1 - $minPrimerSpacing <= $loopLocation);
                              next if ($innerLocation + $loopMinGap > $middleLocation - $middleLength + 1);
                              next if (abs($loopTm - $midTm) > $maxTmDiff);
                          } else {
                              next if (abs($innerTm - $midTm) > $maxTmDiff);
                          }
                          
                          my $outerStartAt = $middleLocation + $minPrimerSpacing;
                          my $outerEndAt = $searchEndAt;
                          
                          my $loopToMiddleDistance = $includeLoopPrimers ? ($middleLocation - $middleLength + 1) - ($loopLocation + $loopLength) : 0;
                          my $innerToMiddleDistance = ($middleLocation - $middleLength) - $innerLocation;
                          
                          my $o_start = binary_search_first_ge($masterOuterR_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterR_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_r, $o_start, $o_end) * $outerPenaltyWeight;
                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($includeLoopPrimers ? $loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight 
                                                                      : $innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight);
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_R >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterR_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterR_data_r->[$k]};
                                      
                                      next if ($outerLocation - $outerLength + 1 - $minPrimerSpacing <= $middleLocation);
                                      next if (abs($midTm - $outTm) > $maxTmDiff);
                                      
                                      my $middleToOuterDistance = ($outerLocation - $outerLength) - $middleLocation;
                                      
                                      my $spacingPenalty = 0;
                                      my $primer3Penalty = 0;
                                      my $detailStr = "";
                                      
                                      if ($includeLoopPrimers) {
                                          $spacingPenalty = ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight) +
                                                            ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $loopPenalty * $loopPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_L:%.1f L_M:%.1f M_O:%.1f] Thm[I:%.1f L:%.1f M:%.1f O:%.1f]", 
                                                ($innerToLoopPenalties_r->[$innerToLoopDistance] * $innerToLoopPenaltyWeight),
                                                ($loopToMiddlePenalties_r->[$loopToMiddleDistance] * $loopToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($loopPenalty * $loopPenaltyWeight),
                                                ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $spacingPenalty = ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
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
"""
    else:
        mid_out_rev_target = r'my \$middleCount = scalar\(\@\{\$masterMiddleR_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{.*?\} \s*# End reverse middle iteration'
        
        bb_mid_out_rev = """
              my $m_start = binary_search_first_ge($masterMiddleR_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleR_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_r, $m_start, $m_end) * $middlePenaltyWeight;
                  
                  my $innerToStemDistance = $includeStemPrimers ? $innerLocation - ($stemLocation + $stemLength) : 0;
                  if($innerToStemDistance < 0) { $innerToStemDistance = 0; }
                  
                  my $base_penalty = ($innerPenalty * $innerPenaltyWeight) + 
                                     ($includeStemPrimers ? $stemPenalty * $stemPenaltyWeight : 0) + 
                                     ($includeStemPrimers ? $innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight : 0);
                  my $min_S_to_mid = $includeStemPrimers ? $minS_stemToMiddle_R : 0;
                  
                  if ($base_penalty + $min_P_mid_range + $min_P_outer_R + $min_S_to_mid + $minS_middleToOuter_R >= $bestSetPenalty) {
                      $chunk_pruned += ($m_end - $m_start + 1);
                  } else {
                      for(my $j = $m_start; $j <= $m_end; $j++)
                      {
                          my $middleInfo = $masterMiddleR_r->[$j];
                          my ($middleLocation, $middleLength, $middlePenalty, $midTm) = @{$masterMiddleR_data_r->[$j]};
                          
                          if ($includeStemPrimers) {
                              next if ($stemLocation - $stemLength + 1 - $minPrimerSpacing <= $middleLocation);
                              next if (abs($stemTm - $midTm) > $maxTmDiff);
                          } else {
                              next if (abs($innerTm - $midTm) > $maxTmDiff);
                          }
                          
                          my $outerStartAt = $middleLocation + $minPrimerSpacing;
                          my $outerEndAt = $searchEndAt;
                          
                          my $innerToMiddleDistance = ($middleLocation - $middleLength) - $innerLocation;
                          if($innerToMiddleDistance < 0) { $innerToMiddleDistance = 0; }
                          
                          my $o_start = binary_search_first_ge($masterOuterR_data_r, $outerStartAt);
                          my $o_end = binary_search_last_le($masterOuterR_data_r, $outerEndAt);
                          if ($o_start != -1 && $o_end != -1 && $o_start <= $o_end) {
                              my $min_P_out_range = query_rmq($rmq_outer_r, $o_start, $o_end) * $outerPenaltyWeight;
                              my $part_penalty = $base_penalty + ($middlePenalty * $middlePenaltyWeight) + 
                                                 ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight);
                              
                              if ($part_penalty + $min_P_out_range + $minS_middleToOuter_R >= $bestSetPenalty) {
                                  $chunk_pruned += ($o_end - $o_start + 1);
                              } else {
                                  for(my $k = $o_start; $k <= $o_end; $k++)
                                  {
                                      $chunk_evaluated++;
                                      my $outerInfo = $masterOuterR_r->[$k];
                                      my ($outerLocation, $outerLength, $outerPenalty, $outTm) = @{$masterOuterR_data_r->[$k]};
                                      
                                      next if ($outerLocation - $outerLength + 1 - $minPrimerSpacing <= $middleLocation);
                                      next if (abs($midTm - $outTm) > $maxTmDiff);
                                      
                                      my $middleToOuterDistance = ($outerLocation - $outerLength) - $middleLocation;
                                      if($middleToOuterDistance < 0) { $middleToOuterDistance = 0; }
                                      
                                      my $spacingPenalty = 0;
                                      my $primer3Penalty = 0;
                                      my $detailStr = "";
                                      
                                      if ($includeStemPrimers) {
                                          $spacingPenalty = ($innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight) +
                                                            ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $stemPenalty * $stemPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_S:%.1f I_M:%.1f M_O:%.1f] Thm[I:%.1f S:%.1f M:%.1f O:%.1f]", 
                                                ($innerToInnerPenalties_r->[$innerToStemDistance] * $innerToStemPenaltyWeight),
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($stemPenalty * $stemPenaltyWeight),
                                                ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      } else {
                                          $spacingPenalty = ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight) +
                                                            ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight);
                                          $primer3Penalty = $innerPenalty * $innerPenaltyWeight + $middlePenalty * $middlePenaltyWeight + $outerPenalty * $outerPenaltyWeight;
                                          $detailStr = sprintf("Spc[I_M:%.1f M_O:%.1f] Thm[I:%.1f M:%.1f O:%.1f]", 
                                                ($innerToMiddlePenalties_r->[$innerToMiddleDistance] * $innerToMiddlePenaltyWeight),
                                                ($middleToOuterPenalties_r->[$middleToOuterDistance] * $middleToOuterPenaltyWeight),
                                                ($innerPenalty * $innerPenaltyWeight), ($middlePenalty * $middlePenaltyWeight), ($outerPenalty * $outerPenaltyWeight));
                                      }
                                      
                                      my $currentSetPenalty = $spacingPenalty + $primer3Penalty;
                                      if ($currentSetPenalty < $bestSetPenalty) {
                                          $chunk_hits++ unless exists $chunk_infos{$innerIndex};
                                          $chunk_infos{$innerIndex} = $includeStemPrimers ? [$stemInfo, $middleInfo, $outerInfo] : [$middleInfo, $outerInfo];
                                          $chunk_penalties{$innerIndex} = [$spacingPenalty, $primer3Penalty, $detailStr];
                                          $bestSetPenalty = $currentSetPenalty;
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
"""

    content = re.sub(mid_out_rev_target, bb_mid_out_rev, content, flags=re.DOTALL)

    if is_loop:
        finish_rev_block_target = r'(\} \s*# End Loop\n)\s*(\} \s*# End Inner chunk loop\n)\s*\$pm_rev->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*\}\);\s*\} \s*# End chunks\s*\$pm_rev->wait_all_children\(\);'
    else:
        finish_rev_block_target = r'(\} \s*# End reverse STEM iteration\n)\s*(\} \s*# End reverse inner chunk loop\n)\s*\$pm_rev->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*\}\);\s*\} \s*# End chunks\s*\$pm_rev->wait_all_children\(\);'
        
    finish_rev_block_repl = r"""\1          
          # Intra-chunk progress reporting
          $chunk_done++;
          if ($chunk_done % 5 == 0 || $chunk_done == $rev_chunk_size) {
              if (open(my $fh, '>>', $rev_prog_file)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n";
                  close($fh);
                  
                  if (open(my $fh_read, '<', $rev_prog_file)) {
                      my $total_done = 0;
                      my $total_hits = 0;
                      while(<$fh_read>) { 
                          chomp; 
                          next unless $_;
                          my ($d, $h) = split /,/, $_;
                          $total_done += $d;
                          $total_hits += $h;
                      }
                      close($fh_read);
                      
                      if ($_LAVA_IS_TTY || 1) {
                          my $elapsed = time() - $_sig_rev_t0 + 0.001;
                          my $eta = ($total_done < $innerReverseCount) ? int(($innerReverseCount - $total_done) / ($total_done / $elapsed)) : 0;
                          my $rate = $total_done / $elapsed;
                          printf("[LAVA-PROGRESS] Signatures Reverse|%d|%d|Sig: %d|%.1f it/s|%d\r", $total_done, $innerReverseCount, $total_hits, $rate, $eta);
                          my $old_h = select(STDOUT); $| = 1; select($old_h);
                      }
                  }
              }
          }
\2      
      $pm_rev->finish(0, {
          infos => \%chunk_infos,
          penalties => \%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_rev->wait_all_children();
  unlink $rev_prog_file if -e $rev_prog_file;
  
  if ($_sig_rev_evaluated > 0) {
      my $pct = ($_sig_rev_pruned / $_sig_rev_evaluated) * 100;
      printf("  [Reverse B&B] Elagage: %.2f%% (%d / %d branches evaluees)\\n", $pct, $_sig_rev_pruned, $_sig_rev_evaluated);
  }
"""
    content = re.sub(finish_rev_block_target, finish_rev_block_repl, content, flags=re.DOTALL)

    with open(filename, 'w') as f:
        f.write(content)
    print(f"Successfully processed {filename}")

process_file("lava_loop_primer.pl", True)
process_file("lava_stem_primer.pl", False)
