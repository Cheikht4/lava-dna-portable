#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;

sub process_file {
    my ($file) = @_;
    my $content = read_file($file);

    # Add RMQ initialization before Forward loops
    my $rmq_init = <<'EOF';
  # --- B&B Initialization ---
  my $min_val_f = sub { my $m = $_[0]; for(@_) { $m = $_ if $_ < $m } return $m; };
  my $minS_innerToLoop_F = @$innerToLoopPenalties_r ? $min_val_f->(@$innerToLoopPenalties_r) * $innerToLoopPenaltyWeight : 0;
  my $minS_loopToMiddle_F = @$loopToMiddlePenalties_r ? $min_val_f->(@$loopToMiddlePenalties_r) * $loopToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_F = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $minS_innerToMiddle_F = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $rmq_middle_f = build_rmq($masterMiddleF_data_r, 2);
  my $rmq_outer_f  = build_rmq($masterOuterF_data_r, 2);
  my $rmq_loop_f   = build_rmq($masterLoopF_data_r, 2);
  my $min_P_middle_F = @$masterMiddleF_data_r ? query_rmq($rmq_middle_f, 0, scalar(@$masterMiddleF_data_r)-1) * $middlePenaltyWeight : 0;
  my $min_P_outer_F  = @$masterOuterF_data_r ? query_rmq($rmq_outer_f, 0, scalar(@$masterOuterF_data_r)-1) * $outerPenaltyWeight : 0;
  
EOF
    $content =~ s/(my \$_sig_fwd_t0   = time\(\);)/$rmq_init$1/s;

    # Same for Reverse
    my $rmq_init_rev = <<'EOF';
  # --- B&B Initialization Reverse ---
  my $minS_innerToLoop_R = @$innerToLoopPenalties_r ? $min_val_f->(@$innerToLoopPenalties_r) * $innerToLoopPenaltyWeight : 0;
  my $minS_loopToMiddle_R = @$loopToMiddlePenalties_r ? $min_val_f->(@$loopToMiddlePenalties_r) * $loopToMiddlePenaltyWeight : 0;
  my $minS_middleToOuter_R = @$middleToOuterPenalties_r ? $min_val_f->(@$middleToOuterPenalties_r) * $middleToOuterPenaltyWeight : 0;
  my $minS_innerToMiddle_R = @$innerToMiddlePenalties_r ? $min_val_f->(@$innerToMiddlePenalties_r) * $innerToMiddlePenaltyWeight : 0;
  my $rmq_middle_r = build_rmq($masterMiddleR_data_r, 2);
  my $rmq_outer_r  = build_rmq($masterOuterR_data_r, 2);
  my $rmq_loop_r   = build_rmq($masterLoopR_data_r, 2);
  my $min_P_middle_R = @$masterMiddleR_data_r ? query_rmq($rmq_middle_r, 0, scalar(@$masterMiddleR_data_r)-1) * $middlePenaltyWeight : 0;
  my $min_P_outer_R  = @$masterOuterR_data_r ? query_rmq($rmq_outer_r, 0, scalar(@$masterOuterR_data_r)-1) * $outerPenaltyWeight : 0;
  
EOF
    $content =~ s/(my \$_sig_rev_t0   = time\(\);)/$rmq_init_rev$1/s;

    # Forward Middle Loop replacement
    my $fwd_mid = <<'EOF';
              my $m_start = binary_search_first_ge($masterMiddleF_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleF_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_f, $m_start, $m_end) * $middlePenaltyWeight;
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
EOF
    $content =~ s/my \$middleCount = scalar\(\@\{\$masterMiddleF_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{\s*my \$middleInfo = \$masterMiddleF_r->\[\$j\];\s*my \(\$middleLocation, \$middleLength, \$middlePenalty, \$midTm\) = \@\{\$masterMiddleF_data_r->\[\$j\]\};\s*# Fast-Fail\s*next if \(\$middleLocation < \$middleStartAt\);\s*last if \(\$middleLocation > \$middleEndAt\);/$fwd_mid/g;

    # Forward Outer Loop replacement
    my $fwd_out = <<'EOF';
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
EOF
    $content =~ s/my \$outerCount = scalar\(\@\{\$masterOuterF_r\}\);\s*for\(my \$k = 0; \$k < \$outerCount; \$k\+\+\)\s*\{\s*my \$outerInfo = \$masterOuterF_r->\[\$k\];\s*my \(\$outerLocation, \$outerLength, \$outerPenalty, \$outTm\) = \@\{\$masterOuterF_data_r->\[\$k\]\};\s*# Fast-Fail\s*next if \(\$outerLocation < \$outerStartAt\);\s*last if \(\$outerLocation > \$outerEndAt\);/$fwd_out/g;

    # Reverse Middle Loop replacement
    my $rev_mid = <<'EOF';
              my $m_start = binary_search_first_ge($masterMiddleR_data_r, $middleStartAt);
              my $m_end = binary_search_last_le($masterMiddleR_data_r, $middleEndAt);
              if ($m_start != -1 && $m_end != -1 && $m_start <= $m_end) {
                  my $min_P_mid_range = query_rmq($rmq_middle_r, $m_start, $m_end) * $middlePenaltyWeight;
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
EOF
    $content =~ s/my \$middleCount = scalar\(\@\{\$masterMiddleR_r\}\);\s*for\(my \$j = 0; \$j < \$middleCount; \$j\+\+\)\s*\{\s*my \$middleInfo = \$masterMiddleR_r->\[\$j\];\s*my \(\$middleLocation, \$middleLength, \$middlePenalty, \$midTm\) = \@\{\$masterMiddleR_data_r->\[\$j\]\};\s*# Fast-Fail\s*next if \(\$middleLocation < \$middleStartAt\);\s*last if \(\$middleLocation > \$middleEndAt\);/$rev_mid/g;

    # Reverse Outer Loop replacement
    my $rev_out = <<'EOF';
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
EOF
    $content =~ s/my \$outerCount = scalar\(\@\{\$masterOuterR_r\}\);\s*for\(my \$k = 0; \$k < \$outerCount; \$k\+\+\)\s*\{\s*my \$outerInfo = \$masterOuterR_r->\[\$k\];\s*my \(\$outerLocation, \$outerLength, \$outerPenalty, \$outTm\) = \@\{\$masterOuterR_data_r->\[\$k\]\};\s*# Fast-Fail\s*next if \(\$outerLocation < \$outerStartAt\);\s*last if \(\$outerLocation > \$outerEndAt\);/$rev_out/g;

    # We need to close the `if ($o_start != -1...)` blocks!
    # They are inserted right where the `for` loops were. 
    # The `for` loop ended with `} # End Outer` and `} # End Middle`.
    # Let's replace the `} # End Outer` with `} } # End Outer` to close the `if`.
    $content =~ s/(\} \s*# End (forward|reverse) outer iteration)/\} $1/g;
    $content =~ s/(\} \s*# End (forward|reverse) middle iteration)/\} $1/g;
    # For loop script it's `} # End Outer`
    $content =~ s/(\} \s*# End Outer)/\} $1/g;
    $content =~ s/(\} \s*# End Middle)/\} $1/g;

    write_file($file, $content);
    print "Processed $file\n";
}

process_file("lava_loop_primer.pl");
process_file("lava_stem_primer.pl");

