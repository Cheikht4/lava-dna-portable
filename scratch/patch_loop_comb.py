import re

with open("lava_loop_primer.pl", "r") as f:
    content = f.read()

# Fix combination phase
comb_start = """  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      next unless defined $bestForwardInfos[$i]; # Skip if no valid F-half found
      
      my $innerF = $masterInnerF_r->[$i];
      my $f_set_infos = $bestForwardInfos[$i]; # [LoopF, MidF, OutF]
       
      # InnerF (F1c) Location data"""
comb_start_new = """  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      next unless defined $bestForwardInfos[$i];
      my $innerF = $masterInnerF_r->[$i];
      foreach my $f_cand (@{$bestForwardInfos[$i]}) {
          my $f_set_infos = $f_cand->{infos};
          
          # InnerF (F1c) Location data"""

content = content.replace(comb_start, comb_start_new)

comb_rev = """      for(my $j = 0; $j < scalar(@{$masterInnerR_r}); $j++) {
          next unless defined $bestReverseInfos[$j];
          
          my $innerR = $masterInnerR_r->[$j];
          my $r_set_infos = $bestReverseInfos[$j]; # [LoopR, MidR, OutR]"""
comb_rev_new = """      for(my $j = 0; $j < scalar(@{$masterInnerR_r}); $j++) {
          next unless defined $bestReverseInfos[$j];
          my $innerR = $masterInnerR_r->[$j];
          foreach my $r_cand (@{$bestReverseInfos[$j]}) {
              my $r_set_infos = $r_cand->{infos};"""

content = content.replace(comb_rev, comb_rev_new)

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
      }
  }"""

content = content.replace(comb_penalties, comb_penalties_new)

with open("lava_loop_primer.pl.new", "w") as f:
    f.write(content)
