import re
import sys

with open("lava_loop_primer.pl", "r") as f:
    content = f.read()

old_comb = """  #-----------------------------------------------------------------------------
  # 5. COMBINE HALVES & CREATE SIGNATURES
  #-----------------------------------------------------------------------------
  print "Combining Best F/R Halves to create LAMP Signatures...\\n";
  
  my $combinedSignatureCount = 0;

  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {"""

new_comb = """  #-----------------------------------------------------------------------------
  # 5. COMBINE HALVES & CREATE SIGNATURES
  #-----------------------------------------------------------------------------
  print "Combining Best F/R Halves to create LAMP Signatures...\\n";
  
  my $combinedSignatureCount = 0;
  
  my $combine_total = scalar(@{$masterInnerF_r});
  my $combine_done = 0;
  my $combine_t0 = time();

  for(my $i = 0; $i < scalar(@{$masterInnerF_r}); $i++) {
      $combine_done++;
      if (($_LAVA_IS_TTY || 1) && ($combine_done % 100 == 0 || $combine_done == $combine_total)) {
          my $elapsed = time() - $combine_t0 + 0.001;
          my $rate = $combine_done / $elapsed;
          my $eta = ($combine_done < $combine_total) ? int(($combine_total - $combine_done) / $rate) : 0;
          printf("[LAVA-PROGRESS] Combinaison|%d|%d|Sigs: %d|%.1f it/s|%d\\r", 
                 $combine_done, $combine_total, $combinedSignatureCount, $rate, $eta);
          my $old_h = select(STDOUT); $| = 1; select($old_h);
      }
"""

old_comb_end = """  print "Created $combinedSignatureCount complete LAMP signatures.\\n";"""

new_comb_end = """  print "\\n"; # Clear the progress bar line
  print "Created $combinedSignatureCount complete LAMP signatures.\\n";"""

if old_comb in content:
    content = content.replace(old_comb, new_comb)
    content = content.replace(old_comb_end, new_comb_end)
    with open("lava_loop_primer.pl", "w") as f:
        f.write(content)
    print("Patched loop")
else:
    print("Could not find loop target block")

with open("lava_stem_primer.pl", "r") as f:
    content2 = f.read()

old_stem_comb = """    # Now, try to combine forward and reverse primer sets into full signatures
    print "Combining Best F/R Halves to create LAMP Signatures...\\n";
    my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner iteration
    for(my $i = 0; $i < $innerForwardCount; $i++)
    {"""

new_stem_comb = """    # Now, try to combine forward and reverse primer sets into full signatures
    print "Combining Best F/R Halves to create LAMP Signatures...\\n";
    my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner iteration
    
    my $combine_total = $innerForwardCount;
    my $combine_done = 0;
    my $combine_t0 = time();

    for(my $i = 0; $i < $innerForwardCount; $i++)
    {
      $combine_done++;
      if (($_LAVA_IS_TTY || 1) && ($combine_done % 100 == 0 || $combine_done == $combine_total)) {
          my $elapsed = time() - $combine_t0 + 0.001;
          my $rate = $combine_done / $elapsed;
          my $eta = ($combine_done < $combine_total) ? int(($combine_total - $combine_done) / $rate) : 0;
          my $current_sigs = scalar(@{$allFoundSignatures_r});
          printf("[LAVA-PROGRESS] Combinaison|%d|%d|Sigs: %d|%.1f it/s|%d\\r", 
                 $combine_done, $combine_total, $current_sigs, $rate, $eta);
          my $old_h = select(STDOUT); $| = 1; select($old_h);
      }
"""

old_stem_comb_end = """  print "Found " .
    scalar(@{$allFoundSignatures_r}) .
    " total signatures across all iterations\\n";"""

new_stem_comb_end = """  print "\\n"; # Clear the progress bar line
  print "Found " .
    scalar(@{$allFoundSignatures_r}) .
    " total signatures across all iterations\\n";"""

if old_stem_comb in content2:
    content2 = content2.replace(old_stem_comb, new_stem_comb)
    content2 = content2.replace(old_stem_comb_end, new_stem_comb_end)
    with open("lava_stem_primer.pl", "w") as f:
        f.write(content2)
    print("Patched stem")
else:
    print("Could not find stem target block")
