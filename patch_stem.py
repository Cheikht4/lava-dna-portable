import re
import sys

with open("lava_stem_primer.pl", "r") as f:
    content = f.read()

# Replace run_on_finish
old_run_on_finish = """  my %validation_results;

  $val_pm->run_on_finish(sub {
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_r) = @_;
      if (defined($data_r) && ref($data_r) eq 'ARRAY') {
          foreach my $res (@$data_r) {
              my ($idx, $cov, $status, $target_count) = @$res;
              $validation_results{$idx} = {
                  coverage => $cov,
                  status   => $status,
                  target_count => $target_count
              };
          }
      }
  });"""

new_run_on_finish = """  my %validation_results;

  my $verbose_val = $options_r->{"verbose_validation"} ? 1 : 0;
  my $verbose_base = $options_r->{"output_file"} . "_validation_detail";

  my $val_t0 = time();
  my $val_done = 0;
  my $val_passed = 0;
  my $val_rejected = 0;
  my %val_distribution = ("<20%"=>0, "20-40%"=>0, "40-60%"=>0, "60-80%"=>0, ">=80%"=>0);
  my $max_rejected_cov = -1;

  $val_pm->run_on_finish(sub {
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_r) = @_;
      if (defined($data_r) && ref($data_r) eq 'ARRAY') {
          foreach my $res (@$data_r) {
              my ($idx, $cov, $status, $target_count) = @$res;
              $validation_results{$idx} = {
                  coverage => $cov,
                  status   => $status,
                  target_count => $target_count
              };
              $val_done++;
              if ($status eq "VALIDEE") {
                  $val_passed++;
              } else {
                  $val_rejected++;
                  $max_rejected_cov = $cov if $cov > $max_rejected_cov;
              }
              
              if ($cov < 20) { $val_distribution{"<20%"}++; }
              elsif ($cov < 40) { $val_distribution{"20-40%"}++; }
              elsif ($cov < 60) { $val_distribution{"40-60%"}++; }
              elsif ($cov < 80) { $val_distribution{"60-80%"}++; }
              else { $val_distribution{">=80%"}++; }
          }
          
          if ($_LAVA_IS_TTY || 1) {
              my $elapsed = time() - $val_t0 + 0.001;
              my $rate = $val_done / $elapsed;
              my $eta = ($val_done < $total_sigs_to_validate) ? int(($total_sigs_to_validate - $val_done) / $rate) : 0;
              printf("[LAVA-PROGRESS] Validation|%d|%d|Valid: %d / Rejet: %d|%.1f it/s|%d\\r", 
                     $val_done, $total_sigs_to_validate, $val_passed, $val_rejected, $rate, $eta);
              my $old_h = select(STDOUT); $| = 1; select($old_h);
          }
      }
  });"""

content = content.replace(old_run_on_finish, new_run_on_finish)


old_chunk_loop = """  # Traitement par chunk avec ForkManager (un processus enfant par chunk)
  foreach my $chunk (@val_chunks) {
      $val_pm->start and next;
      
      my @results_for_chunk;
      my ($start, $end) = @$chunk;
      for(my $idx = $start; $idx <= $end; $idx++) {
          my $signature = $allFoundSignatures_r->[$idx];
          
          my ($final_ids_r, $coverage, $status) = calculateSignatureIntersection(
              $signature, 
              $inputMSA->num_sequences(), 
              $signatureCommonTargetMinPercent,
              $includeStemPrimers,
              "stem"
          );
          
          # NE TRANSMETTRE QUE DES SCALAIRES pour eviter l'explosion memoire
          push @results_for_chunk, [$idx, $coverage, $status, scalar(@$final_ids_r)];
      }
      $val_pm->finish(0, \\@results_for_chunk);
  }"""

new_chunk_loop = """  # Traitement par chunk avec ForkManager (un processus enfant par chunk)
  foreach my $chunk (@val_chunks) {
      $val_pm->start and next;
      
      my $verbose_fh;
      if ($verbose_val) {
          open($verbose_fh, "| gzip >> ${verbose_base}.$$" . ".log.gz") or warn "Cannot open verbose log";
      }
      
      my @results_for_chunk;
      my ($start, $end) = @$chunk;
      for(my $idx = $start; $idx <= $end; $idx++) {
          my $signature = $allFoundSignatures_r->[$idx];
          
          my ($final_ids_r, $coverage, $status) = calculateSignatureIntersection(
              $signature, 
              $inputMSA->num_sequences(), 
              $signatureCommonTargetMinPercent,
              $includeStemPrimers,
              "stem",
              $verbose_val,
              $verbose_fh
          );
          
          # NE TRANSMETTRE QUE DES SCALAIRES pour eviter l'explosion memoire
          push @results_for_chunk, [$idx, $coverage, $status, scalar(@$final_ids_r)];
      }
      if ($verbose_val && defined $verbose_fh) {
          close($verbose_fh);
      }
      $val_pm->finish(0, \\@results_for_chunk);
  }"""
content = content.replace(old_chunk_loop, new_chunk_loop)

old_parent_loop = """  $val_pm->wait_all_children;

  # Re-appliquer les resultats scalaires dans le parent
  for(my $idx = 0; $idx < $total_sigs_to_validate; $idx++) {
      my $signature = $allFoundSignatures_r->[$idx];
      if (exists $validation_results{$idx}) {
          my $res = $validation_results{$idx};
          $signature->setTag("signature_coverage_percent", sprintf("%.2f", $res->{coverage}));
          $signature->setTag("validation_status", $res->{status});
          $signature->setTag("signature_target_count", $res->{target_count});
          $validated_count++;
          
          if ($validated_count % 1000 == 0 || $validated_count == $total_sigs_to_validate) {
              print "[LAVA-PROGRESS] Validated $validated_count / $total_sigs_to_validate signatures...\\n";
          }
      }
  }"""

new_parent_loop = """  $val_pm->wait_all_children;
  print "\\n"; # Clear the progress bar line
  
  if ($verbose_val) {
      print "Aggregating verbose logs...\\n";
      system("cat ${verbose_base}.*.log.gz > ${verbose_base}.log.gz 2>/dev/null");
      system("rm -f ${verbose_base}.*.log.gz");
      print "Verbose log saved to ${verbose_base}.log.gz\\n";
  }

  print "============================================================\\n";
  print "RESUME STATISTIQUE DE LA VALIDATION\\n";
  print "============================================================\\n";
  my $pct_val = $total_sigs_to_validate > 0 ? ($val_passed / $total_sigs_to_validate * 100) : 0;
  my $pct_rej = $total_sigs_to_validate > 0 ? ($val_rejected / $total_sigs_to_validate * 100) : 0;
  printf("Total evalue : %d\\n", $total_sigs_to_validate);
  printf("Validees     : %d (%.1f%%)\\n", $val_passed, $pct_val);
  printf("Rejetees     : %d (%.1f%%)\\n", $val_rejected, $pct_rej);
  printf("Seuil requis : %.1f%%\\n", $signatureCommonTargetMinPercent);
  if ($val_rejected > 0) {
      my $ecart = $signatureCommonTargetMinPercent - $max_rejected_cov;
      printf("Couverture MAX (rejetees) : %.2f%% (ecart au seuil : -%.2f%%)\\n", $max_rejected_cov, $ecart);
  }
  print "Distribution des couvertures :\\n";
  foreach my $bin ("<20%", "20-40%", "40-60%", "60-80%", ">=80%") {
      printf("  %s : %d\\n", $bin, $val_distribution{$bin});
  }
  print "============================================================\\n";

  # Re-appliquer les resultats scalaires dans le parent
  for(my $idx = 0; $idx < $total_sigs_to_validate; $idx++) {
      my $signature = $allFoundSignatures_r->[$idx];
      if (exists $validation_results{$idx}) {
          my $res = $validation_results{$idx};
          $signature->setTag("signature_coverage_percent", sprintf("%.2f", $res->{coverage}));
          $signature->setTag("validation_status", $res->{status});
          $signature->setTag("signature_target_count", $res->{target_count});
      }
  }"""
content = content.replace(old_parent_loop, new_parent_loop)


with open("lava_stem_primer.pl", "w") as f:
    f.write(content)
