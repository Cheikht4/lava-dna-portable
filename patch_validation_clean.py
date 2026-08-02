import re
import sys

def patch_file(filepath, primer_type="loop"):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Target regex for loop
    target_regex_loop = re.compile(
        r'  # --- VALIDATION STEP \(Essential for correct tagging\) ---\n'
        r'  print "Validating and calculating coverage for " \. scalar\(@\{\$allFoundSignatures_r\}\) \. " signatures\.\.\.\\n";\n'
        r'  my \$validated_count = 0;\n'
        r'  foreach my \$signature \(@\{\$allFoundSignatures_r\}\) \{\n'
        r'      # Calculate Intersection and Coverage\n'
        r'      # Calculer l\'intersection des sequences compatibles \(PipelineUtils unifie\)\n'
        r'      # Calculate compatible sequence intersection \(unified PipelineUtils\)\n'
        r'      my \(\$amplified_seqs, \$coverage, \$status\) = calculateSignatureIntersection\(\n'
        r'          \$signature, \n'
        r'          \$inputMSA->num_sequences\(\), \n'
        r'          \$signatureCommonTargetMinPercent,\n'
        r'          \$includeLoopPrimers,\n'
        r'          "loop"\n'
        r'      \);\n'
        r'      \$validated_count\+\+;\n'
        r'  \}\n'
        r'  print "Validation complete\.\\n";',
        re.MULTILINE
    )
    
    # Target regex for stem
    target_regex_stem = re.compile(
        r'  # --- VALIDATION PAR SIGNATURE \(Essential for correct tagging\) ---\n'
        r'  # Validation individuelle de chaque signature avant reduction\n'
        r'  # Individual per-signature validation before reduction \(mirrors LOOP behaviour\)\n'
        r'  print "Validating and calculating coverage for " \. scalar\(@\{\$allFoundSignatures_r\}\) \. " signatures\.\.\.\\n";\n'
        r'  my \$validated_count = 0;\n'
        r'  foreach my \$signature \(@\{\$allFoundSignatures_r\}\) \{\n'
        r'      # Recalcule l\'intersection et met a jour les tags de couverture\n'
        r'      # Recalculate intersection and update coverage tags\n'
        r'      my \(\$amplified_seqs, \$coverage, \$status\) = calculateSignatureIntersection\(\n'
        r'          \$signature,\n'
        r'          scalar\(@sequences\),\n'
        r'          \$signatureCommonTargetMinPercent,\n'
        r'          \$includeStemPrimers,\n'
        r'          "stem"\n'
        r'      \);\n'
        r'      \$validated_count\+\+;\n'
        r'  \}\n'
        r'  print "Validation complete\.\\n";',
        re.MULTILINE
    )

    if primer_type == "loop":
        match = target_regex_loop.search(content)
        include_var = "$includeLoopPrimers"
        total_seq = "$inputMSA->num_sequences()"
        comment_title = "  # --- VALIDATION STEP (Essential for correct tagging) ---"
    else:
        match = target_regex_stem.search(content)
        include_var = "$includeStemPrimers"
        total_seq = "scalar(@sequences)"
        comment_title = "  # --- VALIDATION PAR SIGNATURE (Essential for correct tagging) ---\n  # Validation individuelle de chaque signature avant reduction\n  # Individual per-signature validation before reduction (mirrors LOOP behaviour)"
        
    if not match:
        print(f"Could not find target in {filepath}")
        return
        
    replacement = f"""{comment_title}
  my $total_sigs_to_validate = scalar(@{{$allFoundSignatures_r}});
  print "Validating and calculating coverage for $total_sigs_to_validate signatures...\\n";
  my $validated_count = 0;

  my $val_pm = LLNL::LAVA::ForkManager->new($options_r->{{"threads"}});
  my $actual_threads = $val_pm->{{max_processes}};
  my %validation_results;

  $val_pm->run_on_finish(sub {{
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_r) = @_;
      if (defined($data_r) && ref($data_r) eq 'ARRAY') {{
          foreach my $res (@$data_r) {{
              my ($idx, $cov, $status, $final_ids_r, $primer_cov_r) = @$res;
              $validation_results{{$idx}} = {{
                  coverage => $cov,
                  status   => $status,
                  final_ids => $final_ids_r,
                  primer_cov => $primer_cov_r
              }};
          }}
      }}
  }});

  # Chunking
  my $val_chunk_size = POSIX::ceil($total_sigs_to_validate / ($actual_threads * 4)); # Entrelacement
  $val_chunk_size = 100 if $val_chunk_size < 100;
  
  my @val_chunks;
  for(my $i = 0; $i < $total_sigs_to_validate; $i += $val_chunk_size) {{
      my $end = $i + $val_chunk_size - 1;
      $end = $total_sigs_to_validate - 1 if $end >= $total_sigs_to_validate;
      push @val_chunks, [$i, $end];
  }}
  
  # Distribuer en round-robin
  my @val_worker_batches;
  for(my $i = 0; $i < $actual_threads; $i++) {{
      push @val_worker_batches, [];
  }}
  for(my $i = 0; $i < scalar(@val_chunks); $i++) {{
      my $worker_idx = $i % $actual_threads;
      push @{{$val_worker_batches[$worker_idx]}}, $val_chunks[$i];
  }}

  for(my $w = 0; $w < $actual_threads; $w++) {{
      my $batch = $val_worker_batches[$w];
      next if scalar(@$batch) == 0;
      
      $val_pm->start and next;
      
      my @results_for_worker;
      foreach my $chunk (@$batch) {{
          my ($start, $end) = @$chunk;
          for(my $idx = $start; $idx <= $end; $idx++) {{
              my $signature = $allFoundSignatures_r->[$idx];
              
              my ($final_ids_r, $coverage, $status) = calculateSignatureIntersection(
                  $signature, 
                  {total_seq}, 
                  $signatureCommonTargetMinPercent,
                  {include_var},
                  "{primer_type}"
              );
              
              my $primer_cov_r = [];
              eval {{ $primer_cov_r = $signature->getTag("primer_coverage_details"); }};
              
              push @results_for_worker, [$idx, $coverage, $status, $final_ids_r, $primer_cov_r];
          }}
      }}
      $val_pm->finish(0, \\@results_for_worker);
  }}
  
  $val_pm->wait_all_children;

  # Ré-appliquer les résultats dans le parent
  for(my $idx = 0; $idx < $total_sigs_to_validate; $idx++) {{
      my $signature = $allFoundSignatures_r->[$idx];
      if (exists $validation_results{{$idx}}) {{
          my $res = $validation_results{{$idx}};
          $signature->setTag("signature_intersection_ids", $res->{{final_ids}});
          $signature->setTag("signature_coverage_percent", sprintf("%.2f", $res->{{coverage}}));
          $signature->setTag("signature_target_count", scalar(@{{$res->{{final_ids}}}}));
          $signature->setTag("validation_status", $res->{{status}});
          $signature->setTag("primer_coverage_details", $res->{{primer_cov}});
          $validated_count++;
          
          if ($validated_count % 1000 == 0 || $validated_count == $total_sigs_to_validate) {{
              print "[LAVA-PROGRESS] Validated $validated_count / $total_sigs_to_validate signatures...\\n";
          }}
      }}
  }}

  print "Validation complete.\\n";"""

    # Escape $ and @ manually before writing!
    # Wait, in f-strings we can just not escape them if we don't use f-strings for variables, but here we used {total_seq} etc.
    # The actual string is generated correctly, no backslashes needed. We want literal $ and @ in the perl script!
    # Let's verify.
    
    content = content[:match.start()] + replacement + content[match.end():]
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"Patched {filepath}")

patch_file("lava_loop_primer.pl", "loop")
patch_file("lava_stem_primer.pl", "stem")
