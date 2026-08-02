import re
import sys

def patch_file(filepath, primer_type="stem"):
    with open(filepath, 'r') as f:
        content = f.read()
    
    target_regex = re.compile(
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
    
    match = target_regex.search(content)
    if not match:
        print(f"Could not find target in {filepath}")
        return
        
    include_var = "$includeLoopPrimers" if primer_type == "loop" else "$includeStemPrimers"
    
    replacement = f"""  # --- VALIDATION PAR SIGNATURE (Essential for correct tagging) ---
  my \$total_sigs_to_validate = scalar(@{{$allFoundSignatures_r}});
  print "Validating and calculating coverage for \$total_sigs_to_validate signatures...\\n";
  my \$validated_count = 0;

  my \$val_pm = Parallel::ForkManager->new(\$pipeline_threads);
  my %validation_results;

  \$val_pm->run_on_finish(sub {{
      my (\$pid, \$exit_code, \$ident, \$exit_signal, \$core_dump, \$data_r) = \@_;
      if (defined(\$data_r) && ref(\$data_r) eq 'ARRAY') {{
          foreach my \$res (\@\$data_r) {{
              my (\$idx, \$cov, \$status, \$final_ids_r, \$primer_cov_r) = \@\$res;
              \$validation_results{{\$idx}} = {{
                  coverage => \$cov,
                  status   => \$status,
                  final_ids => \$final_ids_r,
                  primer_cov => \$primer_cov_r
              }};
          }}
      }}
  }});

  # Chunking
  my \$val_chunk_size = POSIX::ceil(\$total_sigs_to_validate / (\$pipeline_threads * 4)); # Entrelacement
  \$val_chunk_size = 100 if \$val_chunk_size < 100;
  
  my \@val_chunks;
  for(my \$i = 0; \$i < \$total_sigs_to_validate; \$i += \$val_chunk_size) {{
      my \$end = \$i + \$val_chunk_size - 1;
      \$end = \$total_sigs_to_validate - 1 if \$end >= \$total_sigs_to_validate;
      push \@val_chunks, [\$i, \$end];
  }}
  
  # Distribuer en round-robin
  my \@val_worker_batches;
  for(my \$i = 0; \$i < \$pipeline_threads; \$i++) {{
      push \@val_worker_batches, [];
  }}
  for(my \$i = 0; \$i < scalar(\@val_chunks); \$i++) {{
      my \$worker_idx = \$i % \$pipeline_threads;
      push \@{{\$val_worker_batches[\$worker_idx]}}, \$val_chunks[\$i];
  }}

  for(my \$w = 0; \$w < \$pipeline_threads; \$w++) {{
      my \$batch = \$val_worker_batches[\$w];
      next if scalar(\@\$batch) == 0;
      
      \$val_pm->start and next;
      
      my \@results_for_worker;
      foreach my \$chunk (\@\$batch) {{
          my (\$start, \$end) = \@\$chunk;
          for(my \$idx = \$start; \$idx <= \$end; \$idx++) {{
              my \$signature = \$allFoundSignatures_r->[\$idx];
              
              my (\$final_ids_r, \$coverage, \$status) = calculateSignatureIntersection(
                  \$signature, 
                  scalar(\@sequences), 
                  \$signatureCommonTargetMinPercent,
                  {include_var},
                  "{primer_type}"
              );
              
              my \$primer_cov_r = [];
              eval {{ \$primer_cov_r = \$signature->getTag("primer_coverage_details"); }};
              
              push \@results_for_worker, [\$idx, \$coverage, \$status, \$final_ids_r, \$primer_cov_r];
          }}
      }}
      \$val_pm->finish(0, \\\@results_for_worker);
  }}
  
  \$val_pm->wait_all_children;

  # R\x65-appliquer les r\x65sultats dans le parent
  for(my \$idx = 0; \$idx < \$total_sigs_to_validate; \$idx++) {{
      my \$signature = \$allFoundSignatures_r->[\$idx];
      if (exists \$validation_results{{\$idx}}) {{
          my \$res = \$validation_results{{\$idx}};
          \$signature->setTag("signature_intersection_ids", \$res->{{final_ids}});
          \$signature->setTag("signature_coverage_percent", sprintf("%.2f", \$res->{{coverage}}));
          \$signature->setTag("signature_target_count", scalar(\@{{\$res->{{final_ids}}}}));
          \$signature->setTag("validation_status", \$res->{{status}});
          \$signature->setTag("primer_coverage_details", \$res->{{primer_cov}});
          \$validated_count++;
          
          if (\$validated_count % 1000 == 0 || \$validated_count == \$total_sigs_to_validate) {{
              print "[LAVA-PROGRESS] Validated \$validated_count / \$total_sigs_to_validate signatures...\\n";
          }}
      }}
  }}

  print "Validation complete.\\n";"""

    content = content[:match.start()] + replacement + content[match.end():]
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"Patched {filepath}")

patch_file("lava_stem_primer.pl", "stem")
