import re

for fname in ['lava_loop_primer.pl', 'lava_stem_primer.pl', 'apply_optimizations.py']:
    with open(fname, 'r') as f:
        c = f.read()

    # 1. Init FWD prog
    c = re.sub(
        r'my \$fwd_prog_file\s*=\s*"\$(?:options|options_r->)\{.*?_fwd_prog_\$\$\.txt";\s*if\s*\(open\s*my\s*\$fh,\s*\'>\',\s*\$fwd_prog_file\)\s*\{\s*for\s*\(1\.\.30\)\s*\{\s*print\s*\$fh\s*"0,0,0,0\\\\n";\s*\}\s*close\s*\$fh;\s*\}',
        r'my $fwd_prog_dir = "$options{\'output_file\'}_fwd_prog_$$";\n  $fwd_prog_dir = "$options_r->{\'output_file\'}_fwd_prog_$$" if ref($options_r);\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);',
        c
    )
    
    # 2. Init REV prog
    c = re.sub(
        r'my \$rev_prog_file\s*=\s*"\$(?:options|options_r->)\{.*?_rev_prog_\$\$\.txt";\s*if\s*\(open\s*my\s*\$fh,\s*\'>\',\s*\$rev_prog_file\)\s*\{\s*for\s*\(1\.\.30\)\s*\{\s*print\s*\$fh\s*"0,0,0,0\\\\n";\s*\}\s*close\s*\$fh;\s*\}',
        r'my $rev_prog_dir = "$options{\'output_file\'}_rev_prog_$$";\n  $rev_prog_dir = "$options_r->{\'output_file\'}_rev_prog_$$" if ref($options_r);\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);',
        c
    )

    # 3. FWD chunk report
    fwd_report_target = r'\$chunk_done\+\+;\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$fwd_chunk_size\) \{\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(open\(my \$fh, \'>>\', \$fwd_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$fwd_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*(?:my \$total_pruned = 0;\s*my \$total_eval = 0;\s*)?while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d, \$h\) = split /,/, \$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*\}\s*close\(\$fh_read\);'
    fwd_report_repl = r"""if ($chunk_done % 5 == 0) {
              # --- Intra-chunk progress reporting ---
              my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
              if (open(my $fh, '>', $prog_file_me)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n";
                  close($fh);
              }
              
              my $total_done = 0;
              my $total_hits = 0;
              foreach my $f (glob("$fwd_prog_dir/chunk_*.prog")) {
                  if (open(my $r, '<', $f)) {
                      my $line = <$r>; close($r);
                      next unless defined $line;
                      chomp $line;
                      my ($d,$h) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                  }
              }"""
    c = re.sub(fwd_report_target, fwd_report_repl, c)

    # 4. REV chunk report
    rev_report_target = r'\$chunk_done\+\+;\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$rev_chunk_size\) \{\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(open\(my \$fh, \'>>\', \$rev_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$rev_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*(?:my \$total_pruned = 0;\s*my \$total_eval = 0;\s*)?while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d, \$h\) = split /,/, \$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*\}\s*close\(\$fh_read\);'
    rev_report_repl = fwd_report_repl.replace('fwd_prog_dir', 'rev_prog_dir')
    c = re.sub(rev_report_target, rev_report_repl, c)

    # 5. Final flush FWD
    fwd_finish_target = r'(\$pm_fwd->finish\(0, \{)'
    fwd_finish_repl = r'my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1'
    c = re.sub(fwd_finish_target, fwd_finish_repl, c)

    # 6. Final flush REV
    rev_finish_target = r'(\$pm_rev->finish\(0, \{)'
    rev_finish_repl = r'my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1'
    c = re.sub(rev_finish_target, rev_finish_repl, c)

    # 7. Unlink FWD -> remove_tree FWD
    fwd_unlink_target = r'unlink \$fwd_prog_file if -e \$fwd_prog_file;'
    fwd_unlink_repl = r'use File::Path qw(remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;'
    c = re.sub(fwd_unlink_target, fwd_unlink_repl, c)

    # 8. Unlink REV -> remove_tree REV
    rev_unlink_target = r'unlink \$rev_prog_file if -e \$rev_prog_file;'
    rev_unlink_repl = r'use File::Path qw(remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;'
    c = re.sub(rev_unlink_target, rev_unlink_repl, c)

    # Add pruning display if it is NOT there
    if 'Elagage' not in c:
        fwd_elagage = r'\n  if ($_sig_fwd_evaluated > 0) {\n      my $pct = ($_sig_fwd_pruned / $_sig_fwd_evaluated) * 100;\n      printf("  [Forward B&B] Elagage: %.2f%% (%d / %d branches evaluees)\\n", $pct, $_sig_fwd_pruned, $_sig_fwd_evaluated);\n  }\n'
        c = re.sub(r'(\$pm_fwd->wait_all_children\(\);\s*use File::Path qw\(remove_tree\);\s*remove_tree\(\$fwd_prog_dir\) if -d \$fwd_prog_dir;)', r'\1' + fwd_elagage, c)

        rev_elagage = r'\n  if ($_sig_rev_evaluated > 0) {\n      my $pct = ($_sig_rev_pruned / $_sig_rev_evaluated) * 100;\n      printf("  [Reverse B&B] Elagage: %.2f%% (%d / %d branches evaluees)\\n", $pct, $_sig_rev_pruned, $_sig_rev_evaluated);\n  }\n'
        c = re.sub(r'(\$pm_rev->wait_all_children\(\);\s*use File::Path qw\(remove_tree\);\s*remove_tree\(\$rev_prog_dir\) if -d \$rev_prog_dir;)', r'\1' + rev_elagage, c)

    with open(fname, 'w') as f:
        f.write(c)
