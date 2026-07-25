import re
import os

def fix_all():
    for fname in ['lava_loop_primer.pl', 'lava_stem_primer.pl', 'apply_optimizations.py']:
        with open(fname, 'r') as f:
            c = f.read()

        # 1. FWD INIT
        init_fwd_tgt = r'my \$fwd_prog_file = "\$?[a-zA-Z_>{}\']+_fwd_prog_\$\$\.txt";\n\s*if \(open my \$fh, \'>\', \$fwd_prog_file\) \{\n\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\n"; \}\s*close \$fh;\n\s*\}'
        init_fwd_repl = r'my $fwd_prog_dir = "$options{\'output_file\'}_fwd_prog_$$";\n  $fwd_prog_dir = "$options_r->{\'output_file\'}_fwd_prog_$$" if ref($options_r);\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);'
        c = re.sub(init_fwd_tgt, init_fwd_repl, c)

        # 2. REV INIT
        init_rev_tgt = r'my \$rev_prog_file = "\$?[a-zA-Z_>{}\']+_rev_prog_\$\$\.txt";\n\s*if \(open my \$fh, \'>\', \$rev_prog_file\) \{\n\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\n"; \}\s*close \$fh;\n\s*\}'
        init_rev_repl = r'my $rev_prog_dir = "$options{\'output_file\'}_rev_prog_$$";\n  $rev_prog_dir = "$options_r->{\'output_file\'}_rev_prog_$$" if ref($options_r);\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);'
        c = re.sub(init_rev_tgt, init_rev_repl, c)

        # 3. FWD PROGRESS BLOCK
        fwd_prog_tgt = r'# Intra-chunk progress reporting\n.*?\} # End Inner chunk loop'
        fwd_prog_repl = r"""# Intra-chunk progress reporting
          if ($chunk_done % 5 == 0) {
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
                      my ($d, $h) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                  }
              }
                      
              if ($_LAVA_IS_TTY || 1) {
                  my $elapsed = time() - $_sig_fwd_t0 + 0.001;
                  my $eta = ($total_done < $innerForwardCount) ? int(($innerForwardCount - $total_done) / ($total_done / $elapsed)) : 0;
                  my $rate = $total_done / $elapsed;
                  printf("[LAVA-PROGRESS] Signatures Forward|%d|%d|Sig: %d|%.1f it/s|%d\\r", $total_done, $innerForwardCount, $total_hits, $rate, $eta);
                  my $old_h = select(STDOUT); $| = 1; select($old_h);
              }
          }
} # End Inner chunk loop"""
        
        # 4. REV PROGRESS BLOCK
        rev_prog_tgt = r"""# Intra-chunk progress reporting
          if \(\$chunk_done \% 5 == 0\) \{
              my \$prog_file_me = "\$rev_prog_dir/chunk_\$chunk_id\.prog";"""
        rev_prog_repl = fwd_prog_repl.replace('fwd_prog_dir', 'rev_prog_dir').replace('fwd_chunk_size', 'rev_chunk_size').replace('innerForwardCount', 'innerReverseCount').replace('Forward', 'Reverse').replace('_fwd_', '_rev_')
        
        # First, we need to match and replace the FWD block.
        # But `.*?` might match TOO much if we just run it on the whole file. 
        # Better strategy: Find chunks explicitly.
        
        # Let's split by `# Intra-chunk progress reporting`
        parts = c.split('# Intra-chunk progress reporting\n')
        if len(parts) == 3:
            # We have FWD and REV
            # Fix FWD
            fwd_tail = parts[1].split('} # End Inner chunk loop\n', 1)
            parts[1] = fwd_prog_repl.replace('# Intra-chunk progress reporting\n', '') + '\n' + fwd_tail[1]
            
            # Fix REV
            rev_tail = parts[2].split('} # End Inner chunk loop\n', 1)
            parts[2] = rev_prog_repl.replace('# Intra-chunk progress reporting\n', '') + '\n' + rev_tail[1]
            
            c = '# Intra-chunk progress reporting\n'.join(parts)
        
        # 5. FWD FINISH
        fwd_fin_tgt = r'(\$pm_fwd->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*(?:pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,)?\s*\}\);\s*\} # End chunks\s*\$pm_fwd->wait_all_children\(\);)\s*unlink \$fwd_prog_file if -e \$fwd_prog_file;'
        fwd_fin_repl = r'my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1\n  use File::Path qw(remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;'
        c = re.sub(fwd_fin_tgt, fwd_fin_repl, c)

        # 6. REV FINISH
        rev_fin_tgt = r'(\$pm_rev->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*(?:pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,)?\s*\}\);\s*\} # End chunks\s*\$pm_rev->wait_all_children\(\);)\s*unlink \$rev_prog_file if -e \$rev_prog_file;'
        rev_fin_repl = r'my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1\n  use File::Path qw(remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;'
        c = re.sub(rev_fin_tgt, rev_fin_repl, c)

        with open(fname, 'w') as f:
            f.write(c)

fix_all()
