import re
import os

def process_file(fname, is_stem=False):
    with open(fname, 'r') as f:
        c = f.read()

    opt = "options" if is_stem else "options_r->"

    # 1. Init FWD
    init_fwd_tgt = r'my \$fwd_prog_file = "\$' + opt.replace('->', r'->') + r'\{\'?output_file\'?\}_fwd_prog_\$\$\.txt";\n\s*if \(open my \$fh, \'>\', \$fwd_prog_file\) \{\n\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\n"; \}\s*close \$fh;\n\s*\}'
    init_fwd_repl = f'my $fwd_prog_dir = "${opt}{{\'output_file\'}}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);'
    c = re.sub(init_fwd_tgt, init_fwd_repl, c)

    # 2. Init REV
    init_rev_tgt = init_fwd_tgt.replace('fwd_', 'rev_')
    init_rev_repl = init_fwd_repl.replace('fwd_', 'rev_')
    c = re.sub(init_rev_tgt, init_rev_repl, c)

    # 3. Report FWD
    fwd_rep_tgt = r'# Intra-chunk progress reporting\n\s*\$chunk_done\+\+;\n\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$fwd_chunk_size\) \{.*?# End Inner chunk loop'
    
    fwd_rep_repl = r"""# Intra-chunk progress reporting
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
                  printf("[LAVA-PROGRESS] Signatures Forward|%d|%d|Sig: %d|%.1f it/s|%d\\n", $total_done, $innerForwardCount, $total_hits, $rate, $eta);
                  my $old_h = select(STDOUT); $| = 1; select($old_h);
              }
          }
} # End Inner chunk loop"""
    c = re.sub(fwd_rep_tgt, fwd_rep_repl, c, flags=re.DOTALL)

    # 4. Report REV
    rev_rep_tgt = r'# Intra-chunk progress reporting\n\s*\$chunk_done\+\+;\n\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$rev_chunk_size\) \{.*?# End Inner chunk loop'
    rev_rep_repl = fwd_rep_repl.replace('fwd_prog_dir', 'rev_prog_dir').replace('fwd_chunk_size', 'rev_chunk_size').replace('innerForwardCount', 'innerReverseCount').replace('Forward', 'Reverse').replace('_fwd_', '_rev_')
    c = re.sub(rev_rep_tgt, rev_rep_repl, c, flags=re.DOTALL)

    # 5. Finish FWD
    fwd_fin_tgt = r'(\$pm_fwd->finish\(0, \{.*?\}\);\n\s*\} # End chunks\n\s*\$pm_fwd->wait_all_children\(\);)\n\s*unlink \$fwd_prog_file if -e \$fwd_prog_file;'
    fwd_fin_repl = r'my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1\n  use File::Path qw(remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;'
    c = re.sub(fwd_fin_tgt, fwd_fin_repl, c, flags=re.DOTALL)

    # 6. Finish REV
    rev_fin_tgt = r'(\$pm_rev->finish\(0, \{.*?\}\);\n\s*\} # End chunks\n\s*\$pm_rev->wait_all_children\(\);)\n\s*unlink \$rev_prog_file if -e \$rev_prog_file;'
    rev_fin_repl = r'my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1\n  use File::Path qw(remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;'
    c = re.sub(rev_fin_tgt, rev_fin_repl, c, flags=re.DOTALL)

    with open(fname, 'w') as f:
        f.write(c)

process_file('lava_loop_primer.pl', is_stem=False)
process_file('lava_stem_primer.pl', is_stem=True)
