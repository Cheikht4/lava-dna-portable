import re
import os

def fix_perl_file(fname, is_stem=False):
    with open(fname, 'r') as f:
        c = f.read()

    opt = "options" if is_stem else "options_r->"

    # 1. Init FWD
    c = re.sub(
        r'my \$fwd_prog_file\s*=\s*"' + re.escape(f"${opt}{{") + r'\'?output_file\'?}_fwd_prog_\$\$\.txt";\s*if\s*\(open\s*my\s*\$fh,\s*\'>\',\s*\$fwd_prog_file\)\s*\{\s*for\s*\(1\.\.30\)\s*\{\s*print\s*\$fh\s*"0,0,0,0\\n";\s*\}\s*close\s*\$fh;\s*\}',
        f'my $fwd_prog_dir = "${opt}{{\'output_file\'}}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);',
        c
    )
    
    # 2. Init REV
    c = re.sub(
        r'my \$rev_prog_file\s*=\s*"' + re.escape(f"${opt}{{") + r'\'?output_file\'?}_rev_prog_\$\$\.txt";\s*if\s*\(open\s*my\s*\$fh,\s*\'>\',\s*\$rev_prog_file\)\s*\{\s*for\s*\(1\.\.30\)\s*\{\s*print\s*\$fh\s*"0,0,0,0\\n";\s*\}\s*close\s*\$fh;\s*\}',
        f'my $rev_prog_dir = "${opt}{{\'output_file\'}}_rev_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);',
        c
    )

    # 3. Report FWD
    fwd_rep_tgt = r'\$chunk_done\+\+;\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$fwd_chunk_size\) \{\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(open\(my \$fh, \'>>\', \$fwd_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$fwd_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*(?:my \$total_pruned = 0;\s*my \$total_eval = 0;\s*)?while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d,\s*\$h(?:,\s*\$p,\s*\$e)?\) = split /,/,\s*\$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*(?:\$total_pruned \+= \$p;\s*\$total_eval \+= \$e;\s*)?\}\s*close\(\$fh_read\);'
    
    fwd_rep_repl = r"""if ($chunk_done % 5 == 0) {
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
                      my ($d, $h) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                  }
              }"""
    c = re.sub(fwd_rep_tgt, fwd_rep_repl, c)

    # 4. Report REV
    rev_rep_tgt = r'\$chunk_done\+\+;\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$rev_chunk_size\) \{\s*(?:# --- Intra-chunk progress reporting ---\s*)?if \(open\(my \$fh, \'>>\', \$rev_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$rev_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*(?:my \$total_pruned = 0;\s*my \$total_eval = 0;\s*)?while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d,\s*\$h(?:,\s*\$p,\s*\$e)?\) = split /,/,\s*\$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*(?:\$total_pruned \+= \$p;\s*\$total_eval \+= \$e;\s*)?\}\s*close\(\$fh_read\);'
    rev_rep_repl = fwd_rep_repl.replace('fwd_prog_dir', 'rev_prog_dir')
    c = re.sub(rev_rep_tgt, rev_rep_repl, c)

    # 5. Finish FWD
    fwd_fin_tgt = r'(\$pm_fwd->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*(?:pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,)?\s*\}\);\s*\} # End chunks\s*\$pm_fwd->wait_all_children\(\);)\s*unlink \$fwd_prog_file if -e \$fwd_prog_file;'
    fwd_fin_repl = r'my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";\n      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \1\n  use File::Path qw(remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;'
    c = re.sub(fwd_fin_tgt, fwd_fin_repl, c)

    # 6. Finish REV
    rev_fin_tgt = r'(\$pm_rev->finish\(0, \{\s*infos => \\%chunk_infos,\s*penalties => \\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*(?:pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,)?\s*\}\);\s*\} # End chunks\s*\$pm_rev->wait_all_children\(\);)\s*unlink \$rev_prog_file if -e \$rev_prog_file;'
    rev_fin_repl = fwd_fin_repl.replace('fwd', 'rev').replace('FWD', 'REV')
    c = re.sub(rev_fin_tgt, rev_fin_repl, c)

    with open(fname, 'w') as f:
        f.write(c)

fix_perl_file('lava_loop_primer.pl', is_stem=False)
fix_perl_file('lava_stem_primer.pl', is_stem=True)
