import os

def fix_all():
    for fname in ['lava_loop_primer.pl', 'lava_stem_primer.pl', 'apply_optimizations.py']:
        with open(fname, 'r') as f:
            lines = f.readlines()

        opt = "$options" if 'stem' in fname else "$options_r->"

        new_lines = []
        in_prog_block = False
        in_fwd_fin_block = False
        in_rev_fin_block = False
        
        block_idx = 0
        
        for i, line in enumerate(lines):
            # 1. FWD INIT
            if 'my $fwd_prog_file = "$options' in line and '_fwd_prog_$$.txt' in line:
                new_lines.append(f'  my $fwd_prog_dir = "{opt}{{\'output_file\'}}_fwd_prog_$$";\n')
                new_lines.append(f'  $fwd_prog_dir = "$options_r->{{\'output_file\'}}_fwd_prog_$$" if ref($options_r);\n')
                new_lines.append('  use File::Path qw(make_path remove_tree);\n')
                new_lines.append('  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n')
                new_lines.append('  make_path($fwd_prog_dir);\n')
                continue
            
            # 2. REV INIT
            if 'my $rev_prog_file = "$options' in line and '_rev_prog_$$.txt' in line:
                new_lines.append(f'  my $rev_prog_dir = "{opt}{{\'output_file\'}}_rev_prog_$$";\n')
                new_lines.append(f'  $rev_prog_dir = "$options_r->{{\'output_file\'}}_rev_prog_$$" if ref($options_r);\n')
                new_lines.append('  use File::Path qw(make_path remove_tree);\n')
                new_lines.append('  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n')
                new_lines.append('  make_path($rev_prog_dir);\n')
                continue
                
            if i > 0 and 'if (open my $fh, \'>\', $fwd_prog_file)' in lines[i-1]: continue
            if i > 1 and 'if (open my $fh, \'>\', $fwd_prog_file)' in lines[i-2]: continue
            if i > 2 and 'if (open my $fh, \'>\', $fwd_prog_file)' in lines[i-3]: continue
            if i > 3 and 'if (open my $fh, \'>\', $fwd_prog_file)' in lines[i-4]: continue
            # Handle the literal newline case (6 lines block instead of 5)
            if i > 4 and 'if (open my $fh, \'>\', $fwd_prog_file)' in lines[i-5]: continue

            if i > 0 and 'if (open my $fh, \'>\', $rev_prog_file)' in lines[i-1]: continue
            if i > 1 and 'if (open my $fh, \'>\', $rev_prog_file)' in lines[i-2]: continue
            if i > 2 and 'if (open my $fh, \'>\', $rev_prog_file)' in lines[i-3]: continue
            if i > 3 and 'if (open my $fh, \'>\', $rev_prog_file)' in lines[i-4]: continue
            if i > 4 and 'if (open my $fh, \'>\', $rev_prog_file)' in lines[i-5]: continue

            if 'if (open my $fh, \'>\', $fwd_prog_file)' in line: continue
            if 'if (open my $fh, \'>\', $rev_prog_file)' in line: continue
            
            # Also need to drop the dangling '"; }' if it exists. But we already skip it.

            # 3. & 4. PROGRESS BLOCKS
            if '# Intra-chunk progress reporting' in line:
                in_prog_block = True
                block_idx += 1
                
                is_fwd = (block_idx == 1)
                pdir = '$fwd_prog_dir' if is_fwd else '$rev_prog_dir'
                ic = '$innerForwardCount' if is_fwd else '$innerReverseCount'
                name = 'Forward' if is_fwd else 'Reverse'
                
                new_lines.append("          # Intra-chunk progress reporting\n")
                new_lines.append("          if ($chunk_done % 5 == 0) {\n")
                new_lines.append(f"              my $prog_file_me = \"{pdir}/chunk_$chunk_id.prog\";\n")
                new_lines.append("              if (open(my $fh, '>', $prog_file_me)) {\n")
                new_lines.append("                  flock($fh, 2);\n")
                new_lines.append("                  print $fh \"$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n\";\n")
                new_lines.append("                  close($fh);\n")
                new_lines.append("              }\n")
                new_lines.append("                  \n")
                new_lines.append("              my $total_done = 0;\n")
                new_lines.append("              my $total_hits = 0;\n")
                new_lines.append(f"              foreach my $f (glob(\"{pdir}/chunk_*.prog\")) {{\n")
                new_lines.append("                  if (open(my $r, '<', $f)) {\n")
                new_lines.append("                      my $line = <$r>; close($r);\n")
                new_lines.append("                      next unless defined $line;\n")
                new_lines.append("                      chomp $line;\n")
                new_lines.append("                      my ($d, $h) = split /,/, $line;\n")
                new_lines.append("                      $total_done += $d // 0;\n")
                new_lines.append("                      $total_hits += $h // 0;\n")
                new_lines.append("                  }\n")
                new_lines.append("              }\n")
                new_lines.append("                      \n")
                new_lines.append("              if ($_LAVA_IS_TTY || 1) {\n")
                t0 = '$_sig_fwd_t0' if is_fwd else '$_sig_rev_t0'
                new_lines.append(f"                  my $elapsed = time() - {t0} + 0.001;\n")
                new_lines.append(f"                  my $eta = ($total_done < {ic}) ? int(({ic} - $total_done) / ($total_done / $elapsed)) : 0;\n")
                new_lines.append("                  my $rate = $total_done / $elapsed;\n")
                if 'apply_optimizations.py' in fname:
                    new_lines.append(f"                  printf(\"[LAVA-PROGRESS] Signatures {name}|%d|%d|Sig: %d|%.1f it/s|%d\\r\", $total_done, {ic}, $total_hits, $rate, $eta);\n")
                else:
                    new_lines.append(f"                  printf(\"[LAVA-PROGRESS] Signatures {name}|%d|%d|Sig: %d|%.1f it/s|%d\\r\", $total_done, {ic}, $total_hits, $rate, $eta);\n")
                new_lines.append("                  my $old_h = select(STDOUT); $| = 1; select($old_h);\n")
                new_lines.append("              }\n")
                new_lines.append("          }\n")
                continue
                
            if in_prog_block:
                if '} # End Inner chunk loop' in line or '} # End reverse inner chunk loop' in line or '} # End Inner Chunk Loop' in line or '} # End inner chunk loop' in line or '} # End forward inner chunk loop' in line:
                    in_prog_block = False
                    new_lines.append(line)
                continue
                
            # 5. & 6. FINISH
            if '$pm_fwd->finish(0, {' in line:
                in_fwd_fin_block = True
                new_lines.append('      my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";\n')
                new_lines.append('      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \n')
                new_lines.append(line)
                continue
                
            if in_fwd_fin_block:
                if 'unlink $fwd_prog_file if -e $fwd_prog_file;' in line:
                    in_fwd_fin_block = False
                    new_lines.append('  use File::Path qw(remove_tree);\n')
                    new_lines.append('  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n')
                    continue
                new_lines.append(line)
                continue

            if '$pm_rev->finish(0, {' in line:
                in_rev_fin_block = True
                new_lines.append('      my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";\n')
                new_lines.append('      if (open(my $fh, \'>\', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }\n      \n')
                new_lines.append(line)
                continue
                
            if in_rev_fin_block:
                if 'unlink $rev_prog_file if -e $rev_prog_file;' in line:
                    in_rev_fin_block = False
                    new_lines.append('  use File::Path qw(remove_tree);\n')
                    new_lines.append('  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n')
                    continue
                new_lines.append(line)
                continue
                
            new_lines.append(line)
            
        with open(fname, 'w') as f:
            f.writelines(new_lines)

fix_all()
