import os

def process_file(fname, is_stem=False):
    with open(fname, 'r') as f:
        c = f.read()

    options_var = "$options" if is_stem else "$options_r->"

    # 1. Init FWD
    old_init_fwd = f"""  my $fwd_prog_file = "{options_var}{{'output_file'}}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\\n"; }} 
      close $fh;
  }}"""
    old_init_fwd2 = f"""  my $fwd_prog_file = "{options_var}{{'output_file'}}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\\n"; }}
      close $fh;
  }}"""
    new_init_fwd = f"""  my $fwd_prog_dir = "{options_var}{{'output_file'}}_fwd_prog_$$";
  use File::Path qw(make_path remove_tree);
  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;
  make_path($fwd_prog_dir);"""
    c = c.replace(old_init_fwd, new_init_fwd).replace(old_init_fwd2, new_init_fwd)

    # 2. Init REV
    old_init_rev = f"""  my $rev_prog_file = "{options_var}{{'output_file'}}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\\n"; }} 
      close $fh;
  }}"""
    old_init_rev2 = f"""  my $rev_prog_file = "{options_var}{{'output_file'}}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\\n"; }}
      close $fh;
  }}"""
    new_init_rev = f"""  my $rev_prog_dir = "{options_var}{{'output_file'}}_rev_prog_$$";
  use File::Path qw(make_path remove_tree);
  remove_tree($rev_prog_dir) if -d $rev_prog_dir;
  make_path($rev_prog_dir);"""
    c = c.replace(old_init_rev, new_init_rev).replace(old_init_rev2, new_init_rev)

    # 3. Report FWD
    old_rep_fwd = """          # Intra-chunk progress reporting
          $chunk_done++;
          if ($chunk_done % 5 == 0 || $chunk_done == $fwd_chunk_size) {
              if (open(my $fh, '>>', $fwd_prog_file)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n";
                  close($fh);
                  
                  if (open(my $fh_read, '<', $fwd_prog_file)) {
                      my $total_done = 0;
                      my $total_hits = 0;
                      while(<$fh_read>) { 
                          chomp; 
                          next unless $_;
                          my ($d, $h) = split /,/, $_;
                          $total_done += $d;
                          $total_hits += $h;
                      }
                      close($fh_read);"""
    
    new_rep_fwd = """          # Intra-chunk progress reporting
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
              }"""
    c = c.replace(old_rep_fwd, new_rep_fwd)

    # 4. Report REV
    old_rep_rev = old_rep_fwd.replace('$fwd_prog_file', '$rev_prog_file').replace('$fwd_chunk_size', '$rev_chunk_size')
    new_rep_rev = new_rep_fwd.replace('$fwd_prog_dir', '$rev_prog_dir')
    c = c.replace(old_rep_rev, new_rep_rev)

    # 5. Finish FWD
    old_fin_fwd = """} # End Inner chunk loop
      
      $pm_fwd->finish(0, {
          infos => \\%chunk_infos,
          penalties => \\%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_fwd->wait_all_children();
  unlink $fwd_prog_file if -e $fwd_prog_file;"""
    
    new_fin_fwd = """} # End Inner chunk loop
      
      my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
      if (open(my $fh, '>', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }
      
      $pm_fwd->finish(0, {
          infos => \\%chunk_infos,
          penalties => \\%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_fwd->wait_all_children();
  use File::Path qw(remove_tree);
  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;"""
    c = c.replace(old_fin_fwd, new_fin_fwd)

    # 6. Finish REV
    old_fin_rev = """} # End Inner chunk loop
      
      $pm_rev->finish(0, {
          infos => \\%chunk_infos,
          penalties => \\%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_rev->wait_all_children();
  unlink $rev_prog_file if -e $rev_prog_file;"""
    
    new_fin_rev = """} # End Inner chunk loop
      
      my $prog_file_me = "$rev_prog_dir/chunk_$chunk_id.prog";
      if (open(my $fh, '>', $prog_file_me)) { flock($fh, 2); print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n"; close($fh); }
      
      $pm_rev->finish(0, {
          infos => \\%chunk_infos,
          penalties => \\%chunk_penalties,
          hits => $chunk_hits,
          done => $chunk_done,
          pruned => $chunk_pruned,
          evaluated => $chunk_evaluated,
      });
  } # End chunks
  $pm_rev->wait_all_children();
  use File::Path qw(remove_tree);
  remove_tree($rev_prog_dir) if -d $rev_prog_dir;"""
    c = c.replace(old_fin_rev, new_fin_rev)

    with open(fname, 'w') as f:
        f.write(c)

process_file('lava_loop_primer.pl', is_stem=False)
process_file('lava_stem_primer.pl', is_stem=True)
