import os

def process_perl(fname, is_stem=False):
    with open(fname, 'r') as f:
        c = f.read()
    
    opt = "$options" if is_stem else "$options_r->"

    # FWD INIT
    old_fwd_init = f"""  my $fwd_prog_file = "{opt}{{'output_file'}}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\n"; }}
      close $fh;
  }}"""
    old_fwd_init2 = f"""  my $fwd_prog_file = "{opt}{{'output_file'}}_fwd_prog_$$.txt";
  if (open my $fh, '>', $fwd_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\n"; }} 
      close $fh;
  }}"""
    new_fwd_init = f"""  my $fwd_prog_dir = "{opt}{{'output_file'}}_fwd_prog_$$";
  use File::Path qw(make_path remove_tree);
  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;
  make_path($fwd_prog_dir);"""
    c = c.replace(old_fwd_init, new_fwd_init).replace(old_fwd_init2, new_fwd_init)

    # REV INIT
    old_rev_init = f"""  my $rev_prog_file = "{opt}{{'output_file'}}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\n"; }}
      close $fh;
  }}"""
    old_rev_init2 = f"""  my $rev_prog_file = "{opt}{{'output_file'}}_rev_prog_$$.txt";
  if (open my $fh, '>', $rev_prog_file) {{
      for (1..30) {{ print $fh "0,0,0,0\n"; }} 
      close $fh;
  }}"""
    new_rev_init = f"""  my $rev_prog_dir = "{opt}{{'output_file'}}_rev_prog_$$";
  use File::Path qw(make_path remove_tree);
  remove_tree($rev_prog_dir) if -d $rev_prog_dir;
  make_path($rev_prog_dir);"""
    c = c.replace(old_rev_init, new_rev_init).replace(old_rev_init2, new_rev_init)

    # FWD REPORT
    old_fwd_rep = """          # Intra-chunk progress reporting
          $chunk_done++;
          if ($chunk_done % 5 == 0 || $chunk_done == $fwd_chunk_size) {
              if (open(my $fh, '>>', $fwd_prog_file)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\n";
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
    
    new_fwd_rep = """          # Intra-chunk progress reporting
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
    c = c.replace(old_fwd_rep, new_fwd_rep)

    # REV REPORT
    old_rev_rep = old_fwd_rep.replace('$fwd_prog_file', '$rev_prog_file').replace('$fwd_chunk_size', '$rev_chunk_size')
    new_rev_rep = new_fwd_rep.replace('$fwd_prog_dir', '$rev_prog_dir')
    c = c.replace(old_rev_rep, new_rev_rep)

    # FWD FINISH
    old_fwd_fin = """} # End Inner chunk loop
      
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
    
    new_fwd_fin = """} # End Inner chunk loop
      
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
    c = c.replace(old_fwd_fin, new_fwd_fin)

    # REV FINISH
    old_rev_fin = old_fwd_fin.replace('fwd_prog_file', 'rev_prog_file').replace('$pm_fwd', '$pm_rev')
    new_rev_fin = new_fwd_fin.replace('fwd_prog_dir', 'rev_prog_dir').replace('$pm_fwd', '$pm_rev')
    c = c.replace(old_rev_fin, new_rev_fin)

    with open(fname, 'w') as f:
        f.write(c)

process_perl('lava_loop_primer.pl', False)
process_perl('lava_stem_primer.pl', True)
