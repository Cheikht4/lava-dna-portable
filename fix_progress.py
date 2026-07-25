import re

with open('apply_optimizations.py', 'r') as f:
    c = f.read()

# Fix prog file init fwd Loop
c = c.replace(
    'my $fwd_prog_file = "$options_r->{\'output_file\'}_fwd_prog_$$.txt";\n  if (open my $fh, \'>\', $fwd_prog_file) {\n      for (1..30) { print $fh "0,0,0,0\\\\n"; }\n      close $fh;\n  }',
    'my $fwd_prog_dir = "$options_r->{\'output_file\'}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);'
)

# Fix prog file init fwd Stem
c = c.replace(
    'my $fwd_prog_file = "$options{\'output_file\'}_fwd_prog_$$.txt";\n  if (open my $fh, \'>\', $fwd_prog_file) {\n      for (1..30) { print $fh "0,0,0,0\\\\n"; } \n      close $fh;\n  }',
    'my $fwd_prog_dir = "$options{\'output_file\'}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);'
)

# Fix prog file init rev Loop
c = c.replace(
    'my $rev_prog_file = "$options_r->{\'output_file\'}_rev_prog_$$.txt";\n  if (open my $fh, \'>\', $rev_prog_file) {\n      for (1..30) { print $fh "0,0,0,0\\\\n"; }\n      close $fh;\n  }',
    'my $rev_prog_dir = "$options_r->{\'output_file\'}_rev_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);'
)

# Fix prog file init rev Stem
c = c.replace(
    'my $rev_prog_file = "$options{\'output_file\'}_rev_prog_$$.txt";\n  if (open my $fh, \'>\', $rev_prog_file) {\n      for (1..30) { print $fh "0,0,0,0\\\\n"; } \n      close $fh;\n  }',
    'my $rev_prog_dir = "$options{\'output_file\'}_rev_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);'
)

# Fix chunk finish fwd (remove $chunk_done++ and append loop)
old_fwd_finish = """          $chunk_done++;
          if ($chunk_done % 5 == 0 || $chunk_done == $fwd_chunk_size) {
              # --- Intra-chunk progress reporting ---
              if (open(my $fh, '>>', $fwd_prog_file)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n";
                  close($fh);
              }
              
              if (open(my $fh_read, '<', $fwd_prog_file)) {
                  my $total_done = 0;
                  my $total_hits = 0;
                  my $total_pruned = 0;
                  my $total_eval = 0;
                  while(<$fh_read>) {
                      chomp;
                      my ($d, $h, $p, $e) = split /,/;
                      $total_done += $d // 0;
                      $total_hits += $h // 0;
                      $total_pruned += $p // 0;
                      $total_eval += $e // 0;
                  }
                  close($fh_read);"""

new_fwd_finish = """          if ($chunk_done % 5 == 0) {
              # --- Intra-chunk progress reporting ---
              my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
              if (open(my $fh, '>', $prog_file_me)) {
                  flock($fh, 2);
                  print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n";
                  close($fh);
              }
              
              my $total_done = 0;
              my $total_hits = 0;
              my $total_pruned = 0;
              my $total_eval = 0;
              foreach my $f (glob("$fwd_prog_dir/chunk_*.prog")) {
                  if (open(my $r, '<', $f)) {
                      my $line = <$r>; close($r);
                      next unless defined $line;
                      chomp $line;
                      my ($d,$h,$p,$e) = split /,/, $line;
                      $total_done += $d // 0;
                      $total_hits  += $h // 0;
                      $total_pruned += $p // 0;
                      $total_eval += $e // 0;
                  }
              }"""

c = c.replace(old_fwd_finish, new_fwd_finish)

# Fix chunk finish rev
old_rev_finish = old_fwd_finish.replace('fwd_', 'rev_')
new_rev_finish = new_fwd_finish.replace('fwd_', 'rev_')
c = c.replace(old_rev_finish, new_rev_finish)

# Fix unlink and final finish fwd
old_fwd_unlink = """  } # End Inner chunk loop
  
  $pm_fwd->finish(0, { 
      infos => \\%chunk_infos,
      penalties => \\%chunk_penalties,
      hits => $chunk_hits,
      pruned => $chunk_pruned,
      evaluated => $chunk_evaluated,
      done => $chunk_done,
  });
} # End chunks

$pm_fwd->wait_all_children();

unlink $fwd_prog_file if -e $fwd_prog_file;"""

new_fwd_unlink = """  } # End Inner chunk loop
  
  my $prog_file_me = "$fwd_prog_dir/chunk_$chunk_id.prog";
  if (open(my $fh, '>', $prog_file_me)) {
      flock($fh, 2);
      print $fh "$chunk_done,$chunk_hits,$chunk_pruned,$chunk_evaluated\\n";
      close($fh);
  }
  
  $pm_fwd->finish(0, { 
      infos => \\%chunk_infos,
      penalties => \\%chunk_penalties,
      hits => $chunk_hits,
      pruned => $chunk_pruned,
      evaluated => $chunk_evaluated,
      done => $chunk_done,
  });
} # End chunks

$pm_fwd->wait_all_children();

use File::Path qw(remove_tree);
remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;"""

c = c.replace(old_fwd_unlink, new_fwd_unlink)

# Fix unlink and final finish rev
old_rev_unlink = old_fwd_unlink.replace('fwd_', 'rev_')
new_rev_unlink = new_fwd_unlink.replace('fwd_', 'rev_')
c = c.replace(old_rev_unlink, new_rev_unlink)

with open('apply_optimizations.py', 'w') as f:
    f.write(c)
