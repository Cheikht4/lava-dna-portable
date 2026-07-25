import re

with open('apply_optimizations.py', 'r') as f:
    c = f.read()

# 1. Fix prog file init fwd Loop
c = re.sub(
    r'my \$fwd_prog_file = "\$options_r->\{\'output_file\'\}_fwd_prog_\$\$\.txt";\s*if \(open my \$fh, \'>\', \$fwd_prog_file\) \{\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\\\n"; \}\s*close \$fh;\s*\}',
    r'my $fwd_prog_dir = "$options_r->{\'output_file\'}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);',
    c
)

# 2. Fix prog file init fwd Stem
c = re.sub(
    r'my \$fwd_prog_file = "\$options\{\'output_file\'\}_fwd_prog_\$\$\.txt";\s*if \(open my \$fh, \'>\', \$fwd_prog_file\) \{\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\\\n"; \}\s*close \$fh;\s*\}',
    r'my $fwd_prog_dir = "$options{\'output_file\'}_fwd_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($fwd_prog_dir) if -d $fwd_prog_dir;\n  make_path($fwd_prog_dir);',
    c
)

# 3. Fix prog file init rev Loop
c = re.sub(
    r'my \$rev_prog_file = "\$options_r->\{\'output_file\'\}_rev_prog_\$\$\.txt";\s*if \(open my \$fh, \'>\', \$rev_prog_file\) \{\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\\\n"; \}\s*close \$fh;\s*\}',
    r'my $rev_prog_dir = "$options_r->{\'output_file\'}_rev_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);',
    c
)

# 4. Fix prog file init rev Stem
c = re.sub(
    r'my \$rev_prog_file = "\$options\{\'output_file\'\}_rev_prog_\$\$\.txt";\s*if \(open my \$fh, \'>\', \$rev_prog_file\) \{\s*for \(1\.\.30\) \{ print \$fh "0,0,0,0\\\\n"; \}\s*close \$fh;\s*\}',
    r'my $rev_prog_dir = "$options{\'output_file\'}_rev_prog_$$";\n  use File::Path qw(make_path remove_tree);\n  remove_tree($rev_prog_dir) if -d $rev_prog_dir;\n  make_path($rev_prog_dir);',
    c
)


# 5. Fix chunk finish block FWD
old_fwd = r'\$chunk_done\+\+;\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$fwd_chunk_size\) \{\s*if \(open\(my \$fh, \'>>\', \$fwd_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$fwd_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d, \$h\) = split /,\/, \$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*\}\s*close\(\$fh_read\);'

new_fwd = r"""if ($chunk_done % 5 == 0) {
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
c = re.sub(old_fwd, new_fwd, c)


# 6. Fix chunk finish block REV
old_rev = r'\$chunk_done\+\+;\s*if \(\$chunk_done \% 5 == 0 \|\| \$chunk_done == \$rev_chunk_size\) \{\s*if \(open\(my \$fh, \'>>\', \$rev_prog_file\)\) \{\s*flock\(\$fh, 2\);\s*print \$fh "\$chunk_done,\$chunk_hits,\$chunk_pruned,\$chunk_evaluated\\n";\s*close\(\$fh\);\s*if \(open\(my \$fh_read, \'<\', \$rev_prog_file\)\) \{\s*my \$total_done = 0;\s*my \$total_hits = 0;\s*while\(<\$fh_read>\) \{\s*chomp;\s*next unless \$_;\s*my \(\$d, \$h\) = split /,\/, \$_;\s*\$total_done \+= \$d;\s*\$total_hits \+= \$h;\s*\}\s*close\(\$fh_read\);'
new_rev = new_fwd.replace('fwd_prog_dir', 'rev_prog_dir')
c = re.sub(old_rev, new_rev, c)

# 7. Add final progress file write & fix unlink for FWD
old_unlink_fwd = r'\\2\s*\$pm_fwd->finish\(0, \{\s*infos => \\\%chunk_infos,\s*penalties => \\\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,\s*\}\);\s*\} # End chunks\s*\$pm_fwd->wait_all_children\(\);\s*unlink \$fwd_prog_file if -e \$fwd_prog_file;'
new_unlink_fwd = r"""\\2      
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
c = re.sub(old_unlink_fwd, new_unlink_fwd, c)

# 8. Add final progress file write & fix unlink for REV
old_unlink_rev = r'\\2\s*\$pm_rev->finish\(0, \{\s*infos => \\\%chunk_infos,\s*penalties => \\\%chunk_penalties,\s*hits => \$chunk_hits,\s*done => \$chunk_done,\s*pruned => \$chunk_pruned,\s*evaluated => \$chunk_evaluated,\s*\}\);\s*\} # End chunks\s*\$pm_rev->wait_all_children\(\);\s*unlink \$rev_prog_file if -e \$rev_prog_file;'
new_unlink_rev = new_unlink_fwd.replace('fwd', 'rev').replace('FWD', 'REV')
c = re.sub(old_unlink_rev, new_unlink_rev, c)

# Save
with open('apply_optimizations.py', 'w') as f:
    f.write(c)
