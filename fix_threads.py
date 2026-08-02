import sys

def patch_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Replace Parallel::ForkManager with LLNL::LAVA::ForkManager and get actual_threads
    content = content.replace(
        'my $val_pm = Parallel::ForkManager->new($options_r->{"threads"});\n  my %validation_results;',
        'my $val_pm = LLNL::LAVA::ForkManager->new($options_r->{"threads"});\n  my $actual_threads = $val_pm->{max_processes};\n  my %validation_results;'
    )
    
    # Replace $options_r->{"threads"} with $actual_threads in chunk size
    content = content.replace(
        'my $val_chunk_size = POSIX::ceil($total_sigs_to_validate / ($options_r->{"threads"} * 4));',
        'my $val_chunk_size = POSIX::ceil($total_sigs_to_validate / ($actual_threads * 4));'
    )
    
    # Replace in loop 1
    content = content.replace(
        'for(my $i = 0; $i < $options_r->{"threads"}; $i++) {',
        'for(my $i = 0; $i < $actual_threads; $i++) {'
    )
    
    # Replace in loop 2
    content = content.replace(
        'my $worker_idx = $i % $options_r->{"threads"};',
        'my $worker_idx = $i % $actual_threads;'
    )
    
    # Replace in loop 3
    content = content.replace(
        'for(my $w = 0; $w < $options_r->{"threads"}; $w++) {',
        'for(my $w = 0; $w < $actual_threads; $w++) {'
    )
    
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"Fixed threads in {filepath}")

patch_file("lava_loop_primer.pl")
patch_file("lava_stem_primer.pl")
