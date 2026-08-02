#!/bin/bash
echo "--- TEST DENGUE 2 (Max 130) ---"
perl lava_loop_primer.pl --alignment_fasta fixtures/dengue2_aligned.fasta --output_file scratch/d2_max130 --signature_max_length 130 --signature_min_length 0 > scratch/max130.log 2>&1
grep "Assemblage :" scratch/max130.log || echo "Not found"

echo "--- TEST DENGUE 2 (Min 180) ---"
perl lava_loop_primer.pl --alignment_fasta fixtures/dengue2_aligned.fasta --output_file scratch/d2_min180 --signature_max_length 300 --signature_min_length 180 > scratch/min180.log 2>&1
grep "Assemblage :" scratch/min180.log || echo "Not found"
