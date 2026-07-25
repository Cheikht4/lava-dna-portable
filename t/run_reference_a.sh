#!/bin/bash
export PERL5LIB="./lib:$PERL5LIB"
echo "Running Test A (Stem)..."
perl lava_stem_primer.pl \
  --alignment_fasta t/fixtures/single_400.fasta \
  --threads 1 \
  --signature_max_length 400 \
  --window_size 20 \
  --max_per_window 3 \
  --output_file t/reference_a_stem \
  --loop_min_gap 20 \
  --dna_conc 400.0 --dntp_conc 1.4 --entropy_threshold 1.5 \
  --inner_primer_min_length 15 --inner_primer_target_length 18 --inner_primer_max_length 22 \
  --inner_primer_min_tm 59.0 --inner_primer_target_tm 60.0 --inner_primer_max_tm 65.0 \
  --middle_primer_min_length 15 --middle_primer_target_length 18 --middle_primer_max_length 22 \
  --middle_primer_min_tm 59.0 --middle_primer_target_tm 60.0 --middle_primer_max_tm 65.0 \
  --outer_primer_min_length 15 --outer_primer_target_length 18 --outer_primer_max_length 22 \
  --outer_primer_min_tm 57.0 --outer_primer_target_tm 58.0 --outer_primer_max_tm 59.0 \
  --loop_primer_min_length 15 --loop_primer_target_length 18 --loop_primer_max_length 22 \
  --loop_primer_min_tm 59.0 --loop_primer_target_tm 60.0 --loop_primer_max_tm 61.0 \
  --max_dist_middle_inner 30 --max_dist_outer_middle 30 --min_primer_coverage 70.0 \
  --min_signatures_for_success 1 --resolve_overlap_by coverage \
  --max_tm_diff 5.0 --primer_iupac_min_percent 80.0 --primer_min_match_percent 70.0 \
  --min_base_frequency 0.2 --three_prime_zone_size 3 --salt_monovalent 50.0 \
  --salt_divalent 8.0 --penalty_plateau 0.25 --penalty_slope 0.15 \
  --max_overlap_percent 0.0 --max_primer_gen 10000.0 > t/reference_a_stem.log 2>&1

echo "Test A (Stem) Done."
