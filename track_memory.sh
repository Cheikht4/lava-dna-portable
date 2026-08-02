#!/bin/bash
perl lava_loop_primer.pl --alignment_fasta t/fixtures/rota_canary_40.fasta --threads 4 > /dev/null &
PID=$!
MAX_RSS=0

while kill -0 $PID 2>/dev/null; do
    # Get RSS of all perl processes running lava_loop_primer
    RSS_TOTAL=0
    while read -r rss cmd; do
        if [[ $cmd == *"lava_loop_primer.pl"* ]]; then
            RSS_TOTAL=$((RSS_TOTAL + rss))
        fi
    done < <(ps -eo rss,command | grep lava_loop_primer | grep -v grep)
    
    if [ $RSS_TOTAL -gt $MAX_RSS ]; then
        MAX_RSS=$RSS_TOTAL
    fi
    sleep 0.2
done

echo "Max total RSS during run (KB): $MAX_RSS"
