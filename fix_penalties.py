import re

def add_penalty_at(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    if 'sub penaltyAt {' not in content:
        insert_pos = content.find('use Getopt::Long;')
        if insert_pos != -1:
            insert_pos = content.find('\n', insert_pos) + 1
        
        penalty_sub = '''
sub penaltyAt {
    my ($table_r, $distance, $label) = @_;
    if ($distance < 0) {
        warn "[PENALTY GUARD] distance negative ($distance) sur $label -> penalite max appliquee\\n";
        return 100;
    }
    return $table_r->[$distance] // 100;
}
'''
        content = content[:insert_pos] + penalty_sub + content[insert_pos:]
        
    pattern = r'\$([a-zA-Z0-9_]+Penalties_r)->\[\$([a-zA-Z0-9_]+Distance|innerSpacing)\]'
    
    def repl(m):
        arr = m.group(1)
        dist = m.group(2)
        label = arr.replace('Penalties_r', '')
        return f"penaltyAt(${arr}, ${dist}, '{label}')"
        
    content = re.sub(pattern, repl, content)
    
    with open(filepath, 'w') as f:
        f.write(content)

add_penalty_at('lava_loop_primer.pl')
add_penalty_at('lava_stem_primer.pl')
