use strict;
use warnings;
use Test::More tests => 3;
use File::Path qw(remove_tree);
use lib 'lib';

my $script = 'lava_loop_primer.pl';
my $fasta = 't/synthetic_entropy.fasta';
my $output = 't/canary_entropy_out.txt';

# Executer avec un seuil de 1.5 pour exclure la zone centrale tres variable
# Run with a threshold of 1.5 to exclude the highly variable central region
my $cmd = "perl -Ilib $script --alignment_fasta $fasta --output_file $output --entropy_threshold 1.5";
my $res = `$cmd 2>&1`;

# On verifie qu'aucune amorce n'a ete trouvee (car la zone est exclue et les flancs sont trop petits pour tout caser, ou bien l'exclusion bloque Primer3)
# Verify that no primer was found (as the central region is excluded and flanks are too small)
ok($res =~ /AUCUN OLIGO TROUVÉ/, "LAVA ne trouve aucune amorce quand la zone centrale est exclue");

# Verifier que le diagnostic mentionne la zone exclue
# Verify that the diagnostic mentions the excluded region
ok($res =~ /Régions exclues par entropie/, "Le diagnostic affiche les zones exclues par entropie");

# Executer avec un seuil de 5.0 (tres haut, rien ne sera exclu)
# Run with a threshold of 5.0 (very high, nothing will be excluded)
$cmd = "perl -Ilib $script --alignment_fasta $fasta --output_file $output --entropy_threshold 5.0";
$res = `$cmd 2>&1`;

ok($res !~ /AUCUN OLIGO TROUVÉ/, "LAVA trouve des amorces quand le seuil d'entropie est assez haut pour ne rien exclure");

# Nettoyage des fichiers temporaires generes par le test / Clean up temporary test output files
for my $p (glob("t/canary_entropy_out*")) {
    if (-d $p) {
        remove_tree($p);
    } else {
        unlink $p;
    }
}
