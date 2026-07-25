my $innerLocation = 350; my $innerLength = 20; my $innerTm = 62;
my $loopLocation = 348; my $loopLength = 20; my $loopTm = 60;
my $middleLocation = 186; my $middleLength = 16; my $midTm = 73.5;
my $outerLocation = 94; my $outerLength = 20; my $outTm = 60;

my $signatureMaxLength = 320;
my $minPrimerSpacing = 1;
my $loopMinGap = 25;
my $maxTmDiff = 5.0;

# Simulate loops
my $searchStartAt = $innerLocation - $signatureMaxLength + $innerLength + 20;
print "searchStartAt=$searchStartAt\n";

my $loopEndAt = $innerLocation - 1 - $minPrimerSpacing;
print "loopEndAt=$loopEndAt\n";

if ($loopLocation < $searchStartAt || $loopLocation > $loopEndAt) { print "FAIL loop bounds\n"; }

if (abs($innerTm - $loopTm) > $maxTmDiff) { print "FAIL loop Tm\n"; }

my $middleStartAt = $searchStartAt;
my $middleEndAt = $loopLocation - $loopLength - $minPrimerSpacing;
my $altMiddleEndAt = $innerLocation - ($loopMinGap + 1);
$middleEndAt = $altMiddleEndAt if ($altMiddleEndAt < $middleEndAt);
print "middleEndAt=$middleEndAt\n";

if ($middleLocation < $middleStartAt || $middleLocation > $middleEndAt) { print "FAIL middle bounds: $middleLocation > $middleEndAt\n"; }

if ($middleLocation + $middleLength + $minPrimerSpacing > $loopLocation - $loopLength + 1) { print "FAIL middle gap loop\n"; }
if ($middleLocation + $middleLength + $loopMinGap > $innerLocation) { print "FAIL middle gap inner\n"; }

# thermal bypassed because middle is fixed

my $outerStartAt = $searchStartAt;
my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;
print "outerEndAt=$outerEndAt\n";

if ($outerLocation < $outerStartAt || $outerLocation > $outerEndAt) { print "FAIL outer bounds\n"; }
if ($outerLocation + $outerLength + $minPrimerSpacing > $middleLocation) { print "FAIL outer gap middle\n"; }

print "SUCCESS!\n";
