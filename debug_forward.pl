use strict;
use warnings;
my $innerLocation = 350;
my $signatureMaxLength = 320;
my $innerLength = 20;
my $innerPairTargetLength = 50;

my $searchStartAt = $innerLocation - $signatureMaxLength + $innerLength + 20;
my $minLoopLocation = $innerLocation - $innerPairTargetLength - 250;
$searchStartAt = $minLoopLocation if ($minLoopLocation > $searchStartAt);
$searchStartAt = 0 if $searchStartAt < 0;

my $loopEndAt = $innerLocation - 1 - 1;

my $minMiddleLocation = $innerLocation - $innerPairTargetLength - 250;
my $middleStartAt = $searchStartAt;
$middleStartAt = $minMiddleLocation if ($minMiddleLocation > $middleStartAt);

my $loopLocation = 348;
my $loopLength = 20;
my $middleEndAt = $loopLocation - $loopLength - 1;

my $middleLocation = 186;

print "Inner: $innerLocation\n";
print "Loop bounds: start=$searchStartAt, end=$loopEndAt\n";
print "Loop is at: $loopLocation\n";
print "Middle bounds: start=$middleStartAt, end=$middleEndAt\n";
print "Middle is at: $middleLocation\n";

if ($middleLocation < $middleStartAt) { print "Middle rejected by start bound!\n"; }
if ($middleLocation > $middleEndAt) { print "Middle rejected by end bound!\n"; }
