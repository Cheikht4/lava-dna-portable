#!/usr/bin/env perl
use strict;
use warnings;

my $b1 = "";
vec($b1, 1, 1) = 1;
vec($b1, 5, 1) = 1;
vec($b1, 10, 1) = 1;

my $b2 = "";
vec($b2, 5, 1) = 1;
vec($b2, 10, 1) = 1;
vec($b2, 20, 1) = 1;

my $b3 = $b1 & $b2;
my $count = unpack("%32b*", $b3);

print "Count: $count\n";
