#!/usr/bin/env perl
use strict;
use warnings;

# This is a test script that would work standalone, but thanks to the Snakemake preamble,
# it can obtain all the invocational arguments, and be used as a Snakemake Perl script
# without needing any modification.

my @keep;
my %vars;
my %flags;
for (@ARGV) {
    if (/^--?(.+?)=(.+)$/) {
        $vars{$1} = $2;
    }
    elsif (/^--?(.+)$/) {
        $flags{$1} = 1;
    }
    else {
        push @keep, $_;
    }
}
@ARGV = @keep;

use Test::Simple tests => 1;

# --foo
ok( $flags{foo}, "Flag 'foo' is correctly passed" );

$vars{min} = 12 unless defined $vars{min};
$vars{max} = 50 unless defined $vars{max};

while (<>) {
    my $len = length $_;
    if ( $len >= $vars{min} && $len <= $vars{max} ) {
        print $_;
    }
    else {
        warn
"INFO: filtering line because its length ($len) is not between $vars{min} and $vars{max}\n";
    }
}
