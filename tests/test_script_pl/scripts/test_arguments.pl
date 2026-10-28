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

use Test::Simple tests => 4;

# --foo --bar --min=100 --max=200
ok( $flags{foo},                 "Flag 'foo' is correctly passed" );
ok( $vars{bar} eq 'foo bar baz', "Variable 'bar' is correctly passed" );
ok( $vars{min} == 12,            "Variable 'min' is correctly passed" );
ok( $vars{max} == 50,            "Variable 'max' is correctly passed" );

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
