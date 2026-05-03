#!/usr/bin/env perl
use strict;
use warnings;

# open (my $out, '>', 'perl.out') or die "Could not open file: $!";
# if (@{ $snakemake->{output}->{positional} }) {
#     open ($out, '>', $snakemake->{output}->{positional}->[0]) or die "Could not open file: $!";
# }

# Write to log file if specified otherswise to stderr
print STDERR "Perl script is being run.\n";
warn "This is printed to the log file\n";

# use Data::Dumper;
# print "The snakemake object is:\n" . Dumper($snakemake);

# Test auto input
while (<>) {
    print "Read line: $_";
}

# $snakemake is injected by Snakemake automatically
my @inputs = @{ $snakemake->{input}->{positional} }
  ;    # positional inputs as array, named ones are not included here
my $named_in = $snakemake->{input}->{named};    # named input with 'named' key

my $output   = $snakemake->{output}->{positional}->[0];
my $threads  = $snakemake->{threads};
my $rulename = $snakemake->{rulename};

print "Running rule: '$rulename' with $threads thread(s)\n";
print "Input: $named_in\n";
print "Output: $output\n";

use Test::Simple tests => 11;

# Test config values
ok(
    $snakemake->{config}->{test} == 1,
    "Config value 'test' is correctly passed"
);
ok(
    $snakemake->{config}->{testint} == 123,
    "Config value 'testint' is correctly passed"
);
ok(
    $snakemake->{config}->{testfloat} == 7.65432,
    "Config value 'testfloat' is correctly passed"
);
ok(
    $snakemake->{config}->{foodict}
      && ref( $snakemake->{config}->{foodict} ) eq "HASH"
      && $snakemake->{config}->{foodict}->{key0} eq "val0"
      && $snakemake->{config}->{foodict}->{key1} eq "val1",
    "Config dict 'foodict' is correctly passed as a hash reference"
);
ok(
    $snakemake->{config}->{'foo\' bar'} eq "let's go",
    "Config key/value with single quotes is correctly passed"
);

# Check wildcards
ok(
    $snakemake->{wildcards}->{rulename} eq "perl",
    "Wildcard 'rulename' is correctly passed"
);
ok(
    $snakemake->{wildcards}->{empty} eq "",
    "Wildcard 'empty' is correctly passed"
);

# Check params
ok(
    ref( $snakemake->{params}->{positional} ) eq "ARRAY"
      && $snakemake->{params}->{positional}->[0] eq "positional_param"
      && $snakemake->{params}->{positional}->[1] eq "second_positional_param",
    "Positional parameter is correctly passed"
);
ok(
    $snakemake->{params}->{integer} == 123,
    "Integer parameter is correctly passed"
);
ok(
    $snakemake->{params}->{astring} eq "foo\n'\\\" ",
    "String parameter is correctly passed"
);
ok(
    $snakemake->{params}->{alist}
      && ref( $snakemake->{params}->{alist} ) eq "ARRAY"
      && $snakemake->{params}->{alist}->[0] eq "a"
      && $snakemake->{params}->{alist}->[1] eq "b",
    "List parameter is correctly passed"
);
