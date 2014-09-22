#!/usr/bin/env perl

use strict;
use warnings;

open (my $fh, "GSE29087_L139_expression_tab.txt") or die $!;
my $header = <$fh>;
chomp $header;

my @fields = split(/\t/, $header);

foreach my $field (@fields) {
    print "Field: [$field]\n";
    $field =~ s/\s//g;
    if ($field =~ /^[A-Z](\d+)/) {
        my $num = $1;
        if ($num >= 7) {
            $field = "MEF_$field";
        }
        else {
            $field = "ES_$field";
        }
    }
}

print join("\t", @fields) . "\n";
while (<$fh>) {
    print;
}

exit(0);
