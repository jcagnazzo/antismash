#!/bin/env perl

use strict;
use warnings;

my %c = ();
while(<>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$c{$1} += 1;
    }
}
foreach my $s (sort {$c{$a} <=> $c{$b}} keys %c){
    print join("\t", $c{$s}, $s)."\n";
}
