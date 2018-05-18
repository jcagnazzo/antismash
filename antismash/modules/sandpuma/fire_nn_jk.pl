#!/bin/env perl

use strict;
use warnings;

my @alpha = (10, 1, 0.1, 0.01, 0.001);
foreach my $testfaa (glob("*/*/fullset_smiles_knife*.faa")){
    my @p = split(/\//, $testfaa);
    my $shuffle = $p[-2];
    my $k = $testfaa;
    $k =~ s/^.+knife(\d+).+$/k$1_as_query/;
    my $trainfaa = join('/', @p[0..$#p-1], $k.'_train.faa');
    my $trainaln = $trainfaa;
    $trainaln =~ s/\.faa$/\.afa/;
    #unless(-e $trainaln){
	#system("mafft --quiet --namelength 75 $trainfaa > $trainaln");
    #}
    my $trainnwk = $trainfaa;
    $trainnwk =~ s/\.faa$/\.nwk/;
    #unless(-e $trainnwk){
	#system("fasttree -log $trainnwk.log < $trainaln > $trainnwk")
    #}
    my $trainrefpkg = $trainfaa;
    $trainrefpkg =~ s/\.faa$/\.refpkg/;
    #unless(-d $trainrefpkg){
	#system("taxit create --aln-fasta $trainaln --tree-stats $trainnwk.log --tree-file $trainnwk -P $trainrefpkg -l a_domain");
    #}
    #print join("\t", $shuffle, $k, $testfaa, $trainfaa, $trainaln, $trainnwk, $trainrefpkg)."\n";
    #print join("\t", $shuffle, $k)."\n";
    foreach my $a (@alpha){
	foreach(1..10){
	    system("python sandpuma2.py $shuffle $k $testfaa $trainfaa $trainaln $trainnwk $trainrefpkg $a");
	}
    }
}
