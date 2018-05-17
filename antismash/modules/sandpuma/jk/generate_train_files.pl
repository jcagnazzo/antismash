#!/bin/env perl

use strict;
use warnings;

foreach my $testfaa (glob("*/fullset_smiles_knife*.faa")){
    my $shuffle = $testfaa;
    $shuffle =~ s/^(jk\d+).+/$1/;
    my $k = $testfaa;
    $k =~ s/^.+knife(\d+).+$/k$1_as_query/;
    my $trainfaa = $shuffle.'/'.$k.'_train.faa';
    my $trainaln = $trainfaa;
    $trainaln =~ s/\.faa$/\.afa/;
    unless(-e $trainaln){
	system("mafft --quiet --namelength 75 $trainfaa > $trainaln");
    }
    my $trainnwk = $trainfaa;
    $trainnwk =~ s/\.faa$/\.nwk/;
    unless(-e $trainnwk){
	system("fasttree -log $trainnwk.log < $trainaln > $trainnwk")
    }
    my $trainrefpkg = $trainfaa;
    $trainrefpkg =~ s/\.faa$/\.refpkg/;
    unless(-d $trainrefpkg){
	system("taxit create --aln-fasta $trainaln --tree-stats $trainnwk.log --tree-file $trainnwk -P $trainrefpkg -l a_domain");
    }
    print join("\t", $shuffle, $k, $testfaa, $trainfaa, $trainaln, $trainnwk, $trainrefpkg)."\n";
}
