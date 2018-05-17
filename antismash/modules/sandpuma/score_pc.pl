#!/bin/env perl

use strict;
use warnings;

my $jqf = shift;
my %j = ();
open my $jfh, '<', $jqf or die $!;
while(<$jfh>){
	chomp;
	next if($_ =~ m/^Shuffle/);
	my ($shuf, $kasq, $query, $true, $pc, $fpc, $snn) = split(/\t/, $_);
	$fpc = $pc if($fpc eq 'no_force_needed');
	$j{$shuf}{$kasq}{$query}{'prediCAT'}{'call'} = $pc;
	if($pc eq 'no_confident_result'){
		$j{$shuf}{$kasq}{$query}{'prediCAT'}{'cov'} = 'N';
	}else{
		$j{$shuf}{$kasq}{$query}{'prediCAT'}{'cov'} = 'Y';
	}
	$j{$shuf}{$kasq}{$query}{'forced_prediCAT'}{'call'} = $fpc;
	if($fpc eq 'no_confident_result'){
	    $j{$shuf}{$kasq}{$query}{'forced_prediCAT'}{'cov'} = 'N';
	}else{
	    $j{$shuf}{$kasq}{$query}{'forced_prediCAT'}{'cov'} = 'Y';
	}
	$j{$shuf}{$kasq}{$query}{'snn'} = $snn;
}
close $jfh;
## pid
my %p = ();
open my $pfh, '<', 'pid.tsv' or die $!;
while(<$pfh>){
	chomp;
	my ($s, $k, $q, $pid) = split(/\t/, $_);
	$q =~ s/[\[\]\.\(\)]/-/g;
	$p{$s}{$k}{$q} = $pid;
}
close $pfh;

print join("\t", 'shuffle', 'jackknife', 'query', 'pid', 'spec', 'called_spec', 'method', 'call_made', 'call', 'snn', 'methshuf') . "\n";
foreach my $shuf (keys %j){
	foreach my $kasq (keys %{$j{$shuf}}){
		foreach my $q (keys %{$j{$shuf}{$kasq}}){
			my @id = split(/_+/, $q);
			foreach my $meth ('prediCAT', 'forced_prediCAT'){
				my $call = 'N';
				if($j{$shuf}{$kasq}{$q}{$meth}{'cov'} eq 'Y'){
					my @true = ();
					my @ca = ();
					if($id[-1] =~ m/\|/){
						@true = split(/\|/, $id[-1]);
					}else{
						$true[0] = $id[-1];
					}
					my $ct = $j{$shuf}{$kasq}{$q}{$meth}{'call'};
					$ct =~ s/(\S+)_.+/$1/;
					## Error correct
					if($ct eq 'abu-iva' || $ct eq 'athr-thr' || $ct eq 'dhab-dht' || $ct eq 'dhb-sal' || $ct eq 'dpg-dhpg' || $ct eq 'hpg-hpg2cl' || $ct eq ''){
						$ct =~ s/-/\|/g;
					}elsif($ct eq 'pro-me-pro'){
						$ct = 'pro|me-pro';
					}
					if($ct =~ m/\|/){
						@ca = split(/\|/, $ct);
					}else{
						$ca[0] = $ct;
					}
					$j{$shuf}{$kasq}{$q}{$meth}{'call'} = $ct;
					foreach my $t (@true){
						foreach my $c (@ca){
							$call = 'Y' if($t eq $c);
						}
					}
				}
				print join("\t", $shuf, $kasq, $q);
				print "\t".$p{$shuf}{$kasq}{$q};
				print "\t".$id[-1];
				print "\t".$j{$shuf}{$kasq}{$q}{$meth}{'call'};
				print "\t".$meth;
				print "\t".$j{$shuf}{$kasq}{$q}{$meth}{'cov'};
				print "\t".$call;
				print "\t".$j{$shuf}{$kasq}{$q}{'snn'};
				print "\t".$meth.'_'.$shuf."\n";
			}
		}
	}
}
