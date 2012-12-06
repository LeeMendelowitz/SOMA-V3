#!/usr/bin/perl

($summary_file, $output_file) = @ARGV;
open(IN, $summary_file);
open(OUT, ">$output_file");

while(<IN>) {

    if($_ =~ m/\d+\s+([\d\.]+)\s+([\d\.]+)/) { 

	print OUT "$1 $2\n"; 
    }
    if($_ =~ m/S=\"(\d+)\" STDDEV=\"([\d\.]+)\"/) {

	printf OUT "%.3f $2\n", $1/1e3;
    }

    $done = 1 if($_ =~ m/<\/RESTRICTION_MAP>/);
    print OUT "100.000 0.000\n" if($done && $_ =~ m/<RESTRICTION_MAP/);
}

close(OUT);
close(IN);
