#!/usr/bin/perl

($infile, $outfile) = @ARGV;

$scaff = `cat $infile`;

open(OUT, ">$outfile");
while($scaff =~ m/\>(\S+).*\n([^>]+)/g) {

    $scaff_name = $1; @contigs = split(/\n/, $2);
    $scaff_pos = 1; $part_no = 1;

    foreach $contig (@contigs) {

	($name, $or, $size, $offset) = split(/\s+/, $contig);
	if($name =~ s/\{(\d+),(\d+)\}//) {

	    ($start, $end) = sort($1, $2);
	}
	else {
	    
	    $start = 1; $end = $size;
	}

	
	printf OUT "$scaff_name\t$scaff_pos\t%d\t%d\tW\t$name\t$start\t$end\t%s\n", ($scaff_pos += $end-$start), $part_no++, 
	($or eq "BE" ? "+" : "-");
	$scaff_pos++;

	if($offset > 0) {

	    printf OUT "$scaff_name\t$scaff_pos\t%d\t%d\tN\t$offset\tfragment\tNo\t\n", ($scaff_pos += $offset-1), $part_no++;
	    $scaff_pos++;
	}
    }
}
close(OUT);
