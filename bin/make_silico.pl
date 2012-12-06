#!/usr/bin/perl

#
# The second argument can also be an enzyme name that
# can be found in repbase.staden
#
($infile, $restriction_seq, $outfile) = @ARGV;

if($match = `grep $restriction_seq repbase.staden`) {

    ($name, $restriction_seq) = split(/\//, $match);
    $restriction_seq =~ s/\'//g;
}

print "ERROR: make_silico.pl was given the argument $restriction_seq and couldn't find the enzyme OR\n". 
      "       was given a degenerate restriction sequence\n" if($restriction_seq =~ m/[^acgtACGT]/);

$restriction_seq = uc($restriction_seq);

open(IN, $infile);
$lines = join("", <IN>);
close(IN);

open(OUT, ">$outfile");
while($lines =~ m/\>(\S+).*\n([^\>]+)/g) {

    $name = $1; $seq = $2; $seq =~ s/\n//g; $seq = uc($seq);
    
    @positions = (); $pos = 0;
    while(($pos = index($seq, $restriction_seq, $pos+1)) != -1) {

	push @positions, ($pos+1);
    }

    print OUT ("$name ".length($seq)." ".@positions."\n@positions \n") if(@positions != 0);
}

close(OUT);
