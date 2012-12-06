#!/usr/bin/perl

if(@ARGV < 3) {

    print "\nUSAGE: soma.pl <fasta-file> <map-summary> <output-dir>\n\n";
    print "fasta-file\t= Contig/Scaffold sequences in multi-fasta format\n";
    print "map-summary\t= Summary of genome-wide optical map\n";
    print "output-dir\t= Directory where the input, log and output information will be written\n\n";
    exit;
}

($seqfile, $mapfile, $dir, $circular_flag) = @ARGV;

$MSG = STDOUT;

#
# Extracting the enzyme name.
#
if(($first_line = `head -1 $mapfile`) =~ m/xml/) {

    $map = `cat $mapfile`; $map =~ m/ENZYME=\"([^\"]+)\"/;
    $enzyme = $1;
}
else {

#
# Assumption: Third field of the header in the map file has
# the enzyme name.
#
    ($genome, $dummy, $enzyme) = split(/\s+/, `head -1 $mapfile`);
}

if($circular_flag eq "") {

    $circular_flag = 1;
    $circular_flag = 0 if(`cat $mapfile` =~ m/CIRCULAR=\"false\"/);
}

print $MSG "Step 1: Creating input files (Enzyme: $enzyme, Circular: $circular_flag)\n";

$seq_prefix = "$dir/seq";
$map_prefix = "$dir/map";

`mkdir $dir` if(!(-e $dir));
`cp -f $seqfile $seq_prefix\.fa`;
`cp -f $mapfile $map_prefix\.txt`;
open(LOG, ">$dir/log");
print LOG "Step 1: Creating input files\n";

print `./make_silico.pl $seq_prefix\.fa $enzyme $seq_prefix\.silico`;
print LOG `python make_opt.py $map_prefix\.txt $map_prefix\.opt`;

print $MSG "Step 2: Matching sequences to the map\n";
print LOG "Step 2: Matching sequences to the map\n";
print LOG `./match -s $seq_prefix\.silico -m $map_prefix\.opt -c $circular_flag -p 0 -f 0.01 -o $seq_prefix`;

print $MSG "Step 3: Resolving unique placements\n";
print LOG "Step 3: Resolving unique placements\n";
print LOG `./place_rest.pl $dir`;

print $MSG "Step 4: Finding optimal placement for the rest\n";
print LOG "Step 4: Finding optimal placement for the rest\n";
print LOG `./schedule $dir 0`;

print $MSG "Step 5: Collecting results\n";
print LOG "Step 5: Collecting results\n";
print LOG `./collect_results.pl $dir`;

