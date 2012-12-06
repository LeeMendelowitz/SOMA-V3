#!/usr/bin/perl

sub max {

    ($a, $b) = @_;
    return ($a > $b ? $a : $b);
}

$max_sigma = 12;

sub get_protruding_ends {

    local($forward, $size, $opt_start, $opt_end, $opt_seq, $ctg_seq, $ctg_match_pos, *opt_size, *opt_sd) = @_;
    local(@ctg_matches, @new, $i, $j, @ctg_aln, @opt_aln);

    @ctg_matches = split(/\s+/, $ctg_match_pos);

    @new = ($ctg_matches[0]);
    for($i = 1; $i < @ctg_matches; $i++) {

	push @new, $ctg_matches[$i] if($ctg_matches[$i]-$ctg_matches[$i-1] >= 700);
    }

    @ctg_matches = @new;

    if(!$forward) {
	@ctg_matches = map {$size-$_} (reverse @ctg_matches);	
    }

    $ctg_seq =~ s/;/ ;/g; $ctg_seq =~ s/\s+\n//; @ctg_aln = split(/\s+/, $ctg_seq);

    $j = 0; $ctg_ind = 0; $ctg_start = "";
    for(; $j < @ctg_aln; $j++) {

	if($ctg_aln[$j] ne ";") {
	    
	    $ctg_ind++;
	}
	else {

	    $ctg_start = $ctg_matches[$ctg_ind-1] if($ctg_start eq "");
	    $ctg_end = $size - $ctg_matches[$ctg_ind-1];
	}
    }

    $opt_seq =~ s/;/ ;/g; $opt_seq =~ s/\s+\n//; @opt_aln = split(/\s+/, $opt_seq);

    $j = 0; $opt_ind = 0; $new_opt_start = "";
    for(; $j < @opt_aln; $j++) {

	if($opt_aln[$j] ne ";") {
	    $opt_ind++;
	}
	else {
	    
	    $new_opt_start = $opt_start+$opt_ind if($new_opt_start eq "");
	    $new_opt_end = $opt_start+$opt_ind;
	}
    }

    $opt_start_length = 0; $opt_start_sd = 0;
    for($i = $opt_start+1; $i < $new_opt_start; $i++) {
	
	$opt_start_length += $opt_size[$i%scalar(@opt_sd)];
	$opt_start_sd += $opt_sd[$i%scalar(@opt_sd)]**2;
    }

    $opt_end_length = 0; $opt_end_sd = 0;
    for($i = $new_opt_end; $i < $opt_end; $i++) {
	
	$opt_end_length += $opt_size[$i%scalar(@opt_sd)];
	$opt_end_sd += $opt_sd[$i%scalar(@opt_sd)]**2;
    }

    return (max(int($ctg_start - ($opt_start_length + $max_sigma*sqrt($opt_start_sd))), 0), 
	    max(int($ctg_end - ($opt_end_length + $max_sigma*sqrt($opt_end_sd))), 0));
}


($dir) = @ARGV;

$opt_data = `cat $dir/map\.opt`;
@opt_size = ($opt_data =~ m/(\S+) \S+\n/g); 
@opt_sd = ($opt_data =~ m/\S+ (\S+)\n/g); 
map { $_ *= 1000 } @opt_size; map { $_ *= 1000 } @opt_sd; 

@data = split(/\n/, `cat $dir/seq\.unique_match`);
for($i = 0; $i < @data; $i += 5) {

    ($contig_id, $size, $forward, $opt_start, $opt_end) = split(/\s+/, $data[$i]);
    ($opt_seq, $ctg_seq, $ctg_match_pos) = ($data[$i+2], $data[$i+3], $data[$i+4]);

    ($start_hang, $end_hang) = get_protruding_ends($forward, $size, $opt_start, $opt_end,
						   $opt_seq, $ctg_seq, $ctg_match_pos, *opt_size, *opt_sd);

    $opt_size[$opt_start] = max($opt_size[$opt_start]-$start_hang, 0);
    $opt_size[$opt_end%scalar(@opt_sd)] = max($opt_size[$opt_end%scalar(@opt_sd)]-$end_hang, 0);
    
    for($j = $opt_start+1; $j < $opt_end; $j++) {

	$opt_sd[$j%scalar(@opt_sd)] = 0;
    } 
}

@matches = split(/\n/, `cat $dir/seq\.match`);

open(PLACED, ">$dir/seq\.placed");
$first_time = 1;

do {

    %copies = 0; $placed = 0;
    
    for($i = 0; $i < @matches; $i += 5) {
	
	next if($matches[$i] eq "");

	($contig_id, $size, $forward, $opt_start, $opt_end) = split(/\s+/, $matches[$i]);
	($opt_seq, $ctg_seq, $ctg_match_pos) = ($matches[$i+2], $matches[$i+3], $matches[$i+4]);
	
	($start_hang, $end_hang) = get_protruding_ends($forward, $size, $opt_start, $opt_end,
						       $opt_seq, $ctg_seq, $ctg_match_pos, *opt_size, *opt_sd);
	
	$conflict = 0;
	for($j = $opt_start; $j <= $opt_end; $j++) {
	    
	    $conflict = 1 if($opt_sd[$j%scalar(@opt_sd)] == 0);
	}

	$conflict = 1 if($opt_size[$opt_start]+$max_sigma*$opt_sd[$opt_start] < $start_hang);
	$conflict = 1 if($opt_size[$opt_end%scalar(@opt_sd)]+$max_sigma*$opt_sd[$opt_end%scalar(@opt_sd)] < $end_hang);
	
	$matches[$i] = "" if($conflict);
	$copies{$contig_id}++ if(!$conflict);
    }
    $first_time = 0;

    for($i = 0; $i < @matches; $i += 5) {
	
	next if($matches[$i] eq "");

	($contig_id, $size, $forward, $opt_start, $opt_end) = split(/\s+/, $matches[$i]);
	($opt_seq, $ctg_seq, $ctg_match_pos) = ($matches[$i+2], $matches[$i+3], $matches[$i+4]);
	
	if($copies{$contig_id} == 1) {
	    
	    print PLACED "$matches[$i]\n$matches[$i+1]\n$matches[$i+2]\n$matches[$i+3]\n$matches[$i+4]\n";
	    $matches[$i] = "";
	    
	    ($start_hang, $end_hang) = get_protruding_ends($forward, $size, $opt_start, $opt_end,
							   $opt_seq, $ctg_seq, $ctg_match_pos, *opt_size, *opt_sd);
	    
	    $opt_size[$opt_start] = max($opt_size[$opt_start]-$start_hang, 0);
	    $opt_size[$opt_end%scalar(@opt_sd)] = max($opt_size[$opt_end%scalar(@opt_sd)]-$end_hang, 0);
	    
	    for($j = $opt_start+1; $j < $opt_end; $j++) {
		
		$opt_sd[$j%scalar(@opt_sd)] = 0;
	    } 

	    $placed = 1;
	    last;
	}
    }

} while($placed);
close(PLACED);

open(REST, ">$dir/seq\.rest");
for($i = 0; $i < @matches; $i += 5) {
	
    next if($matches[$i] eq "");
    print REST "$matches[$i]\n$matches[$i+1]\n$matches[$i+2]\n$matches[$i+3]\n$matches[$i+4]\n";
}
close(REST);
