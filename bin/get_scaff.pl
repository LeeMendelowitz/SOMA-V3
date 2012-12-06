#!/usr/bin/perl

sub max {

    ($a, $b) = @_;
    return ($a > $b ? $a : $b);
}


($input_file, $map_file) = @ARGV;
open(IN, $input_file);
@data = <IN>;

$opt_data = `cat $map_file`;
@opt_size = ($opt_data =~ m/(\S+) \S+\n/g); 
@opt_sd = ($opt_data =~ m/\S+ (\S+)\n/g); 
map { $_ *= 1000 } @opt_size; map { $_ *= 1000 } @opt_sd; 

for($i = 0; $i < @data; $i += 5) {

    ($contig_id, $size, $forward, $opt_start, $opt_end, $type) = split(/\s+/, $data[$i]);
    push @ctg_info, [$contig_id, $size, $forward, $opt_start, $opt_end, $type, 
		     $data[$i+2], $data[$i+3], $data[$i+4]];
}

@ctg_info = sort {$a->[3] <=> $b->[3]} @ctg_info;    

$scaff_num = 0; $all_scaffs = "";
$first = 1; $gap = 0;

foreach $ctg (@ctg_info) {

    ($contig_id, $size, $forward, $opt_start, $opt_end, $type,
     $opt_seq, $ctg_seq, $ctg_match_pos) = @$ctg;

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

    if(!$first) {

	$sz = 0; $sd = 0;
	foreach $i ($last_end_site..$new_opt_start-1) {

	    $sz += $opt_size[$i];
	    $sd += $opt_sd[$i]*$opt_sd[$i];

	    if($opt_sd[$i] == 0) {

		$bases = $offset - $total_gap;
		$all_scaffs .= ">scaff$scaff_num $num_contigs $bases $offset\n$output 0 0 $last_type\n"; 
		$scaff_num++; $offset = $total_gap = $num_contigs = 0; $output = "";
		goto REST;
	    }		
	}

	$sd = sqrt($sd);
	$gap = max(1, $sz - $last_end_off - $ctg_start);
	$output .= sprintf(" %d %d $last_type\n", $gap, $sd);
    }
    else {

	$first = 0;
    }

REST:

    $output .= "$contig_id ".($forward ? "BE" : "EB"). " $size";
    $last_type = $type;
    $last_end_site = $new_opt_end; $last_end_off = $ctg_end; 
    $offset += $gap + $size;
    $total_gap += $gap;
    $num_contigs++;
}

$bases = $offset - $total_gap;
$all_scaffs .= ">scaff$scaff_num $num_contigs $bases $offset\n$output 0 0 $last_type\n"; 
print $all_scaffs;
