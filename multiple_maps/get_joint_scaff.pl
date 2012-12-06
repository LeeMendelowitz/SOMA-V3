#!/usr/bin/perl

sub max {

    local($a, $b) = @_;
    return ($a > $b ? $a : $b);
}

sub abs {

    local($a) = @_;

    return ($a > 0 ? $a : -$a);
}

sub read_optical {

    local($dir, *opt_size, *opt_sd) = @_;
    
    $opt_data = `cat $dir/map\.opt`;
    @opt_size = ($opt_data =~ m/(\S+) \S+\n/g); 
    @opt_sd = ($opt_data =~ m/\S+ (\S+)\n/g); 
    map { $_ *= 1000 } @opt_size; map { $_ *= 1000 } @opt_sd; 
}

sub read_placement_info {

    local($dir, $rest, *ctg_info, *ctg_dir, *ctg_loc, *size) = @_;
    local(@data, $i);

    if($rest) {

	@data = split(/\n/, `cat $dir/seq\.rest`);
    }
    else {

	@data = split(/\n/, `cat $dir/seq\.unique_match $dir/seq\.placed`);
    }

    $j = 0;
    for($i = 0; $i < @data; $i += 5) {

	local($contig_id, $size, $forward, $opt_start, $opt_end) = split(/\s+/, $data[$i]);
	push @ctg_info, [$contig_id, $size, $forward, $opt_start, $opt_end, 
			 $data[$i+2], $data[$i+3], $data[$i+4]];

	$ctg_dir{$contig_id} = $forward;
	$size{$contig_id} = $size;
    }

    for($i = 0; $i <= $#ctg_info; $i++) {

	local($contig_id) = @{$ctg_info[$i]};
	$ctg_loc{$contig_id} = $i;
    }
}

sub get_locations {

    local($start_loc, $start_dir, *ctg_info, *opt_size, *opt_sd, *id, *dir, *loc, *sd) = @_;
    local($first, $gap) = (1, 0);
    local($loc) = ($start_loc); 
    local($i, $j);
    local($first_ctg_outhang, $first_site);

    $dir = $start_dir;

    do {

	local($ctg) = ($ctg_info[$loc]);

	local($contig_id, $size, $forward, $opt_start, $opt_end,
	      $opt_seq, $ctg_seq, $ctg_match_pos) = @$ctg;

	local(@ctg_matches, @new);
	@ctg_matches = split(/\s+/, $ctg_match_pos);

	@new = ($ctg_matches[0]);
	for($i = 1; $i < @ctg_matches; $i++) {

	    push @new, $ctg_matches[$i] if($ctg_matches[$i]-$ctg_matches[$i-1] >= 800);
	}

	@ctg_matches = @new;

	if(!$forward) {
	    @ctg_matches = map {$size-$_} (reverse @ctg_matches);	
	}

	local(@ctg_aln, $ctg_ind, $opt_aln, $opt_ind);
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
		
		$new_opt_start = ($opt_start+$opt_ind) % @opt_size if($new_opt_start eq "");
		$new_opt_end = ($opt_start+$opt_ind) % @opt_size;
	    }
	}

	if(!$first) {

	    ($sz, $sd) = (0, 0);

	    $i = $first_site;
	    $end = ($dir > 0 ? $new_opt_start : ($new_opt_end == 0 ? $#opt_size : $new_opt_end-1)); 

	    while($i != $end) {

		$sz += $opt_size[$i];
		#$sd += $opt_sd[$i]*$opt_sd[$i];
		$sd += 250 + $opt_size[$i]*$opt_size[$i]*0.0004;

		$i += $dir;
		$i = 0 if($i == @opt_size);
		$i = $#opt_size if($i == -1);
	    }

	    $sd = sqrt($sd);
	    $gap = max(0, $sz - ($dir > 0 ? $ctg_start : $ctg_end) - $first_ctg_outhang);
	}
	else {

	    $gap = -$size; $sd = 0;

	    if($dir > 0) {
		
		$first_ctg_outhang = $ctg_end; 
		$first_site = $new_opt_end; 
	    }
	    else {

		$first_ctg_outhang = $ctg_start; 
		$first_site = $new_opt_start-1;
	    }

	    $first_site = 0 if($first_site == @opt_size);
	    $first_site = $#opt_size if($first_site == -1);
	    $first = 0;
	}

	$forward = -1 if($forward != 1);

	$id[$loc] = $contig_id;
	$dir[$loc] = ($forward*$dir == 1 ? "BE" : "EB");
	$loc[$loc] = $gap;
	$sd[$loc] = $sd;

	$loc += $dir;
	$loc = 0 if($loc == @ctg_info);
	$loc = $#ctg_info if($loc == -1);

    } while($loc != $start_loc);    
}

if(@ARGV < 4) {

    print ("\nUSAGE: get_joint_scaff.pl <dir1> <dir2> <output> <anchor>\n\n".
	   "dir1\t= Directory containing SOMA output for first enzyme\n".
	   "dir2\t= Directory containing SOMA output for second enzyme\n".
	   "output\t= Directory to which the files \"joint.scaff\" and \"log\" will be written\n".
	   "anchor\t= Name of contig to use as \"anchor\" for merging placements. Set this to 1 to try all anchors.\n\n");
	
    exit;
}

($directory1, $directory2, $outdir, $tryall, $sd_threshold) = @ARGV;
$sd_threshold  = 12 if($sd_threshold eq "");

open(LOG, ">$outdir/log");
%size = ();

print "Step 1: Reading in reliable contig placements\n";

@opt_size1 = @opt_sd1 = @ctg_info1 = %ctg_dir1 = %ctg_loc1 = ();
&read_placement_info($directory1, 0, *ctg_info1, *ctg_dir1, *ctg_loc1, *size);
&read_optical($directory1, *opt_size1, *opt_sd1);

@opt_size2 = @opt_sd2 = @ctg_info2 = %ctg_dir2 = %ctg_loc2 = ();
&read_placement_info($directory2, 0, *ctg_info2, *ctg_dir2, *ctg_loc2, *size);
&read_optical($directory2, *opt_size2, *opt_sd2);

if($tryall eq "1") {

    map { push @list_ids, $_ if($ctg_loc2{$_} ne ""); } keys(%ctg_loc1);
    @list_ids = sort {$size{$b} <=> $size{$a}} @list_ids;
    print LOG "List of possible anchors: @list_ids\n";
}
else {

    @list_ids = ($tryall);
}

print "Step 2: Using an anchor to find a common orientation & order for placements\n";

foreach $anchor_id (@list_ids) {

    $dir1 = ($ctg_dir1{$anchor_id} == 1 ? 1 : -1); $start_loc1 = $ctg_loc1{$anchor_id};
    $dir2 = ($ctg_dir2{$anchor_id} == 1 ? 1 : -1); $start_loc2 = $ctg_loc2{$anchor_id};

    print LOG "Trying anchor: $anchor_id\n";
    print LOG "Orientations for maps: $dir1 $dir2\n";

    @id1 = @dir1 = @loc1 = @sd1 = ();
    &get_locations($start_loc1, $dir1, *ctg_info1, *opt_size1, *opt_sd1, *id1, *dir1, *loc1, *sd1);

    @id2 = @dir2 = @locs2 = @sd2 = ();
    &get_locations($start_loc2, $dir2, *ctg_info2, *opt_size2, *opt_sd2, *id2, *dir2, *loc2, *sd2);

    %dir1 = map { $id1[$_] => $dir1[$_] } (0..$#loc1);
    %dir2 = map { $id2[$_] => $dir2[$_] } (0..$#loc2);
    %loc1 = map { $id1[$_] => $loc1[$_] } (0..$#loc1);
    %loc2 = map { $id2[$_] => $loc2[$_] } (0..$#loc2);
    %sd1 = map { $id1[$_] => $sd1[$_] } (0..$#loc1);
    %sd2 = map { $id2[$_] => $sd2[$_] } (0..$#loc2);

    %loc = %sd = %dir = ();
    foreach $id (keys(%loc1)) {

	$loc{$id} = $loc1{$id};
	$sd{$id} = $sd1{$id};
	$dir{$id} = $dir1{$id};
    }

    $mismatched_directions = 0; $mismatched_locations = 0;
    foreach $id (keys(%loc2)) {

	if($loc{$id} eq "") {

	    $loc{$id} = $loc2{$id};
	    $sd{$id} = $sd2{$id};
	}
	else {

	    if($dir1{$id} ne $dir2{$id}) {

		print LOG "Directions don't match for $id"; $mismatched_directions++;
		if(&abs($loc1{$id}-$loc2{$id}) > $sd_threshold*($sd1{$id}+$sd2{$id})) {

		    print LOG " and locations don't match: $loc1{$id} $loc2{$id} $sd1{$id} $sd2{$id}\n";
		    $mismatched_locations++;
		}
		else {

		    print LOG " locations match $loc1{$id},$loc2{$id}\n";
		}

		delete($loc{$id});
	    }
	    else {

		if(&abs($loc1{$id}-$loc2{$id}) <= $sd_threshold*($sd1{$id}+$sd2{$id})) {

		    $loc{$id} = $loc2{$id};
		    $sd{$id} = $sd2{$id};
		}
		else {
		
		    print LOG "Locations don't match $id $loc1{$id} $loc2{$id} $sd1{$id} $sd2{$id}\n";
		    $mismatched_locations++; delete($loc{$id});
		}
	    }
	}

	$dir{$id} = $dir2{$id};
    }

    print "Step 3: Reading in possibly ambiguous placements\n" if($tryall ne "1");

    @new_ctg_info1 = ($ctg_info1[$start_loc1]); @new_ctg_dir1 = (); @new_ctg_loc1 = ();
    &read_placement_info($directory1, 1, *new_ctg_info1, *new_ctg_dir1, *new_ctg_loc1, *size);

    @new_ctg_info2 = ($ctg_info2[$start_loc2]); @new_ctg_dir2 = (); @new_ctg_loc2 = ();
    &read_placement_info($directory2, 1, *new_ctg_info2, *new_ctg_dir2, *new_ctg_loc2, *size);

    @new_id1 = @new_dir1 = @new_loc1 = @new_sd1 = ();
    &get_locations(0, $dir1, *new_ctg_info1, *opt_size1, *opt_sd1, *new_id1, *new_dir1, *new_loc1, *new_sd1);

    @new_id2 = @new_dir2 = @new_locs2 = @new_sd2 = ();
    &get_locations(0, $dir2, *new_ctg_info2, *opt_size2, *opt_sd2, *new_id2, *new_dir2, *new_loc2, *new_sd2);
    
    shift @new_id1; shift @new_dir1; shift @new_loc1; shift @new_sd1;   
    shift @new_id2; shift @new_dir2; shift @new_loc2; shift @new_sd2;   

    print "Step 4: Resolving ambiguous placements by merging placement information from the maps\n" if($tryall ne "1");

    for $i (0..$#new_id1) {
  
	for $j (grep {$new_id2[$_] eq $new_id1[$i]} (0..$#new_id2)) {

	    next if($loc{$new_id1[$i]} ne "");

	    if($new_dir1[$i] eq $new_dir2[$j] &&
	       abs($new_loc1[$i]-$new_loc2[$j]) <= $sd_threshold*($new_sd1[$i]+$new_sd2[$j])) {

		print LOG "Found matching placements $new_id1[$i] $new_dir1[$i] $new_loc1[$i],$new_loc2[$j] $loc{$new_id1[$i]}\n"; 

		$dir{$new_id1[$i]} = $new_dir1[$i];
		$loc{$new_id1[$i]} = $new_loc2[$j];
		$sd{$new_id1[$i]} = $new_sd2[$j];

		$dir1{$new_id1[$i]} = $new_dir1[$i];
		$loc1{$new_id1[$i]} = $new_loc1[$i];
		$sd1{$new_id1[$i]} = $new_sd1[$i];
		
		$dir2{$new_id2[$j]} = $new_dir2[$j];
		$loc2{$new_id2[$j]} = $new_loc2[$j];
		$sd2{$new_id2[$j]} = $new_sd2[$j];
	    }
	}
    }

    $output = ""; $first = 1; $gap = 0;
    $offset = 0; $total_gap = 0; $num_contigs = 0;
    $overlapping_placements = 0;
    foreach $id (sort {$loc{$a} <=> $loc{$b}} keys(%loc)) {

	if(!$first) {

	    if($loc{$id} < $loc{$last}+$size{$last}) {

		if($loc{$id}+$size{$id} > $last_end && $loc{$id}+$sd_threshold*($sd{$id}+$sd{$last}) > $last_end) {

		    #$loc{$id} = $last_end;
		}
		else {

		    printf LOG "Overlapping placement $id,$dir{$id},$size{$id},$sd{$id} at $loc1{$id},$loc2{$id}:$loc{$id} with $last,$last_end,%s,%s\n", 
		    ($loc1{$last} eq "" ? "" : "$loc1{$last}-".($loc1{$last}+$size{$last})), 
		    ($loc2{$last} eq "" ? "" : "$loc2{$last}-".($loc2{$last}+$size{$last}));
		    $overlapping_placements++; next;
		}
	    }

	    $gap = $loc{$id}-$last_end; $gap = 1 if($gap <= 0);
	    $output .= " $gap\n";
	}
	else {

	    $first = 0;
	}

	$last_end = $loc{$id}+$size{$id}; $last = $id;
	$output .= "$id $dir{$id} $size{$id}";
	$output .= " ".($loc1{$id} ne "" ? "*" : "").($loc2{$id} ne "" ? "+" : "") if($debug);
	$offset += $gap + $size{$id};
	$total_gap += $gap;
	$num_contigs++;
    }

    $bases = $offset - $total_gap;

    if($tryall ne "1") {
	
	open(OUT, ">$outdir/joint\.scaff");
	print OUT ">all_seq $num_contigs $bases $offset\n$output 0\n";
	close(OUT);
    }
    else {

	print "Placement using anchor: $anchor_id (size=$size{$anchor_id})\n";
	print "Contigs=$num_contigs, Bases=$bases, ";
	print "Direction-Mismatch=$mismatched_directions, Location-Mismatch=$mismatched_locations, ";
	print "Overlaps=$overlapping_placements\n"; 
    }
}

