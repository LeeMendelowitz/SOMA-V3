#!/usr/bin/perl

use xfig;

($input_file, $map_file) = @ARGV;
open(IN, $input_file);
@data = <IN>;
close(IN);

print_xfig_header(1200);

$min_frag_size = 800;

$opt_y = 5000;
$x_offset = 1000;
$tick_size = 100;
$scale = 10;

$lbl_rows = 2;
$cur_lbl_row = 0;
$lbl_off = 300;

$pos = $x_offset; 

$opt_data = `cat $map_file`;
@opt_size = ($opt_data =~ m/(\S+) \S+\n/g); 
@opt_sd = ($opt_data =~ m/\S+ (\S+)\n/g); 
map { $_ *= 1000 } @opt_size; map { $_ *= 1000 } @opt_sd; 

for($i = 0; $i < $#opt_size; $i++) {

    ($size, $sd) = ($opt_size[$i], $opt_sd[$i]);
    $sc_size = int($size/$scale);
	
    print_horiz_text($pos+int($sc_size/2), $opt_y+($cur_lbl_row+1)*$lbl_off, 
		     sprintf("%.1f,%.1f", $size/1000, $sd/1000), 0);

    $opt_pos[$i] = ($pos+=$sc_size);

    print_tick($pos, $opt_y, 0, $tick_size);
    $cur_lbl_row = ++$cur_lbl_row % $lbl_rows; 
}

($size, $sd) = ($opt_size[$i], $opt_sd[$i]);
$sc_size = int($size/$scale);
print_horiz_text($pos+int($sc_size/2), $opt_y+($cur_lbl_row+1)*$lbl_off, 
		 sprintf("%.1f,%.1f", $size/1000, $sd/1000), 0);
$pos+=$sc_size;

print_solid_line($x_offset, $opt_y, $pos, $opt_y);

for($i = 0; $i < @data; $i += 5) {

    ($contig_id, $size, $forward, $opt_start, $opt_end) = split(/ /, $data[$i]);
    push @ctg_info, [$contig_id, $size, $forward, $opt_start, $opt_end, 
		     $data[$i+2], $data[$i+3], $data[$i+4]];
}

@ctg_info = sort {$a->[3] <=> $b->[3]} @ctg_info;    

$ctg_rows = 3;
$cur_ctg_row = 0;
$ctg_off = 600;
$ctg_lbl_off = 200;

foreach $ctg (@ctg_info) {

    ($contig_id, $size, $forward, $opt_start, $opt_end,
     $opt_seq, $ctg_seq, $ctg_match_pos) = @$ctg;

    $ctg_match_pos =~ s/\s+\n//;
    @ctg_matches = split(/\s+/, $ctg_match_pos);

    @new = ($ctg_matches[0]);
    for($i = 1; $i < @ctg_matches; $i++) {

	push @new, $ctg_matches[$i] if($ctg_matches[$i]-$ctg_matches[$i-1] >= $min_frag_size);
	$size += $ctg_matches[$i]-$ctg_matches[$i-1] if($ctg_matches[$i]-$ctg_matches[$i-1] < $min_frag_size);
    }

    @ctg_matches = @new;

    if(!$forward) {
	@ctg_matches = map {$size-$_} (reverse @ctg_matches);	
    }

    $opt_seq =~ s/;/ ;/g; $opt_seq =~ s/\s+\n//; @opt_aln = split(/\s+/, $opt_seq);
    $ctg_seq =~ s/;/ ;/g; $opt_seq =~ s/\s+\n//; @ctg_aln = split(/\s+/, $ctg_seq);

    $j = 0; $opt_ind = 0; $ctg_ind = 0; $ctg_start = "";
    for($i = 0; $i < @opt_aln; $i++) {
	if($opt_aln[$i] ne ";") {
	    $opt_ind++;
	}
	else {

	    for(; $j < @ctg_aln; $j++) {

		if($ctg_aln[$j] ne ";") {
		    $ctg_ind++;
		}
		else {
	
		    $ctg_start = $opt_pos[$opt_start+$opt_ind-1] - int($ctg_matches[$ctg_ind-1]/$scale) if($ctg_start eq ""); 		    
		    $j++;

		    if($opt_start+$opt_ind-1 >= @opt_pos) {

			if($opt_start+$opt_ind-1 == @opt_pos) {

			    $opt_x = $x_offset;
			}
			else {

			    $opt_x = $opt_pos[(($opt_start+$opt_ind-1) % @opt_pos)-1];
			}

			$ctg_x = $ctg_start+int($ctg_matches[$ctg_ind-1]/$scale);
			$ctg_y = $opt_y-($cur_ctg_row+1)*$ctg_off;
			print_dotted_line($opt_x, $opt_y, $opt_x, $opt_y-($ctg_rows+1)*$ctg_off);
			print_dotted_line($opt_x, $opt_y-($ctg_rows+1)*$ctg_off, $x_offset, $opt_y-($ctg_rows+3)*$ctg_off);
			print_dotted_line($ctg_x, $ctg_y, $ctg_x, $opt_y-($ctg_rows+1)*$ctg_off);
			print_dotted_line($ctg_x, $opt_y-($ctg_rows+1)*$ctg_off, $ctg_start+int($size/$scale), $opt_y-($ctg_rows+3)*$ctg_off);	      
			last; 
		    }

		    print_dotted_line($ctg_start+int($ctg_matches[$ctg_ind-1]/$scale),
				      $opt_y-($cur_ctg_row+1)*$ctg_off, 
				      $opt_pos[$opt_start+$opt_ind-1], $opt_y);
		    last;
		}
	    }
	}
    }

    map { print_tick($ctg_start+int($_/$scale), 
		     $opt_y-($cur_ctg_row+1)*$ctg_off, 1, $tick_size) } @ctg_matches;
	
    print_horiz_text($ctg_start, $opt_y-($cur_ctg_row+1)*$ctg_off-$ctg_lbl_off, $contig_id, 0);

    print_solid_line($ctg_start, $opt_y-($cur_ctg_row+1)*$ctg_off, 
		     $ctg_start+int($size/$scale), $opt_y-($cur_ctg_row+1)*$ctg_off);

    if ($forward) {

	print_circle($ctg_start, $opt_y-($cur_ctg_row+1)*$ctg_off, 1, 3);
	print_diamond($ctg_start+int($size/$scale), $opt_y-($cur_ctg_row+1)*$ctg_off, 1, 2);
    }
    else {

	print_diamond($ctg_start, $opt_y-($cur_ctg_row+1)*$ctg_off, 1, 2);
	print_circle($ctg_start+int($size/$scale), $opt_y-($cur_ctg_row+1)*$ctg_off, 1, 3);
    }
    
    $cur_ctg_row = ++$cur_ctg_row % $ctg_rows;  
}
