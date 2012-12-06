#!/usr/bin/perl

($dir) = @ARGV;

$file_prefix = "$dir/seq";
if(-e "$file_prefix\.exhaustive") {

    $schedule_file = "$file_prefix\.exhaustive";
}
else {

    $schedule_file = "$file_prefix\.greedy";
}

$data = `cat $file_prefix\.unique_match`; $data =~ s/([^\d\s]\S+ \d+ \d+ \d+ \d+)/\1 MATCH/g; `echo -n \"$data\" > $file_prefix\.all_placed`;
$data = `cat $file_prefix\.placed`; $data =~ s/([^\d\s]\S+ \d+ \d+ \d+ \d+)/\1 FILTER/g; `echo -n \"$data\" >> $file_prefix\.all_placed`;
$data = `cat $schedule_file`; $data =~ s/([^\d\s]\S+ \d+ \d+ \d+ \d+)/\1 SCHEDULE/g; `echo -n \"$data\" >> $file_prefix\.all_placed`;

`./make_fig.pl $file_prefix\.all_placed $dir/map\.opt > $file_prefix\.fig`;
`cat $file_prefix\.fig | fig2dev -L ps > $file_prefix\.ps`;
`./get_scaff.pl $file_prefix\.all_placed $dir/map\.opt > $file_prefix\.scaff`;
`./scaff2agp.pl $file_prefix\.scaff $file_prefix\.agp`;

$summary = `cat $file_prefix\.match-summary`;
$summary =~ m/\>1: (\d+) Placed: (\d+)/;
$total = $1; $prefig = $2;
$placed_count = `wc -l $file_prefix\.placed`; chomp $placed_count; $placed_count /= 5;
$scheduled_count = `wc -l $schedule_file`; chomp $scheduled_count; $scheduled_count /= 5;
$scaff = `cat $file_prefix\.scaff`; @lengths = ($scaff =~ m/\>\S+\s+\S+\s+(\d+)/g); $length = 0; map {$length += $_} @lengths; $length /= 1e6;

$opt = `cat $dir/map\.opt`; $opt =~ s/100.000 0.000\s*\n//g; $opt =~ s/\s+\S+\n/\n/g; @opt = split(/\s+/, $opt);
map {$opt_size += $_ } @opt; $ratio = $length*1e3/$opt_size;
($dummy, $count) = split(/\s+/, `grep \">\" $file_prefix\.fa | wc`);

open(SUMMARY, ">$file_prefix\.summary");
print SUMMARY ("Summary of results\n".('-'x30)."\n\n");
print SUMMARY "Number of sequences: $count\n";
print SUMMARY "Number of sequences with > 1 site: $total\n";
print SUMMARY "Sequences placed by matching algorithm: $prefig\n";
print SUMMARY "Sequences placed by matching-filtering: $placed_count\n";
print SUMMARY "Sequences placed by scheduling algorithm: $scheduled_count\n";
printf SUMMARY "Total length of placed sequences: %.2f Mbp\n", $length;
printf SUMMARY "(Total length)/(Size of optical map): %.2f\n\n", $ratio;
close(SUMMARY);     
