#!/usr/bin/perl -w
### Writtern by Albert K Tai, PhD
### Research Assistant Professor, Department of Integrative Phyisology and Pathobiology
### Tufts University School of Medicine, Boston, MA

use strict;
use Math::Round;

### Sequence information 

#my $ad_seq = "CTGTAGGCACCATCAATCGTTACGTAG";
my $ad_seq = "CTGTAGGCACCATCAATCGTTAC";
#my $TERC_seq = "AGGCCGCAGGAAGAGGAACGGAGCGAGTCCCCGCGCGCGGCGCGAT"."TCCCTGAGCTGTGGGACGTGCACCCAGGACTCGGCTCACACATGC";
my $TERC_seq = "TCCCTGAGCTGTGGGACGTGCACCCAGGACTCGGCTCACACATGC";
#my $TERCjunc = "TCCCTGAGCTGTGGGACGTGCACCC AGGACTCGGCTCACACATGCAGTTCGCTTTCCTGTTGGTGGGGGGAACGCCGATCGTGCGCATCCGTCACCCCTCGCCGGCA
#ATGGGGGCTTGTGAACCCCCAAACCTGACTGACTGGGCCAGTGTGCTGCAAATTGGCAGGAGACGTGAAGGCACCTCCAAAGTCGGCCAAAA";
my $TERCjunc = "CACCCAGGACTCGGCTCACACATGCAGTTCGCTTTCCTGTTGGTGGGGGGAACGCCGATCGTGCGCATCCGTCACCCCTCGCCGGCA
ATGGGGGCTTGTGAACCCCCAAACCTGACTGACTGGGCCAGTGTGCTGCAAATTGGCAGGAGACGTGAAGGCACCTCCAAAGTCGGCCAAAA";

## Open input SAM file which was passed as an arguement when the perl script is invoked.
unless(open(SAM, "$ARGV[0]")){die;}
print "Opening $ARGV[0] for processing..\n";


## Senitize and reformatting input file name for output
(my $filename = $ARGV[0]) =~ m/.sam$/;
$filename =~ s/.sam/.txt/;
$filename =~ s/trim_25_25_//;
$filename = "tail_counted_results_".$filename;
unless(open(OUT, ">$filename")){die;}
$filename =~ s/.txt/.log/;
unless(open(LOG, ">$filename")){die;}
print LOG "Opening $ARGV[0] for processing..\n";

### This part of the script aims to "collapse" information from SAM file(s)
### by merging reads that have the same chromosome, mapped coordinate,
### primary cigar code and DNA sequence. The use of multiple fields to compare
### is likely reductant and probably does not add much to the stringency of using
### sequence alone. Yet, this script is writtern with intention to use for 
### mapped results to reference genome containing multiple chromosomes and potential
### multi-copy genes. Thus, the author has decided at the end to maintain the four fields
### comparison approach.

my @c_array;
my $c_array_count = 0;
my $read_count = 0;
while (<SAM>) {
	
	### To filter out the SAM file header for process, 
	### which always begin with an "@" sign. Those lines 
	### that do not have an "@' at the begining will be 
	### processed further.
	if ($_ =~ m/^@/) {
		print $_;
	} else {
		### Counter to report progress to the standard output (e.g. monitor)
		### Useful for those that are inpatient or want to keep close tab on the progress, such as myself.


		chomp $_;
		my @sam_line = split(/\t/, $_);
		my $line_entry = "$sam_line[2]\t$sam_line[3]\t$sam_line[5]\t$sam_line[9]";
		if (length ($sam_line[9]) > 149){  #&& length ($sam_line[9])){ ### modification per user request.
			$read_count++;
			if (($read_count % 10000) == 0){
				print "Processed Read Count: $read_count Unique read count: $c_array_count\n";
				print LOG "Processed Read Count: $read_count Unique read count: $c_array_count\n";
			}
			if ($c_array_count == 0) {
				$c_array[0][0] = $line_entry;
				$c_array[0][1] = 1;
				$c_array_count++;
				#print "$c_array[0][0]\t$c_array[0][1]\n";
			} else {
				my $match_flag = 0;
				for (my $i = 0; $i < $c_array_count; $i++){			
					if ($line_entry eq $c_array[$i][0]){
						$c_array[$i][1]++;
						#print "$c_array[$i][0]\t$c_array[$i][1]\n";
						$match_flag = 1;
						$i = $c_array_count + 1;
					}
				}
				if ($match_flag == 0){
					$c_array[$c_array_count][0] = $line_entry;
					$c_array[$c_array_count][1] = 1;
					#print "$c_array[$c_array_count][0]\t$c_array[$c_array_count][1]\n";
					#print "$c_array_count\n";
					$c_array_count++;
				}			
			}
		}
	}
}

print "Total number of $read_count raw reads collapsed to $c_array_count after primary consolidation.\n";
print LOG "Total number of $read_count raw reads collapsed to $c_array_count after primary consolidation.\n";

print "Done collapsing and consolating sequence from SAM files.\nStart removing adaptor sequence and locate gene end.\n";
print LOG "Done collapsing and consolating sequence from SAM files.\nStart removing adaptor sequence and locate gene end.\n";

my @t_array;
my $t_array_count = 0;
my $t_process_count = 0;
for (my $i = 0; $i < $c_array_count; $i++){
#for (my $i = 0; $i < 10; $i++){
	$t_process_count++;
	if (($t_process_count % 1000) == 0) {
		print "Trimming processd: $t_process_count reads\n";
		print LOG "Trimming processd: $t_process_count reads\n";
	}
	#print "line $i \n";
	#print LOG "$c_array[$i][0]\t$c_array[$i][1]\n";
	#print "$c_array[$i][0]\t$c_array[$i][1]\n";
	my $collapsed_count = $c_array[$i][1];
	my @split_line = split(/\t/, $c_array[$i][0]);
	my @output = split(/\t/, find_adaptor($split_line[3]));
	my $first_half = substr ($split_line[3], 0, $output[1]);
	my $adaptor = substr ($split_line[3], $output[1]);
	my @output2 = split(/\t/, find_end($first_half));
	my $tail = substr ($first_half, ($output2[1] + length($TERC_seq))-25); ## 25 was 20.
	#$c_array[$i][2] = $tail; ### Not sure if it is still needed

	
	#if (($tail =~ m/^AGGA/) == 0 ){ print "$i\t$tail\n";} ### For troubleshooting, comment out when done.
	#print "$i\t$tail\n"; if (($tail =~ m/^AGGA/) == 0 ){ sleep 5; } ### For troubleshooting, comment out when done.
	#if (($tail =~ m/^AGG..T....T......TGC/) == 0 ){ print "$i\t$tail\n"; } ### For troubleshooting, comment out when done.
	if ($t_array_count == 0){
		$t_array[0][0] = $tail;
		$t_array[0][1] = $collapsed_count;
		$t_array[0][2] = $first_half;
		#print "Match:\t$tail\t$t_array[0][0] $t_array[0][1]\n"; ### For troubleshooting, comment out when done.
		$t_array_count++
	} else {
		my $match_flag = 0;
		for (my $j = 0; $j < $t_array_count; $j++){
			if ($t_array[$j][0] eq $tail) {
				$t_array[$j][1] += $collapsed_count;
				#print "Match:\t$tail\t$t_array[$j][0] $t_array[$j][1] (+$collapsed_count)\n"; ### For troubleshooting, comment out when done.
				$match_flag = 1;
				$j += $t_array_count;
			}
		}
		
		if ($match_flag == 0){
			$t_array[$t_array_count][0] = $tail;
			$t_array[$t_array_count][1] = $collapsed_count;
			$t_array[$t_array_count][2] = $first_half;
			#print "New:\t$tail\t$t_array[$t_array_count][0] $t_array[$t_array_count][1]\n"; ### For troubleshooting, comment out when done.
			$t_array_count++;
		}
	}
}

print "Done removing adaptor sequence and identify the end of gene.\n";
print LOG "Done removing adaptor sequence and identify the end of gene.\n";
print "Start splitting up polyA tail and cleaning up data...\n";
print LOG "Start splitting up polyA tail and cleaning up data...\n";

my $t_array_read_count_sum = 0;
my @m_array;
my $m_array_count = 0;
for (my $i = 0 ; $i < $t_array_count; $i++){
	$t_array_read_count_sum += $t_array[$i][1];
	my @split_result = split(/\t/, split_tail($t_array[$i][0]));
	my $genomic_match = substr ($TERCjunc, 0, (length($t_array[$i][0])+5));
	print "Line# $i\n";
	print "SAM_seq $t_array[$i][2]\n";
	print "Genome_seq:\t$genomic_match\nRead_seq:\t$t_array[$i][0]\tCount: $t_array[$i][1]\n";
	print "Match_pattern:\t$split_result[0]\nRaw_split:\t$split_result[1]  $split_result[2]\nCorr_split:\t$split_result[3]  $split_result[4]\n";
	print "\n";

	print LOG "Line# $i\n";
	print LOG "SAM_seq $t_array[$i][2]\n";
	print LOG "Genome_seq:\t$genomic_match\nRead_seq:\t$t_array[$i][0]\tCount: $t_array[$i][1]\n";
	print LOG "Match_pattern:\t$split_result[0]\nRaw_split:\t$split_result[1]  $split_result[2]\nCorr_split:\t$split_result[3]  $split_result[4]\n";
	print LOG "\n";


	### Remove clean up tags
	$split_result[4] =~ s/\#//g;
	$split_result[4] =~ s/\*//g;
	my $tail_flag = 0;
	if ($split_result[4] eq "_") {
		$tail_flag = 0;
	} else {
		$tail_flag = 1;
	}

	
	if ($m_array_count == 0){
		$m_array[0][0] = $split_result[3];
		$m_array[0][1] = $split_result[4];
		$m_array[0][2] = $t_array[$i][1];
		$m_array[0][3] = $tail_flag;
		$m_array[0][4] = length($split_result[3]);
		$m_array[0][5] = length($split_result[4]);
		$m_array[0][6] = "$i ";
		$m_array_count++;
	} else {
		my $m_match_flag = 0;
		for (my $k = 0; $k < $m_array_count; $k++){
			if ($split_result[3] eq $m_array[$k][0] && $split_result[4] eq $m_array[$k][1]){
				$m_array[$k][2] += $t_array[$i][1];
				$m_array[$k][6] .= "$i ";
				$m_match_flag = 1;
				$k += $m_array_count;
			}
		}

		if ($m_match_flag == 0){
			$m_array[$m_array_count][0] = $split_result[3];
			$m_array[$m_array_count][1] = $split_result[4];
			$m_array[$m_array_count][2] = $t_array[$i][1];
			$m_array[$m_array_count][3] = $tail_flag;
			$m_array[$m_array_count][4] = length($split_result[3]);
			$m_array[$m_array_count][5] = length($split_result[4]);
			$m_array[$m_array_count][6] .= "$i ";
			$m_array_count++;
		}
	}
}
print "Checksum: Total number of $t_array_read_count_sum reads.\n";
print "Number of unique read after secondary collapse/consolidation: $t_array_count\n";
print LOG "Checksum: Total number of $t_array_read_count_sum reads.\n";
print LOG "Number of unique read after secondary collapse/consolidation: $t_array_count\n\n";
print OUT "Input read count from file $ARGV[0] : $t_array_read_count_sum\n";

my @sorted_m_array = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] || $a->[5] <=> $b->[5] } @m_array;

my $m_array_read_count = 0;
for (my $i = 0; $i < $m_array_count ; $i++){
	$m_array_read_count += $sorted_m_array[$i][2];
	print "$sorted_m_array[$i][2]\t$sorted_m_array[$i][0]\t$sorted_m_array[$i][1]\t$sorted_m_array[$i][4] $sorted_m_array[$i][5]\n";
	print LOG "$sorted_m_array[$i][2]\t$sorted_m_array[$i][0]\t$sorted_m_array[$i][1]\t$sorted_m_array[$i][4] $sorted_m_array[$i][5]\n";
}

print "Chechsum: $m_array_read_count\n";

my $unfilter_results = filter_reads(\@sorted_m_array,$m_array_count,"");
$unfilter_results = final_collapse($unfilter_results,"");

my $CG_filter_results = filter_reads(\@sorted_m_array,$m_array_count,"GC");
$CG_filter_results = final_collapse($CG_filter_results,"CG");

my $CGT_filter_results = filter_reads(\@sorted_m_array,$m_array_count,"GCT");
$CGT_filter_results = final_collapse($CGT_filter_results,"GCT");

print "\n$unfilter_results";
print "\n$CG_filter_results";
print "\n$CGT_filter_results";

print LOG "\n$unfilter_results";
print LOG "\n$CG_filter_results";
print LOG "\n$CGT_filter_results";

print OUT "\n$unfilter_results";
print OUT "\n$CG_filter_results";
print OUT "\n$CGT_filter_results";

		

#print "$no_tail_filter_output\n$AT_tail_filter_output\n$A_tail_filter_output";

close LOG;
close OUT;
close SAM;
exit;

###########################
### Sub routine section ###
###########################

sub final_collapse {
	my $sub_results;
	my @line_array = split(/\|/, $_[0]);

	### Second arhument of the subroutine recording the filtered base(s)."
	### Formating the tail display according the filtered base from tail
	my $unfiltered_base = 'ACGT';
	if ($_[1] =~ m/C/) { $unfiltered_base =~ s/C//g; }
	if ($_[1] =~ m/G/) { $unfiltered_base =~ s/G//g; }	
	if ($_[1] =~ m/T/) { $unfiltered_base =~ s/T//g; }	
	my @final_tail = split (//, $unfiltered_base);
	my $tail_display = join (':', @final_tail);
	$tail_display = '_['.$tail_display.']n';

	### The foreach look that take the input list, collapse base on the genomics sequence
	### sum the read count, perform basic statisitc calculation.

	#print "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStdev\n";
	$sub_results .= "Base from tail filtered : $_[1]\n";
	$sub_results .= "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStdev\n";
	my @with_tail_array;
	my $with_tail_count = 0;
	my $sum_of_reads = 0;
	foreach (@line_array) {
		my @temp_line = split(/\t/, $_);
		#print "@temp_line\n";

		### [0] = sequence
		### [1] = total read count
		### [2] = min tail length
		### [3] = max tail length
		### [4] = sum tail length
		### [5] = mean tail length
		### [6] = sum of square
		### [7] = stdev of tail length


		if ($temp_line[1] eq '_') {
			#print "$temp_line[0]"."_\t$temp_line[2]\tN/A\tN/A\tN/A\tN/A\tN/A\n";
			$sub_results .= "$temp_line[0]"."_\t$temp_line[2]\tN/A\tN/A\tN/A\tN/A\tN/A\n";
			$sum_of_reads += $temp_line[2];
		} else {
			$sum_of_reads += $temp_line[2];
			if ($with_tail_count == 0){
				$with_tail_array[$with_tail_count][0] = $temp_line[0];
				$with_tail_array[$with_tail_count][1] = $temp_line[2];
				$with_tail_array[$with_tail_count][2] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][3] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][4] = length ($temp_line[1]) * $temp_line[2];
				$with_tail_array[$with_tail_count][5] = $with_tail_array[$with_tail_count][4] / $with_tail_array[$with_tail_count][1];
				$with_tail_array[$with_tail_count][6] = sq($with_tail_array[$with_tail_count][4]); #* $with_tail_array[$with_tail_count][4];
				$with_tail_count++;
			} elsif ( (length $temp_line[0]) == (length $with_tail_array[$with_tail_count-1][0]) ) {
				$with_tail_array[$with_tail_count-1][1] += $temp_line[2];				
				$with_tail_array[$with_tail_count-1][3] = (length $temp_line[1]);
				$with_tail_array[$with_tail_count-1][4] += (length $temp_line[1]) * $temp_line[2];
				$with_tail_array[$with_tail_count-1][5] = $with_tail_array[$with_tail_count-1][4] / $with_tail_array[$with_tail_count-1][1];
				$with_tail_array[$with_tail_count-1][6] += sq($with_tail_array[$with_tail_count-1][4]); #* $with_tail_array[$with_tail_count-1][4];

			} else {
				$with_tail_array[$with_tail_count][0] = $temp_line[0];
				$with_tail_array[$with_tail_count][1] = $temp_line[2];
				$with_tail_array[$with_tail_count][2] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][3] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][4] = length ($temp_line[1]) * $temp_line[2];	
				$with_tail_array[$with_tail_count][5] = $with_tail_array[$with_tail_count][4] / $with_tail_array[$with_tail_count][1];		
				$with_tail_array[$with_tail_count][6] = sq($with_tail_array[$with_tail_count][4]); #* $with_tail_array[$with_tail_count][4];
				$with_tail_count++;
			}
		}
	}

	### Reformatting the results as final output

	
	for (my $i = 0 ; $i < $with_tail_count; $i++) {
		my $stdev;
		if ($with_tail_array[$i][1] == 1){
			$stdev = "N/A";
		} else {
			$stdev = sqrt(($with_tail_array[$i][6] - sq($with_tail_array[$i][5]*$with_tail_array[$i][1]))/($with_tail_array[$i][1]-1));	
			$stdev = nearest(0.0001, $stdev);
		}
		#print "$with_tail_array[$i][0]$tail_display\t$with_tail_array[$i][1]\t$with_tail_array[$i][2]\t$with_tail_array[$i][3]\t".
		#"$with_tail_array[$i][4]\t$with_tail_array[$i][5]\t$stdev\n";
		$with_tail_array[$i][5] = nearest(0.001, $with_tail_array[$i][5]);
		$sub_results .= "$with_tail_array[$i][0]$tail_display\t$with_tail_array[$i][1]\t$with_tail_array[$i][2]\t$with_tail_array[$i][3]\t".
		"$with_tail_array[$i][4]\t$with_tail_array[$i][5]\t$stdev\n";
	}
	#print "Total read count: $sum_of_reads\n";
	$sub_results .= "Total read count: $sum_of_reads\n";
}

sub sq {
	my $square = $_[0]*$_[0];
	return $square;
}

sub filter_reads {
	my @subarray = @{$_[0]};
	my $subarray_count = $_[1];
	my @sub_filter = split (//, $_[2]);
	my $filter_count = length($_[2]);

	my @output_array;
	my $output_array_count = 0;
	my $filter_flag;
	my $filtered_read_count = 0;
	my $discarded_read_count = 0;
	my $results_text = '';
	for (my $i = 0 ; $i < $subarray_count ; $i++){
		$filter_flag = 0;
		for (my $j = 0; $j < $filter_count; $j++){
			$filter_flag += ($subarray[$i][1] =~ m/$sub_filter[$j]/);
		}

		if (length($subarray[$i][0]) < 5) {
			$filter_flag = 1;
		}
	
		if ($filter_flag == 0){
			for (my $k = 0 ; $k < 3; $k++){
				$output_array[$output_array_count][$k] = $subarray[$i][$k];
				#print "$output_array[$output_array_count][$k]\t";
				$results_text .= "$output_array[$output_array_count][$k]\t";
			}
			$filtered_read_count += $subarray[$i][2];
			#print "\n";
			$results_text .= '|';
			$output_array_count++;
		} else {
			$discarded_read_count += $subarray[$i][2];
		}
	}
	print "Filtered read count: $filtered_read_count\n";
	print "Discarded read count: $discarded_read_count\n";

	return $results_text;
	
}


sub split_tail {
	my $input_length = length ($_[0]);
	my @tail = split(//, $_[0]);
	my @junc = split (//, substr ($TERCjunc, 0, $input_length));
	#print "$_[0] $input_length\n";
	
	my $match_pattern='';
	for (my $i = 0; $i < $input_length ; $i++) {
		#my $a = shift @tail;
		#my $b = shift @junc;
		my $a = $tail[$i];
		my $b = $junc[$i];
		my $next_a; my $prev_a;
		unless($i == 0){$prev_a = $tail[$i-1];} else {$prev_a = "N";}
		unless($i == ($input_length-1)){$next_a = $tail[$i+1];} else {$next_a = "N";}

		if ($a eq $b ) {			
			$match_pattern .= "1";
		} else {
			$match_pattern .= "0";
		}
		
	}

	my $best_split_coordinate = -1;
	my $best_split_score = -1;
	$match_pattern .= '0';
	for (my $i = 0; $i < $input_length ; $i++) {
		my $front_count = 0;
		my $end_count = 0;
	
		my $front = substr ($match_pattern, 0, $input_length - $i); my $front_length = length ($front);
		my $end = substr ($match_pattern, $input_length - $i, $i + 1); my $end_length = length ($end);
		while ($front =~ /1/g) {$front_count++}; $front_count = $front_count / ($front_length+$end_length);
		while ($end =~ /0/g) {$end_count++}; $end_count = $end_count / ($end_length+$front_length);
		my $delta = $front_count + $end_count;

		if ($delta == 1) {
			$best_split_coordinate = $input_length - $i;
			$best_split_score = 1;
			$i = $input_length + 1;
		} elsif ($delta >= $best_split_score ) {
			$best_split_coordinate = $input_length - $i;
			$best_split_score = $delta;
		} 
		my $rounded_delta = nearest(0.01,$delta);
		my $rounded_best_split_score = nearest(0.01,$best_split_score);
		#print "$front"."|$end\t$rounded_delta\t$rounded_best_split_score\t$best_split_coordinate\n";

	}


	my $input_front_split = substr ($_[0], 0, $best_split_coordinate);
	my $input_end_split = substr ($_[0], $best_split_coordinate);
	if (length($input_end_split) == 0){
			$input_end_split = '_';
	}
	#my $input_front_split_corrected = substr ($TERCjunc, 0, $best_split_coordinate); 

	### Perform some clean up steps
	### If the tail is less 5 bases long or less, and does not contain any A, merge with the genomics sequence
	#if ($input_end_split ne "_" && length($input_end_split) < 6 && ($input_end_split =~ /A/) == 0 ){
	#	$input_front_split .= $input_end_split;
	#	$input_end_split = "_";
	#}

	### Run additional clean up step on those reads with tails.
	### Experimental Step!!
	if ($input_end_split ne "_"){	

		### Clean up 1, remove residues adaptor sequence that was fail to remove during the first pass.
		### Recycling the find_adaptor subroutine, with a cut off value of 0.33 or above as good match.
		### The tail is further turnacted if the adaptor match coordinate is not 0.
		### An "*" is added to the tail to mark secondary processing.
		my @adaptor_test = split(/\t/, find_adaptor($input_end_split));
		#print "$input_end_split\t@adaptor_test\n";
		if ($adaptor_test[0] >= 0.33 && $adaptor_test[1] == 0) {
			$input_end_split = '_*';
		} elsif ($adaptor_test[0] >= 0.33 && $adaptor_test[1] != 0){
			$input_end_split = substr($input_end_split,0,$adaptor_test[1]);
			$input_end_split .= '*';
		}

		### Clean up 2, merge short tail (2 nucleotides or less) that does not contain any "A" with the genomic sequence.
		
		### Disabled per user request
		#if ( $input_end_split ne "_" && $input_end_split ne "_*"  && length($input_end_split) < 4 && ($input_end_split =~ /A/) == 0 && ($input_end_split =~ /T/) == 0 ){
		#	$input_front_split .= $input_end_split;
		#	$input_end_split = "_#";
		#}
	}

	my $input_front_split_corrected = substr ($TERCjunc, 0, length($input_front_split)); 	
	$match_pattern .= "\t$input_front_split\t$input_end_split\t$input_front_split_corrected\t$input_end_split";
	#print "$match_pattern\n";
	return $match_pattern;

}



sub find_end{
	my $seq = $_[0];
	my $hi_score = 0;
	my $hi_score_coordinate = -1;
	my $test_length = length ($seq);

	#for (my $i = (length ($seq) - 0) ; $i > (length ($seq) - length ($ad_seq) - 10); $i--) {
	for (my $i = length ($seq) ; $i > 0; $i--) {
		my $score = comp_end_score(substr ($seq, $i));
		#print "$score\n";
		if ($score > $hi_score){
			$hi_score = $score;
			$hi_score_coordinate = $i;
		}
	}

	my $results = "$hi_score\t$hi_score_coordinate";
	#print "$results\n";
	#exit;

	return $results;
}


sub comp_end_score {

	my $seq1 = $_[0];
	my $test_length = length($TERC_seq);
	my $score = 0;	

	if (length($seq1) < length($TERC_seq)){
		$test_length = length($seq1);
	}

	if ($seq1 =~ m/$TERC_seq$/) {
		$score = 45;
	} else {

		my @array_a = split(//, $seq1);
		my @array_b = split(//, $TERC_seq);

		for (my $i = 0; $i < $test_length; $i++){
			my $a = shift @array_a;
			my $b = shift @array_b;
			if ($a eq $b) {
				if ($i < $test_length/4){  ### lowest 25%
					$score += 0.25;
				} elsif ($i >= $test_length/4 && $i < $test_length/2 ) {
					$score += 0.5
				} elsif ($i >= $test_length/2 && $i < ($test_length/4)*3) {
					$score += 1.0;
				} elsif ($i >= ($test_length/4)*3) {
					$score += 2.0;
				} 
			}
		}
	}
	#print "comp_score_input: $_[0]\t$test_length\t$score\n";

	my $return = nearest(0.01, $score / length($TERC_seq));
	return ($return);
}

sub find_adaptor{
	my $seq = $_[0];
	my $hi_score = 0;
	my $hi_score_coordinate = -1;
	my $test_length = nearest(1 , length($ad_seq)/2);

	#for (my $i = (length ($seq) - 0) ; $i > (length ($seq) - length ($ad_seq) - 10); $i--) {
	for (my $i = length ($seq) ; $i > (length ($seq) - length ($ad_seq)- 10); $i--) {
		my $score = comp_adp_score(substr ($seq, $i));
		#print "$score\n";
		if ($score > $hi_score){
			$hi_score = $score;
			$hi_score_coordinate = $i;
		}
	}

	my $results = "$hi_score\t$hi_score_coordinate";
	#print "$results\n";
	#exit;

	return $results;
}

sub comp_adp_score {

	my $seq1 = $_[0];
	my $test_length = length($ad_seq);
	my $score = 0;	

	if (length($seq1) < length($ad_seq)){
		$test_length = length($seq1);
	}

	if ($seq1 eq $ad_seq) {
		$score = 20.25;
	} else {

		my @array_a = split(//, $seq1);
		my @array_b = split(//, $ad_seq);

		for (my $i = 0; $i < $test_length; $i++){
			my $a = shift @array_a;
			my $b = shift @array_b;
			if ($a eq $b) {
				if ($i == 0){
					$score += 1.0;
				} elsif ($i <= 3) {
					$score += 2.0
				} elsif ($i <=6) {
					$score += 0.5;
				} elsif ($i >= $test_length/4) {
					$score += 0.25;
				} else {
					$score += 0.5;
				}
			}
		}
	}
	#print "comp_score_input: $_[0]\t$test_length\t$score\n";

	my $return = nearest(0.01, $score / (length($ad_seq)* 0.75));
	return ($return);
}



