#!/usr/bin/perl -w
### Writtern by Albert K Tai, PhD
### Research Assistant Professor, Department of Integrative Phyisology and Pathobiology
### Tufts University School of Medicine, Boston, MA

use strict;
use Math::Round;

### Sequence information 

#my $ad_seq = "CTGTAGGCACCATCAATCGTTACGTAG";
my $ad_seq = "CTGTAGGCACCATCAATCGTTAC";
#my $TERC_seq = "AGGCCGCAGGAAGAGGAACGGAGCGAGTCCCCGCGCGCGGCGCGAT"."TCCCTGAGCTGTGGGACGTGCACCCAGGACTCGGCTCACACATGC"; ## TERC
# my $Target_seq = "TCCCTGAGCTGTGGGACGTGCACCCAGGACTCGGCTCACACATGC";
my $Target_seq =  "AGCAGATGGAGTGTTTTCTGTTGAACACTAAGTACTGCCACAAGT"; # scaRNA13 2021.01.15
#my $TERCjunc = "TCCCTGAGCTGTGGGACGTGCACCC AGGACTCGGCTCACACATGCAGTTCGCTTTCCTGTTGGTGGGGGGAACGCCGATCGTGCGCATCCGTCACCCCTCGCCGGCA
#ATGGGGGCTTGTGAACCCCCAAACCTGACTGACTGGGCCAGTGTGCTGCAAATTGGCAGGAGACGTGAAGGCACCTCCAAAGTCGGCCAAAA";
#my $Target_junc = "CACCCAGGACTCGGCTCACACATGCAGTTCGCTTTCCTGTTGGTGGGGGGAACGCCGATCGTGCGCATCCGTCACCCCTCGCCGGCA
#ATGGGGGCTTGTGAACCCCCAAACCTGACTGACTGGGCCAGTGTGCTGCAAATTGGCAGGAGACGTGAAGGCACCTCCAAAGTCGGCCAAAA";
my $Target_junc = "TTGAACACTAAGTACTGCCACAAGTTACTTTTTTTTTTTTAAACTTTGAGTATTTTTTTACAATGTTGCTGGAGGTGATCTGTTTAT
GCTTTGAGAGTGTTCGAATTTAAAATCAGAAAATCATGTCAGTGAGTGAGTCTTTCAAATAATCCTTCGGCATGAAACCTGAGCCTAGTAAC"; # scaRNA13 2021.01.15

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

### @c_array_count store the line entry, which contains the reference, mapped position, CIGAR and Sequnece, stored at first column [0]
### It then keep the count of the exact same line entry in the second column [1]

my $echo;
my @c_array;
my $c_array_count = 0;
my $read_count = 0;
my $collapsed_count = 0;
my $too_short_count = 0;
my $raw_read_length_cutoff = 149;
while (<SAM>) {
	
	### To filter out the SAM file header for process, 
	### which always begin with an "@" sign. Those lines 
	### that do not have an "@' at the begining will be 
	### processed further.
	if ($_ =~ m/^@/) {
		print $_;
	} else {
		chomp $_;
		my @sam_line = split(/\t/, $_);
		### [2] reference
		### [3] position
		### [5] CIGAR
 		### [9] Sequence
		my $line_entry = "$sam_line[2]\t$sam_line[3]\t$sam_line[5]\t$sam_line[9]";
		if (length ($sam_line[9]) > $raw_read_length_cutoff){  #&& length ($sam_line[9])){ ### modification per user request.
			$read_count++;

			### Counter to report progress to the standard output (e.g. monitor)
			### Useful for those that are inpatient or want to keep close tab on the progress, such as myself.
			if (($read_count % 10000) == 0){
				$echo = "Processed Read Count: $read_count Unique read count: $c_array_count\n";
				print $echo;
				print LOG $echo; #"Processed Read Count: $read_count Unique read count: $c_array_count\n";
			} 

			
			### @c_array_count 

			if ($c_array_count == 0) {
				$c_array[0][0] = $line_entry;
				$c_array[0][1] = 1;
				$c_array_count++;
				#print "$c_array[0][0]\t$c_array[0][1]\n"; ###
				$collapsed_count++;
			} else {
				my $match_flag = 0;
				for (my $i = 0; $i < $c_array_count; $i++){			
					if ($line_entry eq $c_array[$i][0]){
						$c_array[$i][1]++;  ## if a match to line entry is found, add 1 to the line entry count
						#print "$c_array[$i][0]\t$c_array[$i][1]\n"; ###
						$match_flag = 1;
						$i = $c_array_count + 1; ## exit the loop by increasing i when a match is found.
						$collapsed_count++;
					}
				}
				### if no match is found, add a new line entry and set the count to 1
				if ($match_flag == 0){
					$c_array[$c_array_count][0] = $line_entry;
					$c_array[$c_array_count][1] = 1;
					#print "$c_array[$c_array_count][0]\t$c_array[$c_array_count][1]\n"; ###
					#print "$c_array_count\n"; ##3
					$c_array_count++;
					$collapsed_count++;
				}			
			} 		
		} else {
			$too_short_count++;
		}
	}
}

### Display the entries with most and least count. Number to display based on $list_count.
#my $list_count=20;
#my @temp_array = sort { $b->[1] <=> $a->[1] } @c_array;
#for (my $i = 0; $i < $list_count; $i++){
#	print "Count:$temp_array[$i][1]\t$temp_array[$i][0]\n";
#}
#@temp_array = sort { $a->[1] <=> $b->[1] } @c_array;
#for (my $i = 0; $i < $list_count; $i++){
#	print "Count:$temp_array[$i][1]\t$temp_array[$i][0]\n";
#}

my $total_raw_count = $read_count + $too_short_count;
$echo = "Total number of $total_raw_count raw reads, with $read_count collapsed to $c_array_count after primary consolidation and $too_short_count reads short than $raw_read_length_cutoff bp discarded.\n";
print "$echo";
print LOG $echo;

$echo = "Done collapsing and consolating sequence from SAM files.\nStart removing adaptor sequence and locate gene end.\n";
print $echo;
print LOG $echo;

### @t_array is table for containing adptor trimmed reads
### [0]
### [1] count
### [2]


my @t_array;
my $t_array_count = 0;
my $t_process_count = 0;
for (my $i = 0; $i < $c_array_count; $i++){
#for (my $i = 0; $i < 10; $i++){
	$t_process_count++;
	### Progress counter
	if (($t_process_count % 10000) == 0) {
		$echo = "Trimming processed: $t_process_count reads\n";		
		print $echo;
		print LOG $echo; #"Trimming processed: $t_process_count reads\n";
		#print "line $i\t";
		#print "$c_array[$i][0]\t$c_array[$i][1]\n";
	}

	#print "line $i \n";
	#print LOG "$c_array[$i][0]\t$c_array[$i][1]\n";
	#print "$c_array[$i][0]\t$c_array[$i][1]\n";
	my $collapsed_count = $c_array[$i][1];
	my @split_line = split(/\t/, $c_array[$i][0]);
	### find_adaptor subroutine
	my @output = split(/\t/, find_adaptor($split_line[3]));
	my $first_half = substr ($split_line[3], 0, $output[1]);
	my $adaptor = substr ($split_line[3], $output[1]);
	### find end subroutine
	my @output2 = split(/\t/, find_end($first_half));
	my $tail = substr ($first_half, ($output2[1] + length($Target_seq))-25); ## 25 was 20.
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

### troubleshooting
#my @sorted_t_array = sort { $b->[1] <=> $a->[1] } @t_array;
#for (my $i = 0; $i < $t_array_count; $i++){
	#print "$sorted_t_array[$i][0]\t$sorted_t_array[$i][1]\t$sorted_t_array[$i][2]\n";
#}

$echo = "Done removing adaptor sequence and identify the end of gene.\n".
        "Start splitting up polyA tail and cleaning up data...\n";
print $echo;
print LOG $echo; #"Done removing adaptor sequence and identify the end of gene.\n";
#print LOG "Start splitting up polyA tail and cleaning up data...\n";

### @m_array = merged_array, contain reads from secondary collapses, which correct for SNP
### [1] 

my $t_array_read_count_sum = 0;
my @m_array;
my $m_array_count = 0;
#my $a_count = 0;
#my $b_count = 0;
#my $c_count = 0;
my $low_score_cutoff = 0.5;
my $low_score_count = 0;
my $low_score_record = '';
for (my $i = 0 ; $i < $t_array_count; $i++){
	$t_array_read_count_sum += $t_array[$i][1];

	## split_tail subroutine
	my @split_result = split(/\t/, split_tail($t_array[$i][0]));  
	my $genomic_match = substr ($Target_junc, 0, (length($t_array[$i][0])+5));
	
	$echo = "Line# $i\n".
	        "SAM_seq $t_array[$i][2]\n".
	        "   Genome_seq:\t$genomic_match\n     Read_seq:\t$t_array[$i][0]\tCount: $t_array[$i][1], Score: $split_result[5]\n".
	        "Match_pattern:\t$split_result[0]\n    Raw_split:\t$split_result[1]  $split_result[2]\n   Corr_split:\t$split_result[3]  $split_result[4]\n".
	        "\n";
	print $echo;
	print LOG $echo;

	### split_results[3] = template-corrected, splited 5'-template 
	### split_results[4] = splited 3'-tail 
	### Remove clean-up tags
	$split_result[4] =~ s/\#//g;
	$split_result[4] =~ s/\*//g;
	my $tail_flag = 0;
	if ($split_result[4] eq "_") {
		$tail_flag = 0;
	} else {
		$tail_flag = 1;
	}

	if ($split_result[5] >= $low_score_cutoff){
		if ($m_array_count == 0){
			$m_array[0][0] = $split_result[3];			# template-corrected, splited 5'-template
			$m_array[0][1] = $split_result[4];			#
			$m_array[0][2] = $t_array[$i][1];			# count
			$m_array[0][3] = $tail_flag;
			$m_array[0][4] = length($split_result[3]);  #template-corrected, splited 5'-template
			$m_array[0][5] = length($split_result[4]);  #
			$m_array[0][6] = "$i ";
			$m_array_count++;
			#$a_count = $t_array[$i][1];
		} else {
			my $m_match_flag = 0;
			for (my $k = 0; $k < $m_array_count; $k++){
				if ($split_result[3] eq $m_array[$k][0] && $split_result[4] eq $m_array[$k][1]){
					$m_array[$k][2] += $t_array[$i][1];
					$m_array[$k][6] .= "$i ";
					$m_match_flag = 1;
					$k += $m_array_count;
					#$b_count = $b_count + $t_array[$i][1];
				}
			}

			if ($m_match_flag == 0){
				$m_array[$m_array_count][0] = $split_result[3];  	### 5'template
				$m_array[$m_array_count][1] = $split_result[4];  	### 3' tail, if exist
				$m_array[$m_array_count][2] = $t_array[$i][1];   	### collapsed count
				$m_array[$m_array_count][3] = $tail_flag;			### with or without tail
				$m_array[$m_array_count][4] = length($split_result[3]);	### length of template
				$m_array[$m_array_count][5] = length($split_result[4]);	### length of tail
				$m_array[$m_array_count][6] .= "$i ";				### pre-collapse line item, after primary collapse but before secondary sequence correction collapse.
				$m_array_count++;
				#$c_count = $c_count + $t_array[$i][1];
			}
		}
	} elsif ($split_result[5] < $low_score_cutoff){
		$low_score_count++;
		$low_score_record .= " $i($t_array[$i][1])"
	}

}

if ($low_score_record eq ""){$low_score_record = "None";}
$echo = "Checksum: Total number of $t_array_read_count_sum reads.\n".
		"Number of unique read after secondary collapse/consolidation: $m_array_count\n".
		"Number of read discarded because of low score $low_score_count, line#(count): $low_score_record\n";
print $echo;
print LOG $echo;
print OUT "Input read count from file $ARGV[0] : $t_array_read_count_sum\n";

my @sorted_m_array = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] || $a->[5] <=> $b->[5] } @m_array;

#for (my $i = 0 ; $i < $m_array_count ; $i++){
#	print "$i\t$sorted_m_array[$i][0]\t$sorted_m_array[$i][1]\t$sorted_m_array[$i][2]\t$sorted_m_array[$i][3]\t$sorted_m_array[$i][4]\t$sorted_m_array[$i][6]\n";
#}

### $sorted_m_array[2] = collapsed count
### $sorted_m_array[0] = 5' template
### $sorted_m_array[1] = 3' tail
### $sorted_m_array[4] = template length
### $sorted_m_array[5] = tail length
my $m_array_read_count = 0;
for (my $i = 0; $i < $m_array_count ; $i++){
	$m_array_read_count += $sorted_m_array[$i][2];
	$echo = "$sorted_m_array[$i][2]\t$sorted_m_array[$i][0]\t$sorted_m_array[$i][1]\t$sorted_m_array[$i][4] $sorted_m_array[$i][5]\n";
	print $echo;
	print LOG $echo;  #"$sorted_m_array[$i][2]\t$sorted_m_array[$i][0]\t$sorted_m_array[$i][1]\t$sorted_m_array[$i][4] $sorted_m_array[$i][5]\n";
}

## Good up to this point!!!
print "Checksum: $m_array_read_count\n";

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

close LOG;
close OUT;
close SAM;

exit;  


###########################
### Sub routine section ###
###########################

### Subroutine to find the RACE adaptor, fist appear in line ~147
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

### line ~150
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

sub split_tail {
	my $input_length = length ($_[0]);
	my @tail = split(//, $_[0]);
	my @junc = split (//, substr ($Target_junc, 0, $input_length));
	#print "$_[0] $input_length\n";
	
	### Create 0 and 1 match pattern for the two input sequence.
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

	
	## The match starts from 3' end, moving 5' one base at a time. The proportion of matches in the front section and mismatch in the end section is calculated
	## The sum of the two proportions will approach 1.000 as it approach the "best" splitting site.  
	my $best_split_coordinate = -1;
	my $best_split_score = -1;
	$match_pattern .= '0';
	for (my $i = 0; $i < $input_length ; $i++) {
		my $front_count = 0;
		my $end_count = 0;
		my $front = substr ($match_pattern, 0, $input_length - $i); my $front_length = length ($front);
		my $end = substr ($match_pattern, $input_length - $i, $i + 1); my $end_length = length ($end);
		while ($front =~ /1/g) {$front_count++}; while ($front =~ /0/g) {$front_count--};
		$front_count = $front_count / ($front_length+$end_length);  ## number of matches on the 5' end as divided by total length. 
		while ($end =~ /0/g) {$end_count++}; while ($end =~ /1/g) {$end_count--};
		$end_count = $end_count / ($end_length+$front_length);
		my $delta = $front_count + $end_count;

		if ($delta == 1) {
			$best_split_coordinate = $input_length - $i;
			$best_split_score = 1;
			$i = $input_length + 1;
		} elsif ($delta > $best_split_score ) {    ### used to be >= instead of >
		### >= prioritize longer tail and > prioritize longer template, when ambiguous
			$best_split_coordinate = $input_length - $i;
			$best_split_score = $delta;
		} 
		my $rounded_delta = nearest(0.001,$delta);
		my $rounded_best_split_score = nearest(0.001,$best_split_score);
		#print "$front"."|$end"."\t$front_count\t$end_count\t$rounded_delta\t$rounded_best_split_score\t$best_split_coordinate\n";  ###
	}


	my $input_front_split = substr ($_[0], 0, $best_split_coordinate);
	my $input_end_split = substr ($_[0], $best_split_coordinate);
	if (length($input_end_split) == 0){
			$input_end_split = '_';
	}

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

	my $input_front_split_corrected = substr ($Target_junc, 0, length($input_front_split)); 	
	$match_pattern .= "\t$input_front_split\t$input_end_split\t$input_front_split_corrected\t$input_end_split\t$best_split_score";
	#print "$match_pattern\n";
	return $match_pattern;

}

sub final_collapse {
	my $sub_results;
	my @line_array = split(/\|/, $_[0]);  ##???????????
	#print "$_[0]\n"; ###
	### Second argument of the subroutine recording the filtered base(s)."
	### Formating the tail display according the filtered base from tail
	my $unfiltered_base = 'ACGT';
	if ($_[1] =~ m/C/) { $unfiltered_base =~ s/C//g; }  ## if the argument contains a C, the C will be removed from $unfiltered_base and becomes "AGT"
	if ($_[1] =~ m/G/) { $unfiltered_base =~ s/G//g; }	
	if ($_[1] =~ m/T/) { $unfiltered_base =~ s/T//g; }	

	## To create display "[A:C:G:T]n" and adds to the sequence with non-templated tails. Visual only
	my @final_tail = split (//, $unfiltered_base);  
	my $tail_display = join (':', @final_tail);
	$tail_display = '_['.$tail_display.']n';

	### The foreach look that take the input list, collapse base on the genomics sequence
	### sum the read count, perform basic statisitc calculation.

	#print "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStdev\n";
	$sub_results .= "Base from tail filtered : $_[1]\n";
	$sub_results .= "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStd_dev\n";
	my @with_tail_array;
	my $with_tail_count = 0;
	my $sum_of_reads = 0;
	my $max_tail_length = 1;  ## added
	foreach (@line_array) {
		my @temp_line = split(/\t/, $_);
		#print "@temp_line\n";
		### $temp_line[0] = sequence
		### $temp_line[1] = "_" or non-templated 3' tail sequence
		### $temp_line[2] = collapsed count		

		### $with_tail_array[0] = sequence
		### $with_tail_array[1] = total read count
		### $with_tail_array[2] = min tail length
		### $with_tail_array[3] = max tail length
		### $with_tail_array[4] = sum tail length
		### $with_tail_array[5] = mean tail length
		### $with_tail_array[6] = sum of square ??
		### $with_tail_array[7] = stdev of tail length

		if ($temp_line[1] eq '_') {
			## If the sequece ended with "_" , i.e. no non-template tail
			#print "$temp_line[0]"."_\t$temp_line[2]\tN/A\tN/A\tN/A\tN/A\tN/A\n";
			$sub_results .= "$temp_line[0]"."_\t$temp_line[2]\tN/A\tN/A\tN/A\tN/A\tN/A\n";
			$sum_of_reads += $temp_line[2];
		} else {
			## if the sequence does NOT ended with "_" , i.e. with non-template 3' tail
			$sum_of_reads += $temp_line[2];
			if ($with_tail_count == 0){
				### Initialized the first record of the array
				$with_tail_array[$with_tail_count][0] = $temp_line[0]; print "0 $temp_line[0]\t$with_tail_array[$with_tail_count][0]\n";# seq
				$with_tail_array[$with_tail_count][1] = $temp_line[2]; # count 
				$with_tail_array[$with_tail_count][2] = length ($temp_line[1]); # min length
				$with_tail_array[$with_tail_count][3] = length ($temp_line[1]); # max length, will be updated in subsequence record
				$with_tail_array[$with_tail_count][4] = length ($temp_line[1]) * $temp_line[2]; # count * tail length
				$with_tail_array[$with_tail_count][5] = $with_tail_array[$with_tail_count][4] / $with_tail_array[$with_tail_count][1];
				for (my $i = 0 ; $i < $temp_line[2] ; $i++){ $with_tail_array[$with_tail_count][6] .= "$with_tail_array[$with_tail_count][3],"; }
				if ($max_tail_length < length ($temp_line[1]) ) {$max_tail_length = length ($temp_line[1]);}
				$with_tail_count++;
			} elsif ( (length $temp_line[0]) == (length $with_tail_array[$with_tail_count-1][0]) ) {
				## if the template has the same length as the previous entry, then assume it is the same template.
				$with_tail_array[$with_tail_count-1][1] += $temp_line[2];  ## add to the collapsed count				
				$with_tail_array[$with_tail_count-1][3] = (length $temp_line[1]);  ## updated max tail length
				$with_tail_array[$with_tail_count-1][4] += (length $temp_line[1]) * $temp_line[2];
				$with_tail_array[$with_tail_count-1][5] = $with_tail_array[$with_tail_count-1][4] / $with_tail_array[$with_tail_count-1][1];
				#$with_tail_array[$with_tail_count-1][6] += sq($with_tail_array[$with_tail_count-1][4]); #* $with_tail_array[$with_tail_count-1][4];
				#if ($max_tail_length < $with_tail_array[$with_tail_count][3]) {$max_tail_length = $with_tail_array[$with_tail_count][3];}
				if ($max_tail_length < length ($temp_line[1]) ) {$max_tail_length = length ($temp_line[1]);}
				for (my $i = 0 ; $i < $temp_line[2] ; $i++){ $with_tail_array[$with_tail_count-1][6] .= "$with_tail_array[$with_tail_count-1][3],"; }
			} else {
				$with_tail_array[$with_tail_count][0] = $temp_line[0];
				$with_tail_array[$with_tail_count][1] = $temp_line[2];
				$with_tail_array[$with_tail_count][2] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][3] = length ($temp_line[1]);
				$with_tail_array[$with_tail_count][4] = length ($temp_line[1]) * $temp_line[2];	
				$with_tail_array[$with_tail_count][5] = $with_tail_array[$with_tail_count][4] / $with_tail_array[$with_tail_count][1];		
				#$with_tail_array[$with_tail_count][6] = sq($with_tail_array[$with_tail_count][4]); #* $with_tail_array[$with_tail_count][4];
				for (my $i = 0 ; $i < $temp_line[2] ; $i++){ $with_tail_array[$with_tail_count][6] .= "$with_tail_array[$with_tail_count][3],"; }
				#$with_tail_array[$with_tail_count-1][6] =~ s/,$//;
				#if ($max_tail_length < $with_tail_array[$with_tail_count][3]) {$max_tail_length = $with_tail_array[$with_tail_count][3];}
				if ($max_tail_length < length ($temp_line[1]) ) {$max_tail_length = length ($temp_line[1]);}
				$with_tail_count++;
			}
		
		#$with_tail_array[$with_tail_count-1][6] =~ s/,$//;
		}
	}

	#my @temp = @with_tail_array;
	#for (my $j = 0 ; $j < $with_tail_count; $j++){
	#	print "$temp[$j][0]\t";
	#	my $stat;
	#	for (my $k = 1 ; $k < 7 ; $k++ ){
	#		print "$k:$temp[$j][$k]\t";
	#		if ($k == 6){ $stat = mean($temp[$j][$k]); }
	#	}
	#	print "$stat\tend\n";
	#}

	### Reformatting the results as final output
		### $with_tail_array[0] = sequence
		### $with_tail_array[1] = total read count
		### $with_tail_array[2] = min tail length
		### $with_tail_array[3] = max tail length
		### $with_tail_array[4] = sum tail length
		### $with_tail_array[5] = mean tail length
		### $with_tail_array[6] = sum of square ??
		### $with_tail_array[7] = stdev of tail length

	#print "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStdev\n";
	$sub_results .= "Base from tail filtered : $_[1]\n";
	$sub_results .= "Sequence\tCount\tMin_len\tMax_len\tSum_len\tMean_len\tStd_dev\tTail_len_Dist\tTail_len_ind";
	for (my $i = 1 ; $i < $max_tail_length + 1; $i++){ $sub_results .= "\t$i"; } 
	$sub_results .= "\n";

	for (my $i = 0 ; $i < $with_tail_count; $i++) {
		$with_tail_array[$i][6] =~ s/,$//;
		my @stat = split(/,/, statcal($with_tail_array[$i][6]));
		my $stdev = $stat[1];
		my $mean  = $stat[0];
		my $sum   = $stat[2];
		#if ($with_tail_array[$i][1] == 1){
			#$stdev = "N/A";
		#} else {
			#$mean = $stat[0]	
			#$stdev = $stat[1];
			#$stdev = sqrt(($with_tail_array[$i][6] - sq($with_tail_array[$i][5]*$with_tail_array[$i][1]))/($with_tail_array[$i][1]-1));	
		$stdev = nearest(0.0001, $stdev);
		#}
		#print "$with_tail_array[$i][0]$tail_display\t$with_tail_array[$i][1]\t$with_tail_array[$i][2]\t$with_tail_array[$i][3]\t".

		
		my @tail_length_dist = split(/;/, consolidate_count($with_tail_array[$i][6],$max_tail_length));
		
		#"$with_tail_array[$i][4]\t$with_tail_array[$i][5]\t$stdev\n";
		$mean = nearest(0.001, $mean);
		$sub_results .= "$with_tail_array[$i][0]$tail_display\t$with_tail_array[$i][1]\t$with_tail_array[$i][2]\t$with_tail_array[$i][3]\t".
		"$sum\t$mean\t$stdev\t$tail_length_dist[0]\t\t$tail_length_dist[1]\n";
	}
	print "Total read count: $sum_of_reads\n";
	$sub_results .= "Total read count: $sum_of_reads\n";
}


###############################
### Sub-sub routine section ###
###############################


### Within find_adaptor subroutine
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

## Within find_end subroutine
sub comp_end_score {

	my $seq1 = $_[0];
	my $test_length = length($Target_seq);
	my $score = 0;	

	if (length($seq1) < length($Target_seq)){
		$test_length = length($seq1);
	}

	if ($seq1 =~ m/$Target_seq$/) {
		$score = 45;
	} else {

		my @array_a = split(//, $seq1);
		my @array_b = split(//, $Target_seq);

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

	my $return = nearest(0.01, $score / length($Target_seq));
	return ($return);
}

## within final_collapse
sub sq {
	my $square = $_[0]*$_[0];
	return $square;
}

## within final_collapse
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

sub statcal {
	my @temp_data = split(/,/, $_[0]);
	my $count = @temp_data;
	my $sum = 0;
	my $var = 0;
	my $num = 0;
	foreach $num (@temp_data) { $sum += $num; }
	my $mean = $sum/$count;
	foreach $num (@temp_data) { $var += sq($num - $mean); }
	$var = $var/$count;
	my $stdev = sqrt($var);
	my $return = "$mean,$stdev,$sum";
	return ($return);
}

sub consolidate_count {
	my @count_matrix = split(/,/, $_[0]);
	my @report_matrix;
	my $report_matrix_count = 0;
	my $max_length = $_[1];
	#my $max_length = $_[1];
	for (my $i = 0; $i < @count_matrix; $i++){
		if ($report_matrix_count == 0){
			$report_matrix[$report_matrix_count][0] = $count_matrix[$i];  ## contains read length
			$report_matrix[$report_matrix_count][1] = 1;                 ## count at this read length  
			$report_matrix_count++;
		} elsif ($count_matrix[$i] == $count_matrix[$i-1]){
			$report_matrix[$report_matrix_count-1][1]++;
		} else {
			$report_matrix[$report_matrix_count][0] = $count_matrix[$i];  ## contains read length
			$report_matrix[$report_matrix_count][1] = 1; 
			$report_matrix_count++;
		} 
	}

	#for (my $i = 0; $i < $max_length; $i++){
	#	if ($report_matrix_count == 0){
	#		$report_matrix[$report_matrix_count][0] = $count_matrix[$i];  ## contains read length
	#		$report_matrix[$report_matrix_count][1] = 1;                 ## count at this read length  
	#		$report_matrix_count++;
	#	} elsif ($count_matrix[$i] == $count_matrix[$i-1]){
	#		$report_matrix[$report_matrix_count-1][1]++;
	#	} else {
	#		$report_matrix[$report_matrix_count][0] = $count_matrix[$i];  ## contains read length
	#		$report_matrix[$report_matrix_count][1] = 1; 
	#		$report_matrix_count++;
	#	} 
	#}

	#my $return_length = "";
	my $return_count = "";
	for (my $i = 1; $i < $max_length + 1 ; $i++){
		my $position_match_flag = 0;
		for (my $k = 0; $k < $report_matrix_count; $k++){
			if ($i == $report_matrix[$k][0]){
				$return_count .= "$report_matrix[$k][1]\t";
				$position_match_flag = 1;
				$k = $report_matrix_count + 2;
			}
		}
		if ($position_match_flag == 0){$return_count .= "0\t"};
	}
	#$return_length =~ s/\t$//;
	$return_count =~ s/\t$//;

	my $return_composite = '';
	for (my $i = 0; $i < $report_matrix_count ; $i++){
		$return_composite .= "$report_matrix[$i][0]($report_matrix[$i][1]),";
	}
	$return_composite =~ s/,$//;
	my $return = $return_composite.';'.$return_count;
	#print "$return\n";
	return ($return);
}

