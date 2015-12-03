#!/usr/bin/perl -w

# Author: Jay Mashl <rmashl@genome.wustl.edu>
# Date: 2015-02-24
# This file is part of BreakRUsT.
#
# version 1: empirical filter for removal of duplicate or near-duplicate Tigra calls. Subsequent calls that are within $tol base pairs of the first call are rejected. Confidence intervals removed.
# version 2: Check SV type and select best quality score.
# version 3: convert TRA translocations to mates
# version 4: user can provide output filenames; add final quality score filter

# Tigra scores can be revised to a hard-coded value (e.g., 60). See code. 

use strict;
use Getopt::Long;


my $tol=10000;
my $qcutoff=0;
my $hadFirst=0;
my ($lastChr, $lastPos);
my ($lastChr2, $lastEnd);
my @pass=();
my $pass_cnt=0;

my $pass_file="";
my $fail_file="";
my $tigra_file="";


my $opt_result;
$opt_result = GetOptions (
    'tigra_file=s' => \$tigra_file,
    'pass_file=s' => \$pass_file,
    'fail_file=s' => \$fail_file,
    'q=i' => \$qcutoff,
    'tol=i'       => \$tol,
    ) || die "Error in command line arguments\n";

unless($opt_result) {
    die help_text();
}
unless ($tigra_file) {   warn "\nYou must provide a filename from which to read tigra variants\n";   die help_text(); }
unless ($pass_file) {   warn "\nYou must provide a filename to which to write variants that pass clustering\n";   die help_text(); }
unless ($fail_file) {   warn "\nYou must provide a filename to which to write variants that fail clustering\n";   die help_text(); }


print "\nClustering variants within tolerance of $tol bp...\n\n";

open(IN, "< $tigra_file");
open (PASS, ">  $pass_file" );
open (FAIL, ">  $fail_file" );

while (my $line=<IN>) {
    if($line =~ /^#/) {
	print PASS $line unless ($line =~ /CIPOS/ || $line =~ /ID=GT/);
	next;
    }
    chomp $line;
    my($this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/,$line;
    my($this_chr2, $this_pos2);
    # Get chr2 from INFO field
    my @a=split /\;/,$info;

    #  map to hash won't work for our data
    foreach (@a) {
        if ($_=~ /^CHR2=/) {
            my @b=split /=/,$_;
            $this_chr2=$b[1];
        } elsif ($_=~ /^END=/) {
            my @b=split /=/,$_;
            $this_pos2=$b[1];
        }
    }

#	print "Read:  $this_chr $this_pos $this_chr2 $this_pos2\n";
    
    if ($hadFirst==0) {  # first line always passes
	$pass[$pass_cnt]{line}   = $line;
	$pass[$pass_cnt]{chr1}   = $this_chr;
	$pass[$pass_cnt]{pos1}   = $this_pos;
	$pass[$pass_cnt]{chr2}   = $this_chr2;
	$pass[$pass_cnt]{pos2}   = $this_pos2;
	$pass[$pass_cnt]{svtype} = $alt;
	$pass[$pass_cnt]{qual}   = $qual;
	++$pass_cnt;
	$hadFirst=1;
#	print "Default retain: $this_chr $this_pos $this_chr2 $this_pos2\n";
    }
    
    else {
	my $reject=0;
	for (my $k=0; $k < $pass_cnt; ++$k) {
	    my ($d1,$d2);
#	    print "Test: ( $k )  $pass[$k]{chr1} $pass[$k]{pos1} $pass[$k]{chr2} $pass[$k]{pos2} \n";





	    if( ($this_chr eq $pass[$k]{chr1}) && ($this_chr2 eq $pass[$k]{chr2}) &&	( $alt eq $pass[$k]{svtype} )    ) {
		$d1=abs(($this_pos+0) - ($pass[$k]{pos1}+0));
		$d2=abs(($this_pos2+0) - ($pass[$k]{pos2}+0));
		if( $d1 <= $tol && $d2 <= $tol) {  # not too different
		    if ($qual > $pass[$k]{qual}) { 
			# Reject old, keep new
			print FAIL $pass[$k]{line}."\n";
			$pass[$k]{line}   = $line;			$pass[$k]{chr1}   = $this_chr;
			$pass[$k]{pos1}   = $this_pos;			$pass[$k]{chr2}   = $this_chr2;
			$pass[$k]{pos2}   = $this_pos2;			$pass[$k]{svtype} = $alt;
			$pass[$k]{qual}   = $qual;
		    } else {
			# Reject
			print FAIL $line."\n";
		    }
		    $reject=1;
		    last;
		}

	    } elsif( ($this_chr eq $pass[$k]{chr2}) && ($this_chr2 eq $pass[$k]{chr1}) &&    ( $alt eq $pass[$k]{svtype} ) ) {
		$d1=abs(($this_pos+0) - ($pass[$k]{pos2}+0));
		$d2=abs(($this_pos2+0) - ($pass[$k]{pos1}+0));
		if( $d1<=$tol && $d2<=$tol) { 
		    if ($qual > $pass[$k]{qual}) { 
			# Reject old, keep new
			print FAIL $pass[$k]{line}."\n";
			$pass[$k]{line}   = $line;			$pass[$k]{chr1}   = $this_chr;
			$pass[$k]{pos1}   = $this_pos;			$pass[$k]{chr2}   = $this_chr2;
			$pass[$k]{pos2}   = $this_pos2;			$pass[$k]{svtype} = $alt;
			$pass[$k]{qual}   = $qual;
		    } else {
			# Reject
			print FAIL $line."\n";
		    }
		    $reject=1;
		    last;
		}
	    }
	}  # for
	    
	if( !$reject ) {
	    $pass[$pass_cnt]{line}   = $line;		$pass[$pass_cnt]{chr1}   = $this_chr;
	    $pass[$pass_cnt]{pos1}   = $this_pos;		$pass[$pass_cnt]{chr2}   = $this_chr2;
	    $pass[$pass_cnt]{pos2}   = $this_pos2;		$pass[$pass_cnt]{svtype} = $alt;
	    $pass[$pass_cnt]{qual}   = $qual;
	    ++$pass_cnt;
	}

    }
}  #while

my $mate_id=1;

for (my $k=0; $k < $pass_cnt; ++$k) {
    my($this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/, $pass[$k]{line};
    
    # Final quality filter on the deemed best
    if( $qual <  $qcutoff ) {
	$filter = "low_quality";
	print FAIL join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
    } else {

	my($this_chr2, $this_pos2);
	# Get chr2 from INFO field
	my @a=split /\;/,$info;
	foreach (@a) {
	    if ($_=~ /^CHR2=/) {
		my @b=split /=/,$_;
		$this_chr2=$b[1];
	    } elsif ($_=~ /^END=/) {
		my @b=split /=/,$_;
		$this_pos2=$b[1];
	    }
	}
	# Expand translocations to mate pairs
	if( $info =~ /SVTYPE=TRA/) {
	    # Tigra-EXT is hard-coded with +/-10 bp. Increasing by 100bp for read length.
	    # Reconstruct info string
	    $info = "SVTYPE=BND;MATEID=".$mate_id."_2;IMPRECISE;CIPOS=-110,110;CIEND=-110,110;SOMATIC";
	    $id=$mate_id."_1";
	    $ref=".";
	    $alt=$this_chr2.":".$this_pos2;
	    print PASS  join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info)."\n";
	    $info = "SVTYPE=BND;MATEID=".$mate_id."_1;IMPRECISE;CIPOS=-110,110;CIEND=-110,110;SOMATIC";
	    $id=$mate_id."_2";
	    $ref=".";
	    $alt=$this_chr.":".$this_pos;
	    print PASS  join("\t", $this_chr2,$this_pos2,$id,$ref,$alt,$qual,$filter,$info)."\n";
	    $mate_id++;
	}
	else {
	    # Tigra-EXT is hard-coded with +/-10 bp. Increasing by 100bp for read length.
	    $info =~ s/CIPOS=\-10,10\;CIEND=\-10,10\;/CIPOS=\-110,110\;CIEND=\-110,110\;/;
	    $id=".";
	    $alt=".";
	    # Optional: hard code value
	    #    $qual="60";
	    print PASS  join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info)."\n";
	}
    }
    
} # for

close(FAIL);
close(PASS);
close(IN);

print "DONE\n\n";


sub help_text {
    return <<HELP;

SYNOPSIS
BP_cluster_filter      [options]  

1. Remove effectively redundant calls within an interval. 
2. Converts translocations to VCF mate pairs (but does not yet provide reference base or strand information)


OPTIONS
--tigra_file           input tigra file, usually filtered
--pass_file            output file containing variants that passed filtering
--fail_file            output file containing variants that failed filtering
--q                    (optional) Tigra quality score cutoff  [0]
--tol                  (optional) variants with same sv type within this distance of the first such call are considered the same  [10000]


--help                 this message

AUTHOR
Jay Mashl              Original code

HELP
}


