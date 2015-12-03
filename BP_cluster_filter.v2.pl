#!/usr/bin/perl -w

# Author: Jay Mashl <rmashl@genome.wustl.edu>
# Date: 2015-02-22
# This file is part of BreakRUsT.
#
# version 1: empirical filter for removal of duplicate or near-duplicate Tigra calls. Subsequent calls that are within $tol base pairs of the first call are rejected. Confidence intervals removed.
# version 2: Check SV type and select best quality score.

# Tigra scores can be revised to a hard-coded value (e.g., 60). See code. 

use strict;
use warnings;

my $tol=10000;


my $hadFirst=0;

my ($lastChr, $lastPos);
my ($lastChr2, $lastEnd);


my @pass=();
my $pass_cnt=0;



open (PASS, ">", "BP_cluster.v2.PASS" );
open (FAIL, ">", "BP_cluster.v2.FAIL" );

while (my $line=<>) {
    if($line =~ /^#/) {
	print PASS $line unless ($line =~ /CIPOS/ || $line =~ /ID=GT/);
	next;
    }
    chomp $line;
    my($this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/,$line;
    my($this_chr2, $this_pos2);
    # Get chr2 from INFO field
    my @a=split /\;/,$info;
    foreach (@a) {
	my @b=();
	if ($_=~m/^CHR2=/) {
	    @b=split /=/,$_;
	    $this_chr2=$b[1];
	} elsif ($_=~m/^END=/) {
	    @b=split /=/,$_;
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

# Note: Seems like read pairs of a tandem duplication could appear to correspond to a small INS. Allow Tigra to override BreakDancer.
	    # A DUP that looks like an INS would be called INS (per Kai).
#		( $alt eq $pass[$k]{svtype} || ($alt=~/(INS|DUP)/ && $pass[$k]{svtype}=~/(INS|DUP)/)) ) {

	    if( ($this_chr eq $pass[$k]{chr1}) && ($this_chr2 eq $pass[$k]{chr2}) &&
		( $alt eq $pass[$k]{svtype} )    ) {
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

	    } elsif( ($this_chr eq $pass[$k]{chr2}) && ($this_chr2 eq $pass[$k]{chr1}) &&
		     ( $alt eq $pass[$k]{svtype} ) ) {
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


for (my $k=0; $k < $pass_cnt; ++$k) {
    
    my($this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/, $pass[$k]{line};
    $id=".";
    $alt=".";
    # Optional: hard code value
    #    $qual="60";

    # Tigra-EXT seems hard-coded with +/-10 bp. Increasing by 100bp for read length.
    $info =~ s/CIPOS=\-10,10\;CIEND=\-10,10\;/CIPOS=\-110,110\;CIEND=\-110,110\;/;

    print PASS  join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info)."\n";

}

close(FAIL);
close(PASS);
