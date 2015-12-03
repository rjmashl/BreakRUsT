#!/usr/bin/perl -w

# Author: Jay Mashl <rmashl@genome.wustl.edu>
# Creation Date: 2015-01-23
# Last Updated: 2015-02-16
# This file is part of BreakRUsT. 
#
# empirical filter for removal of duplicate or near-duplicate Tigra calls. Subsequent calls that are within $tol base pairs of the first call are rejected. Disadvantage is that subsequent putative calls with higher scores get marked as FAIL. Will fix this later.

# Tigra scores can be revised to a hard-coded value (e.g., 60). See code. 

use strict;
use warnings;

my $tol=10000;


my $hadFirst=0;

my ($lastChr, $lastPos);
my ($lastChr2, $lastEnd);


my @pass=();
my $pass_cnt=0;



open (PASS, ">", "BP_cluster.PASS" );
open (FAIL, ">", "BP_cluster.FAIL" );

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

#    $line2 = join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info);

#	print "Read:  $this_chr $this_pos $this_chr2 $this_pos2\n";

    
    if ($hadFirst==0) {  # first line always passes
	$pass[$pass_cnt]{line} = $line;
	$pass[$pass_cnt]{chr1} = $this_chr;
	$pass[$pass_cnt]{pos1} = $this_pos;
	$pass[$pass_cnt]{chr2} = $this_chr2;
	$pass[$pass_cnt]{pos2} = $this_pos2;
	++$pass_cnt;
	$hadFirst=1;
#	print "Default retain: $this_chr $this_pos $this_chr2 $this_pos2\n";
    }

    else {
	
	my $reject=0;
	for (my $k=0; $k < $pass_cnt; ++$k) {
	    my ($d1,$d2);
#	    print "Test: ( $k )  $pass[$k]{chr1} $pass[$k]{pos1} $pass[$k]{chr2} $pass[$k]{pos2} \n";
	    if( ($this_chr eq $pass[$k]{chr1}) && ($this_chr2 eq $pass[$k]{chr2})  ){  # if chrs same in this order
#		print "chrs same as is\n";
         	$d1=abs(($this_pos+0) - ($pass[$k]{pos1}+0));
		$d2=abs(($this_pos2+0) - ($pass[$k]{pos2}+0));
         	if( $d1 <= $tol && $d2 <= $tol) {  # too different then keep
#		    print "REJECT: $this_chr $this_pos $this_chr2 $this_pos2  dists $d1 $d2\n";
		    $reject=1;
		    last;
		}
	    } else {
		if( ($this_chr eq $pass[$k]{chr2}) && ($this_chr2 eq $pass[$k]{chr1})  ){  # if chrs same in opposite order
#		    print "chrs opp order\n";
		    $d1=abs(($this_pos+0) - ($pass[$k]{pos2}+0));
		    $d2=abs(($this_pos2+0) - ($pass[$k]{pos1}+0));
		    if( $d1<=$tol && $d2<=$tol) {  # too different then keep
#			print "REJECT: $this_chr $this_pos $this_chr2 $this_pos2  dists $d1 $d2\n";
			$reject=1;
			last;
		    }
		}
		else { # at least one chr is different; so keep
		    $reject=0;
		}
	    }
	    
	}
	if ($reject) {
#	    print "Rejected\n";
	    print FAIL $line."\n";
	} else {
#	    print "Retained\n";
	    $pass[$pass_cnt]{line} = $line;
	    $pass[$pass_cnt]{chr1} = $this_chr;
	    $pass[$pass_cnt]{pos1} = $this_pos;
	    $pass[$pass_cnt]{chr2} = $this_chr2;
	    $pass[$pass_cnt]{pos2} = $this_pos2;
	    ++$pass_cnt;
	}

    }
#    print "---------------------------------------------------------\n";
#

}  #while

for (my $k=0; $k < $pass_cnt; ++$k) {
    
    my($this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/, $pass[$k]{line};
    $id=".";
    $alt=".";
#    $qual="60";
    $info =~ s/CIPOS=\-10,10\;CIEND=\-10,10\;//;

    print PASS  join("\t", $this_chr,$this_pos,$id,$ref,$alt,$qual,$filter,$info)."\n";

}

close(FAIL);
close(PASS);
