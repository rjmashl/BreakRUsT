#!/usr/bin/perl -w

# Author: Jay Mashl <rmashl@genome.wustl.edu>
# Date: 2015-02-05
# This file is part of BreakRUsT.
#
# Tigra-EXT filter (based on 20141208-downloaded version)
#
# We assume input to Tigra was breakdancer CTX, so we should consider only <TRA> in tigra-ext output.

use strict;
use warnings;

use Getopt::Long;

my @bddb=();
#Defaults
my $sv_file;
my $pass_file;
my $fail_file;
my $bd_file;
my $help;
my $bd_strand=0;
my $match_dist=1000;


my $opt_result;
$opt_result = GetOptions (
    'sv_file=s'          => \$sv_file,
    'pass_file=s'        => \$pass_file,
    'fail_file=s'        => \$fail_file,
    'match_distance=i'   => \$match_dist,
    'breakdancer_file=s' => \$bd_file,
    'strand_filter'      => \$bd_strand, 
    'help'               => \$help,
    ) || die "Error in command line arguments\n";

unless($opt_result) {
    die help_text();
}
if($help){
    print STDOUT help_text();
    exit 0;
}
unless ($sv_file) {
  warn "\nYou must provide an sv file\n";
  die help_text();
}

unless ($pass_file) {
  warn "\nYou must provide a filename for variants that pass filtering\n";
  die help_text();
}

unless ($fail_file) {
  warn "\nYou must provide a filename for variants that fail filtering\n";
  die help_text();
}

unless ($bd_file) {
  warn "\nYou must provide a breakdancer results files\n";
  die help_text();
}


print STDOUT "\nStrand filter will".(($bd_strand)?(""):(" not"))." be applied.\n";
print STDOUT "\n";

#read_breakdancer
open (BD, "<", $bd_file ) || die "Error: cannot open $bd_file";
while (my $line=<BD>) {
    next if($line =~ /^#/);
    my($chr1,$pos1,$reads1,$chr2,$pos2,$reads2) = split /\t/,$line;
    my $tmp="$chr1\t$pos1\t$reads1\t$chr2\t$pos2\t$reads2";
    my $href={ chr1 => $chr1, pos1 => $pos1, reads1 => $reads1, chr2 => $chr2, pos2 => $pos2, reads2 => $reads2, line => $tmp};
    push (@bddb, $href);
}
close(BD);









open (INPUT, "<", $sv_file ) || die "Error: cannot open $sv_file";
open (PASS, ">", $pass_file ) || die "Error: cannot open $pass_file";
open (FAIL, ">", $fail_file ) || die "Error: cannot open $fail_file";

while (my $line=<INPUT>) {
    if($line =~ /^#/) {
	print PASS $line ;
	next;
    }
    chomp $line;
    my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/,$line;
    
    # FILTER: Wrong SV type
    if ($alt !~ m/TRA/  || $tigra_contig !~ m/CTX/ || $info !~ m/SVTYPE=TRA/) {
	$filter="bad_sv_type";
	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	
    } else {
	# Get chr2 from INFO field
	my @a=split /\;/,$info;
	my $chr2="";
	my $pos2="";
	foreach (@a) {
	    if ($_=~m/^CHR2=/) {
		my @b=split /=/,$_;
		$chr2=$b[1];
	    } elsif ($_=~m/^END=/) {
		my @b=split /=/,$_;
		$pos2=$b[1];
	    }
	}
	# Get chrs from tigra contig
	my ($contig_chr1,$contig_pos1,$contig_chr2,$contig_pos2) = split /\^/,$tigra_contig;
	
	# FILTER: reject unmapped chroms.
	if ( $chr !~ m/^[1-9XY]/ || $chr2 !~ m/^[1-9XY]/ || $contig_chr1 !~ m/^[1-9XY]/ || $contig_chr2 !~ m/^[1-9XY]/ ) {
	    $filter="unmapped_chr";
	    print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	}
	else {
	    
	    # FILTER: reject mismatched chromosomes
	    if( (($chr ne $contig_chr1)||($chr2 ne $contig_chr2)) && (($chr ne $contig_chr2)||($chr2 ne $contig_chr1))  ) {
		    $filter="chr_mismatch";
		    print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	    }
	    
	    else { 

		my @entries=();
		match_breakdancer(\@entries, $chr,$pos,$chr2,$pos2);

		# FILTER: tigra break points too different from breakdancer break points
		if (scalar(@entries) < 1) {
		    $filter="breakdancer_miss";
		    print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
		}

		else { 
		    
		    foreach my $v (@entries) {
#			print PASS "# entries = ".scalar(@entries)."\n";

			if ($bd_strand==1) {
			    #Get BD strand readcounts
			    my ($fwd1,$rev1,$fwd2,$rev2);
			    if($bddb[$v]{reads1} =~ m/(\d+)\+(\d+)\-/i) { ($fwd1,$rev1) = ($1,$2); }
			    if($bddb[$v]{reads2} =~ m/(\d+)\+(\d+)\-/i) { ($fwd2,$rev2) = ($1,$2); }
			    if ($fwd1==0 || $rev1==0 || $fwd2==0 || $rev2 == 0) {
				$filter="strandedness";
				print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig, $bddb[$v]{line})."\n";
			    }
			    else {  #PASS
				$filter="PASS";
				print PASS join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig, $bddb[$v]{line})."\n";
			    }
			}
			else {  #PASS
				$filter="PASS";
			    print PASS join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig, $bddb[$v]{line})."\n";
			}
			
		    }
		    
		}
	    }
	    
	}
	
    }
    
}  #while


close(FAIL);
close(PASS);
close(INPUT);


sub help_text {
    return <<HELP;

SYNOPSIS
tigra-ext_filter.pl    [options]

OPTIONS
--sv_file              VCF file from Tigra-sv
--breakdancer_file     breakdancer file
--match_distance       maximum distance in bp that tigra breakpoints may be from breakdancer breakpoints [default: 1000]
--strand_filter        filters out calls lacking reads on either strand  [default: off]
--pass_file            output file containing variants that passed filtering
--fail_file            output file containing variants that failed filtering
--help                 this message

AUTHOR
Jay Mashl              Original code

HELP
}


sub match_breakdancer {
    my ($p_entries,$gchr1,$gpos1,$gchr2,$gpos2) = @_;


    for(my $k=0; $k<scalar(@bddb); $k++) {
	if( ($gchr1 eq $bddb[$k]{chr1}) && ($gchr2 eq $bddb[$k]{chr2})) {
	    if( abs($gpos1 - $bddb[$k]{pos1}) <= $match_dist  &&  abs($gpos2 - $bddb[$k]{pos2}) <= $match_dist) {
		push(@$p_entries, $k);
	    }
	} elsif (($gchr1 eq $bddb[$k]{chr2}) && ($gchr2 eq $bddb[$k]{chr1})) {
	    if( abs($gpos1 - $bddb[$k]{pos2}) <= $match_dist  &&  abs($gpos2 - $bddb[$k]{pos1}) <= $match_dist) {
		push(@$p_entries, $k);
	    }
	}

    }
}

