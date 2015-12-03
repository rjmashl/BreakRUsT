#!/usr/bin/perl -w

# Author: Jay Mashl <rmashl@genome.wustl.edu>
# Date: 2015-02-24
# This file is part of BreakRUsT.
#
# Tigra-EXT filter (based on 20141208-downloaded version)
# 
# Fixed: The Tigra-ext vcf reassembled SVtype often does not match that of the original BD call on which it was based. 
#        In a sense, Tigra can reassemble the wrong SV. To account for this, we can simply
#        require the original BD input file to contain a call that is sufficiently close to the reassembled call.
#

use strict;
use Getopt::Long;

my @bddb=();
my %seen=();
#Defaults
my $sv_file;
my $pass_file;
my $fail_file;
my $bd_file;
my $help;
my $bd_strand=0;
my $match_dist=1000;
my $cnt_int=5000;
my $key="";
my $bMappedOnly=0;
my $lib_key="";
my $min_libreads=0;
my $keycol=-1;
my $bKeyCol=0;


my $opt_result;
$opt_result = GetOptions (
    'sv_file=s'          => \$sv_file,
    'pass_file=s'        => \$pass_file,
    'fail_file=s'        => \$fail_file,
    'match_distance=i'   => \$match_dist,
    'breakdancer_file=s' => \$bd_file,
    'strand_filter'      => \$bd_strand, 
    'mapped_only'        => \$bMappedOnly, 
    'keycol=i'           => \$keycol,
    'lib_key=s'          => \$lib_key,
    'min_libreads=i'     => \$min_libreads,
    'help'               => \$help,
    ) || die "Error in command line arguments\n";

unless($opt_result) {
    die help_text();
}
if($help){    print STDOUT help_text();    exit 0; }
unless ($sv_file) {   warn "\nYou must provide an sv file\n";   die help_text();  }
unless ($pass_file) {   warn "\nYou must provide a filename for variants that pass filtering\n";   die help_text(); }
unless ($fail_file) {   warn "\nYou must provide a filename for variants that fail filtering\n";   die help_text(); }
unless ($bd_file) {   warn "\nYou must provide a breakdancer results files\n";  die help_text(); }
# Set filter options
if ($lib_key ne "" && $keycol > -1) {   warn "\nSpecify only one of --lib_key or --keycol\n";   die help_text(); }

if (($lib_key ne "" || $keycol > -1)) {
    if ($min_libreads < 1) {   
	warn "\nNumber of --min_libreads not possible\n";   die help_text(); 
    }
} else {
    warn "\nNeither --lib_key nor --keycol provided. No libread filtering will be performed.\n";
    $min_libreads=0; #enforce for flagging
}


print STDOUT "\nStrand filter will".(($bd_strand)?(""):(" not"))." be applied.\n";
print STDOUT "\n";

#read_breakdancer
my $lc=`wc -l < $bd_file`;
open (BD, "<", $bd_file ) || die "Error: cannot open $bd_file";
print STDOUT "\nReading BreakDancer file $bd_file ...\n";
my $cnt=0;
my @headers=();
my @colname=();

while (my $line=<BD>) {
    chomp $line;
    $cnt++;
    printf STDOUT "%d %%\n", int(($cnt/$lc)*100.0) if ($cnt % $cnt_int == 0);
    if($line =~ /^#/) {
	if($line =~ /^#Chr/) {
	    @headers = split /\t/, $line;
	    if ($keycol >= 0 ) { # extract string from col  0 -> 11, 1 -> 12, ....		
		$lib_key = $headers[ $keycol+11 ];
	    } else {  # determine col given string
		for(my $i=11; $i<scalar(@headers); $i++) {
		    $keycol = $i-11 if ($headers[$i] =~ /$lib_key/);
		}
	    }
	    print STDOUT "Column $keycol with header $lib_key will be used for library read suppport.\n";
	    print STDOUT "Reading BD file...\n";
	    $bKeyCol=1;  #ensure flag is set for filtering
	}
	next;
    }

    
    if (($lib_key ne ""  && $keycol < 0) || ($lib_key eq "" && $keycol >= 0) ) {
	die "\nAn error has occured. No header found. Please check BD input file for missing #Chr line\n\n";
    } 

    my($chr1,$pos1,$reads1,$chr2,$pos2,$reads2,$vartype,$bdsize,$bdscore,$bdnumreads,$bdlibreads) = split /\t/,$line;
    my $tmp="$chr1\t$pos1\t$reads1\t$chr2\t$pos2\t$reads2\t$vartype";
    # Extract library num reads for the requested library (string match is the only possbility)
    my $libreads;
    my @a = split /\:/, $bdlibreads;
    for( my $i=0; $i<scalar(@a); $i++) {
	if( $a[$i] =~ /$lib_key/) {
	    my @b = split /\|/, $a[$i];
	    $libreads = $b[1];
#	    print STDOUT "$chr1 $pos1 has $libreads lib reads.\n";
	    last;
	}
    }
    my $href={ chr1 => $chr1, pos1 => $pos1, reads1 => $reads1, chr2 => $chr2, pos2 => $pos2, reads2 => $reads2, vartype => $vartype, numreads => $bdnumreads,  libreads => $libreads, line => $tmp};
#    print STDOUT $tmp."\n";
    push (@bddb, $href);
}
printf STDOUT "%d %%\n", 100;
close(BD);
printf STDOUT "There were ".scalar(@bddb)." BreakDancer entries read.\n";

# Main
$lc=`wc -l < $sv_file`;
chomp $lc;
open (INPUT, "<", $sv_file ) || die "Error: cannot open $sv_file";
open (PASS, ">", $pass_file ) || die "Error: cannot open $pass_file";
open (FAIL, ">", $fail_file ) || die "Error: cannot open $fail_file";
print STDOUT "\nProcessing $lc lines of Tigra file $sv_file ...\n";
$cnt=0;
while (my $line=<INPUT>) {
    $cnt++;
    printf STDOUT "%d %%\n", int(($cnt/$lc)*100.0) if ($cnt % $cnt_int == 0);
    if( $line =~ /^#/ ) {
	print PASS $line ;
	next;
    }
#    print STDOUT $line."\n";
   chomp $line;
    my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig) = split /\t/,$line;
    $alt =~ s/<//;
    $alt =~ s/>//;
    # Get chr2 from INFO field
    my @a=split /\;/,$info;
    my $chr2="";
    my $pos2="";
# Works only if all member can split on /= /
#    my %h = map {split /=/} @a;
#    $chr2 = $h{'CHR2'};
#    $pos2 = $h{'END'};
    foreach (@a) {
	if ($_=~ /^CHR2=/) {
	    my @b=split /=/,$_;
	    $chr2=$b[1];
	} elsif ($_=~ /^END=/) {
	    my @b=split /=/,$_;
	    $pos2=$b[1];
	}
    }

    # Get chrs from tigra contig
    my ($contig_chr1,$contig_pos1,$contig_chr2,$contig_pos2,$contig_type) = split /\^/,$tigra_contig;

    # FILTER: reject calls that did not assemble
# The following line did not overspecify, but it did have a redundant test
#    if( $contig_type=~/(DEL|INS|INV)/ && $info =~ /UNKNOWN_LEN/) {
    if( $info =~ /UNKNOWN_LEN/) {
	$filter="assembly_failed";
	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	next;
    }

    # FILTER: reject unmapped chroms.
#    if ( $bMappedOnly==1 && ( $chr !~ /^[1-9XY]/ || $chr2 !~ /^[1-9XY]/ || $contig_chr1 !~ /^[1-9XY]/ || $contig_chr2 !~ /^[1-9XY]/ )) {
    if ( $bMappedOnly==1 && ( $chr !~ /^[1-9XY]/ || $chr2 !~ /^[1-9XY]/  )) {
	$filter="unmapped_chr";
	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	next;
    }
	
    # FILTER: inconsistent SV types. This may need more attention.
#    if( ($contig_type=~/INS/       && $alt!~/(INS|DUP)/) || 
#	($contig_type=~/(DEL|INV)/ && $alt!~/$contig_type/) || 
#	($contig_type=~/[CI]TX/    && $alt!~/TRA/) ) { 
#	$filter="sv_type_mismatch";
#	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
#	next;
#    } 
	    
    # FILTER: reject mismatched chromosomes
#    if( ($contig_type=~/(DEL|INS|INV|ITX)/ && ($chr ne $chr2 || $chr ne $contig_chr1 || $chr ne $contig_chr2)) ||
#	($contig_type=~/CTX/ && ($chr ne $contig_chr1 || $chr2 ne $contig_chr2) && ($chr ne $contig_chr2 || $chr2 ne $contig_chr1)) ) {
#	$filter="chr_mismatch";
#	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
#	next;
 #   }



    my @entries=();

    # FILTER: proximity of tigra breakpoints across all BD breakpoints
    for(my $k=0; $k<scalar(@bddb); $k++) {
	# tigra produces DEL,DUP,INS,INV,TRA.  These are to match BreakDancer CTX,ITX,DEL,INS,INV
	
	if( $alt eq $bddb[$k]{vartype} ) {  # DEL,INS,INV
	    if( ($chr eq $bddb[$k]{chr1}) &&  $chr eq $chr2 ) {
		if( abs($pos - $bddb[$k]{pos1}) <= $match_dist  &&  abs($pos2 - $bddb[$k]{pos2}) <= $match_dist) { #implicitly assumes pos < pos2
		    push(@entries, {idx => $k, dist =>  abs($pos - $bddb[$k]{pos1}) +  abs($pos2 - $bddb[$k]{pos2}), str_pass => "-1"});
		    next;
		}
	    }
	}

	if( ($alt eq "DUP" &&  $bddb[$k]{vartype} eq "INS") || ( $alt eq "TRA" && $bddb[$k]{vartype} eq "ITX" ) ) {
	    if( ($chr eq $bddb[$k]{chr1}) &&  $chr eq $chr2 ) {
		if( abs($pos - $bddb[$k]{pos1}) <= $match_dist  &&  abs($pos2 - $bddb[$k]{pos2}) <= $match_dist) { #implicitly assumes pos < pos2
		    push(@entries, {idx => $k, dist =>  abs($pos - $bddb[$k]{pos1}) +  abs($pos2 - $bddb[$k]{pos2}), str_pass => "-1"});
		    next;
		}
	    }
	}

	if( $alt eq "TRA" && $bddb[$k]{vartype} eq "CTX" ) {
	    if( ($chr eq $bddb[$k]{chr1}) && ($chr2 eq $bddb[$k]{chr2})) {
		if( abs($pos - $bddb[$k]{pos1}) <= $match_dist  &&  abs($pos2 - $bddb[$k]{pos2}) <= $match_dist) {
		    push(@entries, {idx => $k, dist =>  abs($pos - $bddb[$k]{pos1}) +  abs($pos2 - $bddb[$k]{pos2}), str_pass => "-1"});
		    next;
		} 
	    } elsif (($chr eq $bddb[$k]{chr2}) && ($chr2 eq $bddb[$k]{chr1})) {
		if( abs($pos - $bddb[$k]{pos2}) <= $match_dist  &&  abs($pos2 - $bddb[$k]{pos1}) <= $match_dist) {
		    push(@entries, {idx => $k, dist => abs($pos - $bddb[$k]{pos2}) + abs($pos2 - $bddb[$k]{pos1}), str_pass => "-1"});
		    next;
		}
	    }
	}

    }

    if (scalar(@entries) < 1) {
	$filter="no_matching_bd";
	print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	next;
    }

    # Get strand counts and record whether reads are balanced
    for (my $i=0; $i<scalar(@entries); $i++) {
	$entries[$i]{libreads} = $bddb[ $entries[$i]{idx} ]{libreads};

	if ($bd_strand) {
	    my ($fwd1,$rev1,$fwd2,$rev2);
	    if($bddb[ $entries[$i]{idx} ]{reads1} =~ /(\d+)\+(\d+)\-/) { ($fwd1,$rev1) = ($1,$2); }
	    if($bddb[ $entries[$i]{idx} ]{reads2} =~ /(\d+)\+(\d+)\-/) { ($fwd2,$rev2) = ($1,$2); }

	    $entries[$i]{str_data} = $fwd1."+".$rev1."-:".$fwd2."+".$rev2."-\n";

	    if ($fwd1==0 || $rev1==0 || $fwd2==0 || $rev2==0)  {
		$entries[$i]{str_pass} = 0;
	    }   else {
		$entries[$i]{str_pass} = 1;
	    }
	} else { 
	    $entries[$i]{str_pass} = 1;
	}
    }
    

    # FILTER: strand / check there is data to work with
    if ($bd_strand) {
	@entries = grep { $_->{str_pass} } @entries;
	if (scalar(@entries) < 1) {
	    $filter="strandedness";
	    print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	    next;
	}
    }

    # FILTER: meets sufficient lib reads
    if ($min_libreads > 0) {
	@entries = grep { $_->{libreads} >= $min_libreads } @entries ;
	if (scalar(@entries) < 1) {
	    $filter="min_libreads";
	    print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
	    next;
	}
    }
    
    # At this point, the Tigra call matches one or more BD calls.
#    @entries = sort { $a->{dist} <=> $b->{dist} } @entries;
#    my $bFoundMatch=0;
#    foreach my $v (@entries) {
#	$key = $chr."_".$pos."_".$alt."_".$chr2."_".$pos2."_";
#	if( ! exists( $seen{$key} )) {
#	    $filter="STATUS";
#	    $seen{$key}= { qual => $qual, line => join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)};
#	} else {
#	    if ($qual > $seen{$key}{qual}) {
#		$seen{$key}{line} =~ s/STATUS/score_was_bettered/;
#		print FAIL $seen{$key}{line};
#		$filter="STATUS";
#		$seen{$key}= { qual => $qual, line => join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)};
#	    } else {
#		$filter="score_unimproved";
#		print FAIL join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$gt,$gt_dat,$tigra_contig)."\n";
#	    }
#	}
#    } #foreach

    # Save the call;annotate info field with number of agreeing BD calls
    $filter="PASS";
    print PASS join("\t", $chr,$pos,$id,$ref,$alt,$qual,$filter, $info.";BDSUP=".scalar(@entries), $gt,$gt_dat,$tigra_contig)."\n";


    
} #while input

printf STDOUT "%d %%\nDone.\n", 100;
close(FAIL);


#foreach my $k (keys %seen) {
#    $seen{$k}{line} =~ s/STATUS/PASS/;
#    print PASS  $seen{$k}{line}."\n";
#}

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
--mapped_only          filters out chromosomes other than 1-22,X,Y

--lib_key | --keycol   (optional, select one) BD column tag (e.g., tumor bam filename) or 
                            support library column number (0-based)

--min_libreads         (optional) value to use for supporting reads in library corresponding to previous column 
                            (e.g., 3 supporting reads from tumor libraries used to make the call)

--pass_file            output file containing variants that passed filtering
--fail_file            output file containing variants that failed filtering
--help                 this message

AUTHOR
Jay Mashl              Original code

HELP
}


