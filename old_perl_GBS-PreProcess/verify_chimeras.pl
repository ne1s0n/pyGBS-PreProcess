#!/usr/bin/perl

my $usage = "$0 reads1 reads2 cutsite cutsite_remnant output_directory\n
 This script will go through pair-end GBS data and validate the pair-ends
themselves. It will assume data has been trimmed of adapter and reverse-
complement barcode, and that it shouldn't have any internal cut sites.

 If there is an internal cutsite, the reads will be considered chimeras, and
the two ends will be placed in an 'unpaired' list. Note that just because two
reads don't have any cut sites inside there's no guarantee they are not
chimeras. the ONLY way to be sure they're not chimeras is if the two reads
overlap without any cut sites.

 Use an overlap-detection tool on the clean pairs to separate those that overlap
from those that don't. Then use the following mapping stragegy:

1. Overlaps: Map in single-end mode (duh).
2. No overlap, no cutsites: Use in pair-end mode with independent pair mapping.
For example, NextGenMap offers the --fast-pair option, which will do just that.
This way, reads that are in the wrong place are reported as broken pairs.
3. Individual reads: As the name implies, map in single-end mode.

 The cutsite will be put in a regular expression, so it can be written as
a regex. For example, ApeKI remnant cut site would be GC[AT]GC. For the
remnant, it'll be C[AT]GC. Note that both cut site and remnant have to
be included even if they are the same!
\n";

unless ($#ARGV == 4) {
    die $usage;
}

open (R1, $ARGV[0]) || die;
open (R2, $ARGV[1]) || die;
my $cut = $ARGV[2];
my $rmn = $ARGV[3];
open (W1, ">$ARGV[4]/ver_$ARGV[0]") || die;
open (W2, ">$ARGV[4]/ver_$ARGV[1]") || die;
open (W3, ">$ARGV[4]/ver_$ARGV[0].singles.fq") || die;

my $discarded = 0;
my $paired;
my $unpaired;

LOOP: while (<R1>) {
    chomp;

    my $cutsites = 0;
    
    my @x = ($_);
    my @y;
    for (1..3) {
        my $i = <R1>;
        chomp ($i);
        push (@x, $i);
    }
    for (1..4) {
        my $i = <R2>;
        chomp ($i);
        push (@y, $i);
    }

    # First, we count the total number of cutsites, it shouldn't be more than
    # 1 per read! And we only want cut sites in the middle. First thing, check
    # that both reads have a remnant in the beginning. If not, discard them.

    if ($x[1] !~ /^$rmn/ || $y[1] !~ /^$rmn/) {
        $discarded++;
        next LOOP;
    }

    # Okay, so let's count the number of cutsites that both reads have!

    $cutsites += () = $x[1] =~ m/$cut/g;
    $cutsites += () = $y[1] =~ m/$cut/g;

    my $max_cutsites = 0;
    if ($cut eq $rmn) {
	# If the cut site is completely kept as the remnant, the regexp
	# above will match it at the very start of the read. In that case,
	# the minimum number of matches will be 2, because we've already
	# verified that both reads have the remnant in the beginning.
	$max_cutsites = 2;
    }

    if ($cutsites == $max_cutsites) {
        $paired++;
        for my $i (@x) {
	          print W1 "$i\n";
	      }
        for my $j (@y) {
	          print W2 "$j\n";
        }
    }
    else {
        $unpaired++;
        $x[1] =~ s/^($rmn.*?$cut).*/\1/;
        $y[1] =~ s/^($rmn.*?$cut).*/\1/; # The regexp _includes_ the full cutsite
	
	$x[1] =~ s/$rmn$//;
	$y[1] =~ s/$rnm$//; # We remove the remnant from the cutsite at the end
	                    # since it comes from the chimera. The only thing left
			    # should be the non-remnant recognition site bases

	# Note that proper non-chimeric GBS tags that are shorter than the
	# read length may have two cut sites (what they read depends on whether
	# the cut site is a palindrome), so the regexp above might or might
	# not catch them. If it does, they will STILL be placed
	# in the unpaired list, but that makes sense. FWD and REV are
	# identical, so pair-end mapping makes zero sense (and can break some
	# mapping software to boot). In principle the pair-end overlap might
	# be able to error correct the reads (take the base with the highest
	# score at each position), but that's beyond the scope of this script

        $x[3] = substr($x[3],0,length($x[1]));
        $y[3] = substr($y[3],0,length($y[1])); # Quality scores must be trimmed

        for my $i (@x) {
	          print W3 "$i\n";
        }
        for my $j (@y) {
	          print W3 "$j\n";
	}
    }
}

print STDERR "Proper Pairs: $paired\n";
print STDERR "Broken Pairs: $unpaired\n";
print STDERR "Discarded pairs: $discarded\n";
