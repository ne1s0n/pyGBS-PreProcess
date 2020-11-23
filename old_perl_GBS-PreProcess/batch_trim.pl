#!/usr/bin/perl

# Input should be a space/tab separated file containing both ends
# of the sample in separate files (skewer can't do interleaved)

my $usage = "$0 prepared_keyfile\n
 The prepared_keyfile should be a file containing the following data:

 Read1.fq Read2.fq Barcode1 Barcode2

 There are no headers, the order is the thing taken into account. The
script will run the trimming software with the illumina universal adapter
as the contaminant to remove, but will prepend the reverse complement of
the opposite barcode to it too (i.e., barcode2 when trimming read1, and vice
versa).

 The objective is to remove the barcode as well as the universal adapter,
leaving the cut site at the end (if any) clean and intact. You might then
remove more of the read if it happens to have more than one cut site (this
means it's chimeric), but you might not. And doing it this way ensures you
get the DNA that was sequenced, and that sequencing errors in the cut site
won't leave you with reads that have the barcode too.\n\n";

unless ($#ARGV == 0) {
    die $usage;
}

open (READ, $ARGV[0]) || die "Can't open list of files $ARGV[0]\n";
LOOP: while (<READ>) {
    if (/^\#/) {
        next LOOP;
    }
    chomp;
    my @x = split;
    if ($#x != 3) {
        die "Wrong number of columns, line $.\n";
    }
    my $out = `head -n 2 $x[0] | tail -n 1`;
    chomp($out);
    my $readlength = length($out);
    $out = `head -n 2 $x[1] | tail -n 1`;
    chomp($out);
    if (length($out) < $readlength) {
        $readlength = length($out);
    }

    # Length of reads has been established, we won't necessarily use it since
    # adapter removal might leave shorter reads, but it's possible if we
    # decide to, for example, remove the very short double cut tags.

    # Create the adapters containing the barcodes of the opposite read. Note
    # that we're using Illumina universal, so change the code here if you
    # want something else!

    $x[2] = &revcomp($x[2]);
    $x[2] .= "AGATCGGAAGAGC";
    $x[3] = &revcomp($x[3]);
    $x[3] .= "AGATCGGAAGAGC";

    system("trim_galore -q 0 --dont_gzip -e 0.1 --paired -a2 $x[2] -a $x[3] $x[0] $x[1]");
#    system("~/Downloads/TrimGalore-0.4.3/trim_galore -q 0 --dont_gzip -e 0.1 --paired $x[0] $x[1]");
    # Below is in case we want to use Skewer (doesn't work at all when
    # the adapter isn't at the very end, in the case of GBS this happens a
    # lot if you use a short base cutter)
    # The options are: Median Q of 30 or the read is thrown away, throw away any
    # read that's clipped (shouldn't clip much, but it does for some reason)
    # because we just want full length reads. Remember we want good quality
    # SNPs!

    #system("skewer -m pe -Q 30 -q 0 -t 16 -n -l $readlength --quiet $x[0] $x[1]");
}
close READ;

sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACTG/TGAC/;
    return $seq;
}
