#!/usr/bin/perl

use strict;

use Getopt::Long;

use threads;
use threads::shared;
use Thread::Semaphore;

use IO::CaptureOutput qw/capture/;

# Input should be a space/tab separated file containing both ends
# of the sample in separate files (skewer can't do interleaved)

my $usage = "Usage: $0 [parameters] -i input_keyfile

Optional parameters:
 -d Dry run, print a list of commands that will be run but
    run nothing. This is useful if you want to send the tasks
    to your favourite computer cluster or task manager.
 -g Create compressed GZIP output (default: Don't GZIP);
 -h This help menu.
 -t [int] Number of threads (default = 1).


 The prepared_keyfile should be a file containing the following data:

 Read1.fq Barcode1 Barcode2

 There are no headers, the order is the thing taken into account. You
can add headers if you want with a hash # symbol at the start for human-
readable purposes, but any line starting with a # will be ignored (so
you can put comments in the middle of the keyfile if you so wish).

 The script will run the trimming software with the illumina universal
adapter as the contaminant to remove, but will prepend the reverse
complement of the opposite barcode to it too (i.e., barcode2 when
trimming read1).

 The objective is to remove the barcode as well as the universal adapter,
leaving the cut site at the end (if any) clean and intact. You might then
remove more of the read if it happens to have more than one cut site (this
means it's chimeric), but you might not. And doing it this way ensures you
get the DNA that was sequenced, and that sequencing errors in the cut site
won't leave you with reads that have the barcode too.\n\n";

# Parse the command line
my $help;
my $threads = 1;
my $keyfile;
my $gzip;
my $dry_run;
GetOptions ('h' => \$help,
            't=i' => \$threads,
            'i=s' => \$keyfile,
            'g' => \$gzip,
            'd' => \$dry_run);
if ($help) {
    die $usage;
}
unless (-e $keyfile) {
    die "Error: Keyfile not found!\n\nUse $0 -h for help";
}



my @commands;
my @output :shared;

open (READ, $keyfile) || die "Can't open keyfile $keyfile\n";
LOOP: while (<READ>) {
    if (/^\#/) {
        next LOOP;
    }
    chomp;
    my @x = split;
    if ($#x != 2) {
        die "Wrong number of columns, line $.\n";
    }
    my $out;
    if ($x[0] =~ /\.gz/) {
        $out = `zcat $x[0] | head -n 2 | tail -n 1`;
    }
    else {
        $out = `head -n 2 $x[0] | tail -n 1`;
    }
    chomp($out);

    my $readlength = 0;

    if (length($out) < $readlength) {
        $readlength = length($out);
    }

    # Length of reads has been established, we won't necessarily use it since
    # adapter removal might leave shorter reads, but it's possible if we
    # decide to, for example, remove the very short double cut tags.

    # Create the adapters containing the barcodes of the opposite read. Note
    # that we're using Illumina universal, so change the code here if you
    # want something else!

    $x[1] = &revcomp($x[1]);
    $x[1] .= "AGATCGGAAGAGC";
    $x[2] = &revcomp($x[2]);
    $x[2] .= "AGATCGGAAGAGC";

    my $runfile = "trim_galore -q 0 -e 0.1 -a $x[2] ";
    if ($gzip) {
	$runfile .= "--gzip ";
    }
    $runfile .= $x[0];

    push (@commands, $runfile);
    
#    system("trim_galore -q 0 --gzip -e 0.1 -a $x[2] $x[0]");
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

if ($dry_run) {
    for my $i (@commands) {
	print $i, "\n";
    }
    exit;
}

my $threadcount = -1;
my $sem = Thread::Semaphore->new($threads);
my @threads = map{
    $sem->down;
    $threadcount++;
    threads->create(\&run_command, ($_, $threadcount))
} @commands;

print @output;

sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACTG/TGAC/;
    return $seq;
}

sub run_command {
    my ($job, $i) = @_;
    sleep(int(rand(2)));
    my ($stdout, $stderr);
    capture sub {
	system("$job");
    } => \$stdout, \$stderr;
    $output[$i] = $stdout;
    $output[$i] .= $stderr;
    $sem->up; #Slot released and available
}
