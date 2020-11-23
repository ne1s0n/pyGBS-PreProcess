# GBS-PreProcess
Scripts to pre-process Genotyping-by-Sequencing data (GBS) before analysis with STACKS or,
potentially, any other GBS pipeline. The scripts will trim reads down by cut site and the
universal illumina adapter, taking into account the enzyme used and the barcodes.

The code was specifically made for trimming data from double-barcode GBS experiments. If
you use standard barcodes in only the first end (or just single-end reads) you can use the
data directly with TASSEL without fussing about with this software.

Prerrequisites:

    1. Perl, with the following modules:
        Getopt::Long;
        threads;
        threads::shared;
        Thread::Semaphore;
        IO::CaptureOutput;
      Most of these should be part of your standard installation. If not, add them with:
        cpan install [module name]
    2. TrimGalore (https://github.com/FelixKrueger/TrimGalore) installed somewhere in your $PATH
    3. TrimGalore itself requires cutadapt (https://github.com/marcelm/cutadapt),
       which can easily be installed with pip (pip install --user --upgrade cutadapt). If after
       the pip installation you can't run it, add ~/.local/bin/ to your $PATH.
    4. Trimgalore also requires FastQC under regular installation circumstances, but this
       wrapper doesn't generate fastqc reports so it can be skipped.

Other trimming software can be used, but the code will have to be modified. TrimGalore is
useful because it allows for the trimming based on different adapters for each end (which
the scripts take advantage of by using the barcodes as adapters), and is fast.

The scripts should be run in order:

**batch_trim.pl**: This program will take a list of fastq reads and barcodes, and wrap around
TrimGalore (which, in turn, is also a wrapper) to remove as much adapter as possible, as
well as the barcodes from the reverse complement of the reads (i.e., the barcode when the
GBS tag is shorter than the read length of the Illumina instrument). Trimming by Illumina
adapter alone can leave the barcode in place, and the barcode might not be picked up if the
reads are then trimmed by cutsite if there's sequencing errors in the multiple cutsites.

 The script requires a tab-separated "keyfile" containing the following information:

    sampleA_1.fq sampleA_2.fq barcode1 barcode2
    sampleB_1.fq sampleB_2.fq barcode1 barcode2
    (...)

 Headers are optional and only useful for human-reading purpose. The barcodes are read
independently, so they can be the same if the user wants to. The code will create a
cutadapt job for each pair of files with the right barcodes prepared for trimming and
with the Illumina universal adapter.

--**OR**--

**batch_trim_se.pl**: A single-end version of the batch_trim.pl script. It works in the
exact same way as batch_trim.pl: It will take a list of fastq reads and barcodes, and
wrap around TrimGalore (which, in turn, is also a wrapper) to remove as much adapter as
possible, as well as the barcodes from the reverse complement of the reads (i.e., the
barcode when the GBS tag is shorter than the read length of the Illumina instrument).
Trimming by Illumina adapter alone can leave the barcode in place, and the barcode might
not be picked up if the reads are then trimmed by cutsite of there's sequencing errors
in the multiple cutsites.

 The script requires a tab-separated "keyfile" containing the following information:

    sampleA_1.fq barcode1 barcode2
    sampleB_1.fq barcode1 barcode2
    (...)

 Headers are optional and only useful for human-reading purpose. The barcodes are read
independently, so they can be the same if the user wants to. The code will create a
cutadapt job for each pair of files with the right barcodes prepared for trimming and
with the Illumina universal adapter.

This script incorporates additional options:

    -t [int] Number of threads to run in parallel (default = 1)
    -g GZIP the output of TrimGalore
    -d Do a dry run (print all commands without running them)

 Barcode1 is not technically used, but it's useful to leave there for consistency.
 
 Eventually, both trimming scripts will be merged into one with additional options. This
is a work in progress.

**verify_chimeras.pl**: This program takes the output of the previous job and looks for reads
with multiple consecutive cutsites. run the script with no parameters to get the help
on what the software expects, but basically it wants:

verify_chimeras.pl reads1 reads2 cutsite cutsite_remnant output_directory

The first two parameters (reads1 and reads2) are self-explanatory, the Illumina reads coming
out of the trimming.

The first cutsite is the full **cutsite**, not the **remnant** of the enzyme cut site,
and for the sake of making things easier, it takes the form of a perl type regular
expression. So, for example, ApeKI is GC\[AT\]GC and **not** GCWGC. The same goes for the
remnant, it'll be C\[AT\]GC. The software will look for cut sites in a single read,
and remove anything other than the first valid fragment. If there are no cut sites, the
reads are considered valid pairs and output separately. If there are cutsites, the read
is a chimera and it's written as such, with the two ends written as single files because
they really aren't paired, just two random fragments that ended up joined together and
sequenced in one go.

If the reads are missing the initial remnant, they're considered invalid and discarded
completely. If you want to skip this first remnant check completely, just enter a . as the
remnant and all reads will be valid (have a valid "any character, once" at the start of it).
