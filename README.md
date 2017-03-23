# BIOINF545 ATAC-seq Lab

* [Introduction](#intro)
* [Setting up your workspace](#setup)
* [Retrieval of sequence data](#retrieval)
* [Basic quality checking with FastQC](#fastqc)
* [Trimming adapter sequence from reads](#trimming)
* [Aligning the trimmed reads to a reference genome](#aligning)
* [Sifting the aligned reads](#sifting)
* [Calling peaks on the aligned reads](#callingpeaks)
* [Creating a browser track so we can look at the peaks in the UCSC Genome Browser](#browsertrack)
* [Predicting transcription factor binding footprints](#footprinting)

## <a name="intro"></a>Introduction

We're going to work through a basic ATAC-seq data analysis
pipeline. We'll check the quality of the assay data, clean it up, map
it to a reference genome, call peaks on the aligned reads, and create
a browser track for the UCSC Genome Browser. We'll also predict the
locations of bound transcription factors.

The commands you'll run will be set in code blocks with a gray
background, like this:

```bash
echo "This is a command."
echo "This is another command."
```

To run each command, copy the entire line and paste it
into a shell window on the class workstation.

## <a name="setup"></a>Setting up your workspace

Open up a terminal session and create a working directory under your
home directory by copying and pasting the commands below:

```bash
export LAB_DIR=~/bioinf545/labs/atac-seq
mkdir -p $LAB_DIR && cd $LAB_DIR
```

Set a few environment variables to reduce typing later:

```bash
export LAB_DATA="/class/data/bio545/atac-seq-lab"
export REF_DIR=${LAB_DATA}/data
export REF=hg19
export PICARD_JAR="${LAB_DATA}/bin/picard.jar"
export R_LIBS_SITE=${LAB_DATA}/R/%p/%v
```

Make sure we're running the right versions of the tools:

```bash
export PATH=${LAB_DATA}/bin:${LAB_DATA}/ve/bin:$PATH
```

## Retrieval of sequence data

Whether trying to replicate published results or working with your own
data, it's essential to know how high-throughput sequence data was
produced. Paired-end data is processed differently than single-ended,
and many of the tools we use are sensitive to details like read length
and reference genome size.

We're going to work with data from the original ATAC-seq paper by
Buenrostro et al.:

http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959825/

We don't have time to download and analyze the entire sample during lab, so we're
going to start with a subset in the following FASTQ files:

* SRR891268.1.fq.gz -- contains the first reads of the pairs
* SRR891268.2.fq.gz -- contains the second reads of the pairs

They contain reads from sample SRR891268 that aligned (mostly) to
chromosome 1.

Copy them from the class data directory to your working directory:

```bash
cp ${LAB_DATA}/data/*.gz ${LAB_DIR}
cp ${LAB_DATA}/data/*.meme ${LAB_DIR}
```

Now list your directory contents to make sure the files are there:

```bash
ls -lFh
```

You should see files that look like this:

```bash
# NOTE: do not copy/paste this block. This is only an example of output.
total 21M
-r--r--r--. 1 scjp users  19M Mar 23 11:50 chr20.fa.gz
-r--r--r--. 1 scjp users  684 Mar 23 11:50 CTCF_known2.meme
-r--r--r--. 1 scjp users 730K Mar 23 11:50 NA12878.CTCF_known2.chr1.bed.gz
-r--r--r--. 1 scjp users  63K Mar 23 11:50 NA12878.CTCF_known2.scores.chr1.gz
-r--r--r--. 1 scjp users 496K Mar 23 11:50 SRR891268.1.fq.gz
-r--r--r--. 1 scjp users 512K Mar 23 11:50 SRR891268.2.fq.gz
```

If you're curious about how to retrieve published sequence data, read
on, otherwise skip to [Basic quality checking with FastQC](#fastqc).

The authors submitted their data to the NCBI Sequence Read Archive. In
footnotes under the paper's methods section, they provide a link:

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753

Our data comes from the first sample:

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1155957

From that page, we know the sample was extracted from 50,000 human
lymphoblastoid cells from the GM12878 line. The library was sequenced
on an Illumina HiSeq 2000, producing paired-end reads.

The sample page also links to the SRA record:

http://www.ncbi.nlm.nih.gov/sra?term=SRX298000

The run ID at the bottom of the page, `SRR891268`, is what you need to
download the data. You also need to have the NCBI SRA tools installed,
so you can use fastq-dump:

You do not need to do this step for the lab.

```bash
fastq-dump --gzip --split-files SRR891268
```

The `--split-files` argument is used to download into two files, one
containing the first read of the pairs, and the other containing the
second reads (mates). Tools used later in the pipeline often use the presence of the
mate file to determine whether they're processing paired end
data.

## <a name="fastqc"></a>Basic quality checking with FastQC

FastQC performs some basic quality control checks on raw sequence
data. To check the ATAC-seq reads, run:

```bash
fastqc SRR891268.1.fq.gz
```

When it completes, copy the file to your maching and open the HTML report in a web browser.
You can execute a command _similar_ to this from your own Desktop directory:
Note to substitute USERNAME with your own user name.

```bash
scp USERNAME@bcs2.bioinformatics.med.umich.edu:/users/USERNAME/bioinf545/labs/atac-seq/SRR891268.1_fastqc.html .
```

Note the report section on per-base content. Do you think TN5 has integration bias?

Don't worry about the GC content. We're using too little data to get a reliable curve.

Note the section on adapter contamination. How does this look? Is it expected?

When you're done reviewing, exit the browser. In real work you'd of
course check both read 1 and 2 files, but in the interest of time we'll move on.

## <a name="trimming"></a>Trimming adapter sequence from reads

For ATAC-seq data, we trim adapter sequence using Jason Buenrostro's
approach: try to align the paired-end reads to each other, and if that
can be done with a Levenshtein edit distance of one or less, chop off
any sequence outside the alignment. The nice thing about this
technique is that you don't need to know which adapter sequences were
used. Other tools generally need to be told, or can sometimes guess,
using a list of known adapters. The command line:

```bash
trim_adapters SRR891268.1.fq.gz SRR891268.2.fq.gz
```

When it's done, compare the first few reads in the original and
trimmed files of first reads:

```bash
zdiff -u SRR891268.1.fq.gz SRR891268.1.trimmed.fastq.gz | less
```

You should see something like this, where the lines marked with `+`
and `-` in the left margin differ. The `-` lines are the original
versions of the reads or quality lines, and the `+` lines are
trimmed. (Note that there are also comment lines in the FASTQ files
that start with `+` -- ignore those.)

```diff
--- /dev/fd/5   2017-03-21 14:23:27.830504681 -0400
+++ -   2017-03-21 14:23:27.833919209 -0400
@@ -15,9 +15,9 @@
 +SRR891268.38259 HWI-ST281:266:C1LTTACXX:1:1101:14488:7554 length=50
 CCCFFFFFGHHHHJJJJJJJJJJJJIJ@GIHIHJJJBEEHBDFFEEDDDB
 @SRR891268.38633 HWI-ST281:266:C1LTTACXX:1:1101:18609:7598 length=50
-TTTCTCGTGTTACATCGCGCCATCATTGGTATATGGCTGTCTCTTATACA
+TTTCTCGTGTTACATCGCGCCATCATTGGTATATGG
 +SRR891268.38633 HWI-ST281:266:C1LTTACXX:1:1101:18609:7598 length=50
-CCCFFFFFHHHHHJJJJJJJJJJJJJJJJGHIJJJJJJJJJJJJJJIJJJ
+CCCFFFFFHHHHHJJJJJJJJJJJJJJJJGHIJJJJ
 @SRR891268.43221 HWI-ST281:266:C1LTTACXX:1:1101:12315:8330 length=50
 GGGCCGGGCGGTCCCTTTAACGGCGCGGCCCGAGGGGCGCAGGCGGGAGG
 +SRR891268.43221 HWI-ST281:266:C1LTTACXX:1:1101:12315:8330 length=50
```

You can exit from the less program by pressing `q`.

## <a name="aligning"></a>Aligning the trimmed reads to a reference genome

With the adapter cleanup complete, we can finally align the reads to a
reference genome and see where the ATAC-seq transpositions happened.
Note this is a long command. Make sure to copy the entire line.

```bash
bwa mem -t 4 ${REF_DIR}/${REF} SRR891268.1.trimmed.fastq.gz SRR891268.2.trimmed.fastq.gz | samtools sort -@ 4 -O bam -T SRR891268.tmp -o SRR891268.bam -
```

We specify bwa's `mem` algorithm, and give it both files of paired-end
reads. The `mem` algorithm is the latest, and recommended for any
reads longer than 70bp. It also requires just one step, which is why
we're using it with the 50bp reads in this lab, but if you're ever
working with short reads, you'll probably want to at least try the
older "backtrack" algorithm (invoked with `bwa aln` ) for each file of
reads, add the separate `bwa sampe` step to combine the results for
each pair, and compare to the `bwa mem` alignments.

We also pipe bwa's output through `samtools sort` to create the final
BAM file. You'll see a lot of piping in bioinformatics analyses on
Linux. It's generally more efficient, since each command doesn't have
to write its output to disk. Sometimes it is worth preserving the
output of big tasks, though, if you know you'll be feeding it to
multiple downstream processes.

The `-O bam` argument to `samtools sort` requests BAM output, the `-T`
option specifies the basis for the temporary files it creates while
sorting, `-o` names the output file, and `-` specifies that the input
will come from standard input -- the pipe into which bwa sends its
output. You may also see standard input specified as `/dev/stdin`,
usually when a program doesn't recognize `-`; `/dev/stdin` just looks
like a regular file to them.

Finally, note the `-t 4` option: we're telling bwa to use four
processing threads to align the reads faster. We also tell samtools to
use four threads for sorting and compressing with the `-@ 4` option.

On Linux, you can see how many processors are available with the
`lscpu` command. Picking the right number of threads can be
tricky. Too few and your analysis takes longer than it should, but too
many and it could take even longer, as the machine struggles to
balance all the work. You need to know how busy the machine is, and
also how well a program can use multiple processors; some don't scale
well, so there's a point of diminishing returns, after which you're
wasting processors and not getting your results any sooner.

The combination of bwa and samtools is pretty efficient. Running the
above command took about 24 seconds with one thread, and only six
seconds with four. That's on this tiny sample data; with typical
genomic analyses, the difference can be many hours. But running with
eight threads only shaved another second and a half from the run time,
and even with 16 threads, it still took four seconds.

When you have a lot of data to align, it can be more efficient to run
multiple bwa commands concurrently, with a few threads each, than to
run one at a time with a large number of threads.

## <a name="sifting"></a>Sifting the aligned reads

Not all reads map well. We'll use bwa's annotations to sift out the
good ones for subsequent analysis.

Then we'll mark duplicate alignments. Duplication is a complicated
assessment: duplicate reads can come from the same original DNA
fragment, or they can be PCR or optical artifacts of the library prep
or sequencing process. The documentation for the tool we'll use,
[Picard, from the Broad Institute](https://broadinstitute.github.io/picard/command-line-overview.html),
explains how it identifies duplicates.

Here's the command:

```bash
java -Xmx8g -jar ${PICARD_JAR} MarkDuplicates I=SRR891268.bam O=SRR891268.md.bam ASSUME_SORTED=true METRICS_FILE=SRR891268.markdup.metrics VALIDATION_STRINGENCY=LENIENT
```

The output will be a BAM file containing all of bwa's output, with
duplicate reads marked. You could also just have Picard remove them, if you
have no need for them later in the pipeline.

Now we need to index the BAM file with duplicates marked:

```bash
samtools index SRR891268.md.bam
```

Finally, we'll sift out the good alignments -- reads that mapped
uniquely, with good quality, to autosomal references:

```bash
export CHROMOSOMES=$(samtools view -H SRR891268.md.bam | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e chrX -e chrY -e 'VN:' | sed 's/SN://' | xargs echo); samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 SRR891268.md.bam $CHROMOSOMES > SRR891268.pruned.bam
```

Yes, really. You'll see complicated commands strung together like this
all the time. If it helps, this is more complex than average.

The first bit, before the semicolon, creates an environment variable
`CHROMOSOMES` to hold a list of autosomal references obtained from the
header of the BAM file:

> export CHROMOSOMES=$(samtools view -H SRR891268.md.bam | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e chrX -e chrY -e 'VN:' | sed 's/SN://' | xargs echo);

That environment variable is used in the last argument to the
`samtools view command to only retrieve reads that aligned to those
references:

> samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 SRR891268.md.bam $CHROMOSOMES

As for the rest of the arguments:

- `-b`: requests BAM output
- `-h`: requests that the header from the input BAM file be included

The next few use SAM flags to filter alignments. There's a
[detailed specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for SAM files, which describes the flags that
tools can use to annotate aligned reads.

- `-f 3`: only include alignments marked with the SAM flag `3`, which means
  "properly paired and mapped"
- `-F 4`: exclude aligned reads with flag `4`: the read itself did not map
- `-F 8`: exclude aligned reads with flag `8`: their mates did not map
- `-F 256`: exclude alignments with flag `256`, which means
  that bwa mapped the read to multiple places in the reference genome,
  and this alignment is not the best
- `-F 1024`: exclude alignments marked with SAM flag `1024`, which
  indicates that the read is an optical or PCR duplicate (this flag
  would be set by Picard)
- `-F 2048`: exclude alignments marked with SAM flag `2048`,
  indicating chimeric alignments, where bwa decided that parts of the
  read mapped to different regions in the genome. These records are
  the individual aligned segments of the read. They usually indicate
  structural variation. We're not going to base peak calls on them.

Finally, we use a basic quality filter, `-q 30`, to request
high mapping-quality alignments.

## <a name="callingpeaks"></a>Calling peaks on the aligned reads

We'll use [MACS2](https://github.com/taoliu/MACS) to "call peaks" in the aligned reads -- we're looking
for regions with lots of transposition events, which indicate open
chromatin.

```bash
macs2 callpeak -t SRR891268.pruned.bam -n SRR891268.broad -g hs -q 0.05 --nomodel --shift -100 --extsize 200 -B --broad
```

The arguments are:

- `-t`: the "treatment" file -- the input, which is the sifted BAM
  file from the last step
- `-n`: the name of the experiment, which is used to name files
- `-g`: the genome's mappable size; 'hs' is an alias for the human
  genome's mappable size
- `-q`:  the false discovery rate cutoff for significant regions
  (peaks)
- `--nomodel`, `--shift`, and `--extsize`: MACS2 was designed for
  ChIP-seq data, so we're telling it not to use its built-in model,
  but to extend and shift reads in a way appropriate for ATAC-seq.
- `-B`: Create bedGraph files we'll use to create a browser track.
- `--broad`: request that adjacent enriched regions be combined into
  broad regions

## <a name="browsertrack"></a>Creating a browser track so we can look at the peaks in the UCSC Genome Browser

We're going to use a prefab browser track that contains peaks called
on the entire data, not the subset we're working with in the lab, but
this is how you would create a track for the called peaks:

```bash
LC_COLLATE=C sort -k1,1 -k2,2n SRR891268.broad_treat_pileup.bdg > SRR891268.broad_treat_pileup.sorted.bdg
bedGraphToBigWig SRR891268.broad_treat_pileup.sorted.bdg ${REF_DIR}/${REF}.chrom_sizes SRR891268.broad_peaks.bw
```

First open a web browser and navigate to the following custom browser URL:
https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=Scjparker&hgS_otherUserSessionName=bf545%2DATAC%2Dseq

This should open to the GCK locus. GCK is an islet-specific gene that is not "on" in GM12878 cells.

Now scroll down and click on "add custom tracks" and then in the "Paste URL or data" box, paste the following track:

```
track type=bigWig name="GM12878 ATAC-seq peaks" description="GM12878 ATAC-seq peaks" visibility=full color=255,128,0 alwaysZero=on maxHeightPixels=50:50:50 windowingFunction=mean smoothingWindow=3 bigDataUrl=https://theparkerlab.med.umich.edu/gb/tracks/bioinf545/gm12878.broad_treat_pileup.bw
```

Then click on "submit", then click "go" to return to the GCK locus now with GM12878 ATAC-seq data. How does the GM12878 chromatin accessibility look at GCK? At the two flanking genes, POLD2 and YKT6? Can you see the open chromatin promoter regions?


## <a name="footprinting">Predicting transcription factor binding footprints</a>

One of the questions that ATAC-seq helps answer is, "Where might
transcription factor binding be happening in interesting cell types?" If we can
correlate regions of open chromatin revealed by ATAC-seq with
predicted transcription factor (TF) binding sites, we can identify
regions that are potentially bound in a cell.

Transcription factors usually bind to specific sequence motifs, which
can be described with a position weight matrix (PWM). We can use these
PWMs to search the genome for all of a transcription
factor's potential binding sites. For this, we use a program called
FIMO from the MEME Suite (http://meme-suite.org/doc/fimo.html). FIMO
builds a file of motif locations scored according to their PWM match.

As input, FIMO needs a FASTA file containing the reference in which
you want to find motifs, a MEME-format file describing motifs, and
optionally a background file containing the background nucleotide
frequencies in the reference. In the interest of time, we're going to
use chromosome 20 to illustrate FIMO usage.

To create the background file, run `fasta-get-markov`:

```bash
zcat chr20.fa.gz | fasta-get-markov /dev/stdin chr20.bg
```

To scan the CTCF motif across chr20 and find sequence matches:

```bash
zcat chr20.fa.gz | fimo -bgfile chr20.bg CTCF_known2.meme /dev/stdin
```

FIMO's output goes into a `fimo_out` dir. The motifs are in a GFF
file, which you can convert to BED format with:

```bash
gff2bed < fimo_out/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' | gzip > NA12878.CTCF_known2.chr20.bed.gz
```

For TF footprinting, we use a program called CENTIPEDE
(http://centipede.uchicago.edu/), which combines binding site
coordinates and experimental data showing open chromatin to calculate
the probability that a transcription factor is present at a candidate
binding site.

To present the experimental data to CENTIPEDE, we've written a script
called `make_cut_matrix` to count ATAC-seq transposition events around
putative binding sites. Its input is an indexed BAM file (we'll have
to make sure to index the pruned BAM file we just created) and the BED
file of binding site motifs found by FIMO. Its output is a matrix
which can be passed to CENTIPEDE along with the BED file of binding
site motifs found by FIMO. Here's how to run it, using FIMO results
from chromosome 1:

```bash
samtools index SRR891268.pruned.bam
make_cut_matrix -v -d -p 2 -f 3 -F 4 -F 8 -q 30 SRR891268.pruned.bam NA12878.CTCF_known2.chr1.bed.gz | gzip -c > NA12878.CTCF_known2.matrix.gz
```

Finally, we've written an R script to invoke CENTIPEDE called,
predictably enough, `run_centipede.R`. Here's the command to predict
bound TFs in the lab data:

```bash
run_centipede.R NA12878.CTCF_known2.matrix.gz NA12878.CTCF_known2.chr1.bed.gz NA12878.CTCF_known2.centipede.bed.gz 8
```

The arguments are the matrix, the list of motif locations, the name of
the file into which CENTIPEDE should write its predictions, and the
field in the motif input file that contains the motifs' scores.

You can create a genome browser track from the CENTIPEDE output, with
a little cleanup:

```bash
zcat NA12878.CTCF_known2.centipede.bed.gz | awk 'BEGIN{IFS="\t"; OFS="\t"} $NF > 0.99 {print $1,$2,$3,$4,$5,$6,$NF}' | gzip > NA12878.CTCF_known2.bound.bed.gz
```

We have again created a browser track with predicted CTCF binding
sites in the entire hg19 reference, not just chromosome 1. If you add
this track to your earlier Genome Browser session, you should see
several sites around the GCK locus where CTCF is predicted to
be bound. Just paste this into the custom track submission form:

```
browser position chr7:44,116,289-44,266,468
track name="Predicted bound CTCF motifs in GM12878" description="Predicted bound CTCF motifs in GM12878" visibility=full color=0,128,64 alwaysZero=on maxHeightPixels=50:50:50 windowingFunction=mean smoothingWindow=3
https://theparkerlab.med.umich.edu/gb/tracks/bioinf545/SRR891268.CTCF_known2.bound.bed
```

To see all of this put together in the entire human reference genome,
explore this Genome Browser session:

https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=Scjparker&hgS_otherUserSessionName=bf545%2DATAC%2Dseq%2Dall

It includes tracks for CTCF binding sites found by FIMO, bound CTCF
motifs predicted by CENTIPEDE with GM12878 ATAC-seq data, ATAC-seq
peaks called on GM12878 ATAC-seq data, and ChIP-seq signal for CTCF
in GM12878.
