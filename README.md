<h1 align="center">pb-assembly</h1>
<p align="center">PacBio Assembly Tool Suite:
    Subreads in ⇨ Assembly out
</p>

***

# Availability

The latest pre-release, unstable, experts-only linux/mac binaries can be installed via [bioconda](https://bioconda.github.io/).

```bash
conda install pb-assembly
```

Alternatively, if you don't have administrator access you can install the pb-assembly suite into an environment 
in your $HOME directory

```bash
conda create -n denovo_asm
source activate denovo_asm
conda install pb-assembly
```
    
If you are looking for a GUI based long read denovo genome assembler, you are urged to learn about 
[HGAP4](https://www.pacb.com/videos/tutorial-hgap4-de-novo-assembly-application/)
    
These binaries are not ISO compliant.
For research only.
Not for use in diagnostics procedures.

No support for source builds.
No support via mail to developers.
Please *do not* contact a PacBio Field Applications Scientist or PacBio Customer Service for assistance.
Please file GitHub issues for problems and questions!

**This is an early beta!** Expect extreme changes and different output between
versions until release of the first stable release.
Furthermore, the command-line options are not stable yet,
and can change at any point, do not rely on it yet.

# Scope

_pb-assembly_ is the bioconda recipe encompassing all code and dependencies necessary to
run the FALCON assembly pipeline and subsequently attempt to phase your genome with FALCON_unzip
and finally polish with arrow.

Installed package recipes include:

    - pb-falcon
    - pb-dazzler
    - genomicconsensus
    - etc (all other dependencies)
    
# Overview

FALCON and FALCON-Unzip are de novo genome assemblers for PacBio long reads, also known as 
single-molecule real-time (SMRT) sequences. FALCON is a diploid-aware assembler which follows 
the hierarchical genome assembly process (HGAP) and is optimized for large genome assembly (e.g. 
non-microbial). FALCON produces a set of primary contigs (a-contigs), which represent divergent 
allelic variants. Each a-contig is associated with a homologous genomic region on an p-contig.

FALCON-Unzip is a true diploid assembler. It takes the contigs from FALCON and phases the reads 
based on heterozygous SNPs identified in the initial assembly. It then produces a set of partially
phased primary contigs and fully phased haplotigs which represent divergent haplotyes.

# Usage

## Assemble

```bash
fc_run fc_run.cfg
```

## Unzip and polish
```bash   
fc_unzip fc_unzip.cfg
```

## Configuration

Both FALCON and FALCON_unzip take a config file as their only input parameter.

### FALCON Configuration

The FALCON pipeline is complex and has a multitude of configuration options. Many of the tools that comprise 
the FALCON Assembly pipeline were written by Gene Meyers and are extensively documented at his 
[dazzlerblog](http://dazzlerblog.wordpress.com). The FALCON denovo assembly pipeline consists of 3 basic stages,
raw read overlapping and consensus calling, also known as **pre-assembly**, pre-assembled read overlapping
or **pread overlapping**, and finally the assembly itself.

Here is a sample [fc_run.cfg](cfgs/fc_run_200kb.cfg) that was designed to work with the 200kb test case found
below

Here is a sample [fc_run.cfg](cfgs/fc_run_human.cfg) that was used with a recent ~2.9Gb Human assembly.
 
Below is a breakdown of the configuration options available to FALCON

#### Input


```ini
[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=pass
```

Your list of paths to the input fasta files is specified in your `input_fofn` and your `input_type` can be
either `raw` or `pread`. If specifying `pread`, the pipeline will skip the entire `0-rawreads` pre-assembly phase.

By default, dusting is turned on and is run after generating the raw read database with default options as 
recommended by Gene Meyer's. If you wish to modify your 
[dusting parameters](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/) you can set the 
flag `pa_DBdust_option`.

Filtering options for your input data for pre-assembly can also be set with the `pa_fasta_filter_option` flag.
Recognized values are described below

|Value | Setting          
|:-----|----------------
|pass  |The no-op filter - passes every FASTA record to stdout.                                                                           |
|median|Applies the median-length ZMW filter by running two passes over the data. Only one subread per ZMW is output, based on median-length selection.
|streamed-median|     Applies the median-length ZMW filter by running a single-pass over the data. The input subreads should be groupped by ZMW.
|internal-median|     Applies the median-length ZMW filter only on internal subreads (ZMWs with >= 3 subreads) by running two passes over the data. For ZMWs with < 3 subreads, the maximum-length one is selected.
|streamed-internal-median| Applies the median-length ZMW filter only on internal subreads (ZMWs with >= 3 subreads) by running a single pass over the data. The input subreads should be groupped by ZMW. For ZMWs with < 3 subreads, the maximum-length one is selected.

#### Data Partitioning

```ini
pa_DBsplit_option=-x500 -s200
ovlp_DBsplit_option=-x500 -s200
```

For the first and second stages of FALCON, the data needs to be read in to a 
[dazzler DB](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/). The `-x` flag filters 
reads smaller than what's specified while the `-s` flag controls the size of DB blocks. 


#### Repeat Masking

```ini
pa_HPCTANmask_option=
pa_REPmask_code=0,300;0,300;0,300
```

Repeat masking occurs in two phases, **Tandem** and **Interspersed**. Tandem repeat masking is run
with a modified version of `daligner` called `datander` and thus uses a similar 
[parameter set](https://dazzlerblog.wordpress.com/command-guides/damasker-commands/). Whatever settings you use
 for pre-assembly daligner overlapping in the next section will be used here for tandem repeat masking. You can 
 supply additional arguments for tandem repeat masking that will be passed to `HPC.TANmask` with the 
 `pa_HPCTANmask_option`.
 
The second phase of masking dealing with interspersed repeats can be run in up to 3 iterations specified with the
`pa_REPmask_code` option. The parameters needed for each iteration are both the group size and coverage specified 
as `group,coverage` pairs separated by semicolons as seen above. 

For information and theory on how to set up your rounds of repeat masking, consult this 
[blog post](https://dazzlerblog.wordpress.com/2016/04/01/detecting-and-soft-masking-repeats/).


#### Pre-assembly

```ini
genome_size=0
seed_coverage=20
length_cutoff=1000    
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option=-e.8 -l2000 -k18 -h480  -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 1800
falcon_sense_greedy=False
```

If you wish to auto-calculate your seed read coverage, then it's necessary to enter your `genome_size` in base
pairs, the desired `seed_coverage` as well as set `length_cutoff=-1` to force the auto-calculation. Alternatively, 
if you don't know your genome size, are unsure of the `seed_coverage` you would like to use or if you would rather 
just leverage all reads above a specific length, you can use the the `length_cutoff` flag to
manually set that limit. It's important to note that whatever value `length_cutoff` gets set to is a limit that
carries through to the unzipping algorithm, and any reads smaller than that cutoff will not be used for phasing. If
you're just assembling, there is probably no harm in setting `length_cutoff` high. If you are unzipping however,
then you will be artificially limiting your phasing dataset and it's probably in your interest to have a 
lower `length_cutoff`
 
Overlap options for `daligner` are set with the `pa_HPCdaligner_option` and `pa_daligner_option` flags. Previous 
versions of FALCON had a single parameter, however currently this was broken into two flags, one that affects requested 
resources `pa_HPCdaligner_option` and one that affects the overlap search `pa_daligner_option`.
For `pa_HPCdaligner_option`, the `-v` parameter is passed to the `LAsort` and `LAmerge` programs while `-B` and `-M` 
parameters are passed to the `daligner` sub-commands. To understand the theory and how to configure `daligner` see 
[this blog post](https://dazzlerblog.wordpress.com/2014/07/10/dalign-fast-and-sensitive-detection-of-all-pairwise-local-alignments/)
and this [command reference guide](https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide/) 
 
You can configure basic pre-assembly consensus calling options with the `falcon_sense_option` flag. The `--output-multi`
flag is necessary for generating proper fasta headers so should not be removed unless you know what you are doing. The
parameters `--min-idt`, `--min-cov` and `--max-n-read` set the minimum alignment identity, minimum coverage necessary 
and max number of reads, respectively, for calling consensus with the `arrow` consensus algorthim.

By default, `-fo` are the parameters passed to `LA4Falcon`. The option `falcon_sense_greedy` changes this
 parameter set to `-fog` which essentially attempts to maintain relative information between reads that have been
 broken.


#### Pread overlapping

```ini
ovlp_daligner_option=-e.96 -s1000 -h60 -t32
ovlp_HPCdaligner_option=-v -M24 -l500
```

The second phase of overlapping of corrected reads occurs in a similar fashion to the overlapping performed in the
pre-assembly, however no repeat masking is performed as the repeats have already been masked, no consensus is
called and overlaps are simply identified to be fed into the final assembly. The parameter options work the same 
was as described above in the Pre-assembly section.

#### Final Assembly

```ini
overlap_filtering_setting=--max-diff 100 --max-cov 100 --min-cov 2
fc_ovlp_to_graph_option=
length_cutoff_pr=1000
```

The parameter `overlap_filter_setting` allows one to set criteria for filtering corrected read overlaps. 
`--max-diff` filters overlaps that have a coverage difference between the 5' and 3' ends larger than specified. 
`--max-cov` filters highly represented overlaps typically caused by contaminants or repeats and `--min-cov` allows
you to specify a minimum overlap coverage. Setting `--min-cov` too low will allow more overlaps to be detected
at the expense of additional chimeric / mis-assemblies.
    
`length_cutoff_pr` is the minimum length of pre-assembled *preads* used for the final assembly. Typically you
want to set this value to allow for approximately 15-30X coverage of corrected reads in the final assembly.

#### Miscellaneous configuration options

Additional configuration options that don't necessarily fit into one of the previous categories are described here.


```ini
target=assembly
skip_checks=False
LA4Falcon_preload=false
```

FALCON can be configured to stop after any of it's three stages with the `target` flag set to either
`overlapping`, `pre-assembly` or `assembly`. Each option will stop the pipeline at the end of it's corresponding
stage, `0-rawreads`, `1-preads_ovl` or `2-asm-falcon` respectively. The default is full `assembly` pipeline.

The flag `skip_checks` disables `.las` file checks with `LAcheck` which has been known to cause errors on certain
systems in the past.

The parameter `LA4Falcon_preload` passes the `-P` parameter to `LA4Falcon` which causes all the reads to be loaded
into memory. On slow filesystems this can make a huge difference due to the random-access; but it will
dramatically increase the memory requirement for the consensus stage.

#### Job Distribution

```ini
[job.defaults]
job_type=sge
pwatcher_type=blocking
JOB_QUEUE = default
MB = 32768
NPROC = 6
njobs = 32
submit = qsub -S /bin/bash -sync y -V  \
  -q ${JOB_QUEUE}     \
  -N ${JOB_NAME}      \
  -o "${JOB_STDOUT}"  \
  -e "${JOB_STDERR}"  \
  -pe smp ${NPROC}    \
  -l h_vmem=${MB}M    \
  "${JOB_SCRIPT}"

[job.step.da]
NPROC=4
MB=49152
njobs=240
```
 
Default job configuration options are specified in the `[job.defaults]` section of your config file. The first 
option you should set is for the `job_type`. Allowed values are `sge`, `pbs`, `torque`, `slurm`, `lsf` and `local`. 
If running on a cluster, you need to configure the `submit` string to work with your job scheduler. The `submit`
string in the sample above is a tested and working SGE submit string. If you are running in `local` mode on a 
single machine, the submit string should be something like  `submit=bash -C ${CMD} >| ${STDOUT_FILE} 2>| ${STDERR_FILE}`.


Next, you need to tell FALCON how to deal with pipeline flow control by setting your process watcher `pwatcher_type`. 
There two possible values you can set, either `blocking` or `fs_based`. `fs_based` is the default and relies on
the pipeline polling the file system periodically to determine whether a `sentinel` file has appeared that would 
signal the pipeline to continue. The other option is to use a `blocking` process watcher which can help with 
systems that have issues with filesystem latency. In this case, the end of the job is determined by the finishing 
of the system call, rather than by file system polling.

Next you will find your job distribution settings. You will find settings for your default job queue `JOB_QUEUE`, 
memory allocated per job `MB`, number of processors per job `NPROC` as well as number of concurrently running 
jobs `njobs`. You will also find the `submit` option which allows you to set up your job scheduler submission 
string. The example above works with a typical SGE grid. 


Each stage of the assembly pipeline can be given different default parameters with different `[job.step.*]` sections.
There are 6 optional stages you can configure by using different 2-3 letter codes. `da` and `la` refer to 
pre-assembly daligner jobs and pre-assembly LAshow/LAsort jobs respectively. `cns` refers 
to the pread consensus calling stage. `pda` and `pla` refer to the pread daligner overlapping and pread
LAshow/LAsort stages respectively while `asm` refers to the final assembly. If you omit a specific `[job.step.*]` 
section, the `[job.defaults]` will be applied. `[job.step.da]`, `[job.step.la]`, `[job.step.cns]`, `[job.step.pda]`, 
`[job.step.pla]`, and `[job.step.asm]` are the available sections.
 
### FALCON_unzip Configuration

```ini
[General]
max_n_open_files = 1000
[Unzip]
input_fofn=input.fofn
input_bam_fofn=input_bam.fofn
```

Relative to setting FALCON parameters, FALCON_unzip is a breeze as the majority of the options have to do
exclusively with job distribution. The first and only setting in the [General] section is for `max_n_open_files`.
During the read tracking stage the pipeline can be writing to many `.sam` files at the same time. This can cause
problems with certain networked filesystems, so the default is to set `max_n_open_files=300`. Feel free to raise
this number if file system latency is not an issue for you.

Similar to FALCON, the parameter `input_fofn` simply refers to the input file of fasta names. This setting should
be redundant with your fc_run.cfg. Finally, if you wish to polish your unzipped genome, you will need to also 
specify a list of your input bam files with `input_bam_fofn`.

Here is a sample [fc_unzip.cfg](file://cfgs/fc_unzip.cfg) that will need to be tuned to your compute environment.

#### Job Distribution

Configuration of your `[job.defaults]` section is identical to FALCON as described previously. The only difference
are the job specific settings specific to FALCON_unzip. Available sections are `[job.step.unzip_track_reads]`, 
`[job.step.unzip_blasr_aln]`, `[job.step.unzip.phasing]` and `[job.step.unzip.hasm]`

## Example test case

TODO

## FAQ

TODO


## Detailed Description
The hierarchical genome assembly process proceeds in two rounds. The first round of assembly involves
the selection of seed reads, or the longest reads in the dataset (user-defined length_cutoff). All 
shorter reads are aligned to the seed reads, in order to generate consensus sequences with high accuracy. 
We refer to these as pre-assembled reads but they can also be thought of as “error corrected” reads. 
During the pre-assembly process, seed reads may be split or trimmed at regions of low read coverage 
(user-defined min_cov for falcon_sense_option). The performance of the pre-assembly process is captured 
in the pre-assembly stats file.

In the next round of HGAP, the preads, are aligned to each other and assembled into genomic contigs.

<img width="600px" src="img/HGAP.png" alt="HGAP" />

For more complex genomes assembled with FALCON, “bubbles” in the contig-assembly graph that result 
from structural variation between haplotypes may be resolved as associate and primary contigs. 
The unzip process will extend haplotype phasing beyond “bubble” regions, increasing the amount of phased 
contig sequence. It is important to note that while individual haplotype blocks are phased, phasing does 
not extend between haplotigs. Thus, in part C) of the figure below, haplotig_1 and haplotig_2 may 
originate from different parental haplotypes. Additional information is needed to phase the haplotype 
blocks with each other.

<img width="600px" src="img/FALCON_pipeline.png" alt="FALCON pipeline" />

Associate contig IDs contain the name of their primary contig but the precise location of alignment must 
be determined with third party tools such as NUCmer. For example, in a FALCON assembly, 000123F-010-01 
is an associated contig to primary contig 000123F. In a FALCON-Unzip assembly, 000123F_001 is a haplotig 
of primary contig 000123F.

Below are examples of alignments between associate and primary contigs from FALCON, and haplotigs and 
primary contigs from FALCON-Unzip. Alignments were built with NUCmer and visualized with Assemblytics. 
Precise coordinates may be obtained with the show-coords utilty from MUMmer.

<img width="600px" src="img/dotplots.png" alt="Associate contigs VS Haplotigs" />


## Acknowledgements
Thanks to Jason Chin for the original concept and Chris Dunn/Ivan Sovic for their numerous improvements.
## Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
