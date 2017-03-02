# MinVar: automatic detection of drug-resistance mutations in HIV-1

MinVar is a command-line tool to discover mutations conferring drug resistance
in HIV-1 populations using deep sequencing data.

----

#### The simplest example

    [user@host ~]$ minvar -f sample_file.fastq
    ... a few minutes later ...
    [user@host ~]$ column -t -s ',' annotated_DRM.csv
    gene      pos  mut  freq    category
    ...
    RT        238  T    1.0     NNRTI
    RT        250  N    0.9547  unannotated
    RT        272  P    1.0     unannotated
    RT        293  V    1.0     unannotated
    RT        297  A    1.0     unannotated
    RT        333  D    0.9384  unannotated
    RT        333  E    0.0354  unannotated
    RT        335  C    1.0     unannotated
    protease  10   P    0.0223  Other
    protease  10   Q    0.0185  Other
    protease  10   S    0.0741  Other
    protease  10   T    0.0468  Other
    protease  10   V    0.5948  PIMinor
    protease  11   L    1.0     PIMinor
    protease  13   V    1.0     unannotated
    protease  14   R    1.0     unannotated
    protease  15   V    0.7143  unannotated
    protease  20   T    1.0     PIMinor
    protease  32   I    1.0     PIMajor
    ...

## Important features

- MinVar is an _opinionated_ software: it just takes a fastq file as input and
  does not ask the standard user to set any parameter at run time. Nevertheless,
  the experienced user/developer can easily change some of its settings in the
  source code.
- It has been tested on both Illumina MiSeq and Roche/454 sequencing reads.
- It uses state-of-the-art third tools to filter, recalibrate, and align reads
  and to call variants.
- Finally, single nucleotide variants are phased at codon level and
  amino acid mutations are called and annotated.
- Drug-resistance mutations are annotated according to
  [Stanford HIV Drug Resistance Database (HIVDB)](https://hivdb.stanford.edu).
- The annotated mutations are saved in a csv file (see example above) and also
  included in a report in markdown format that is finally converted to PDF.

#### Citation

MinVar has been introduced and validated in  
Huber, Metzner _et al._, (2017) MinVar: A rapid and versatile tool for HIV-1
drug resistance genotyping by deep sequencing _Journal of virological methods_
240:7-13, [doi:10.1016/j.jviromet.2016.11.008](http://dx.doi.org/10.1016/j.jviromet.2016.11.008)
