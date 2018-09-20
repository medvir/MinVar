## MinVar: automatic detection of drug-resistance mutations in HIV-1

MinVar is a command-line tool to discover mutations conferring drug resistance
in HIV-1 and HCV populations using deep sequencing data.

----

#### The simplest example

    [user@host ~]$ minvar -f sample_file.fastq
    ... a few minutes later ...
    [user@host ~]$ column -t -s ',' merged_muts_drm_annotated.csv
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
- It has been tested with HIV-1 on both Illumina MiSeq and Roche/454 sequencing
  reads. HCV has been tested on MiSeq only.
- It uses state-of-the-art third tools to filter, recalibrate, and align reads
  and to call variants.
- Finally, single nucleotide variants are phased at codon level and
  amino acid mutations are called and annotated.
- HIV-1 drug-resistance mutations are annotated according to
  [Stanford HIV Drug Resistance Database (HIVDB)](https://hivdb.stanford.edu).
- The annotated mutations are saved in a csv file (see example above) and also
  included in a report in markdown format that is finally converted to PDF.
- The PDF report can be customized by adding contact information specified in the
  file `~/.minvar/contac.ini` with the following syntax (only change what comes
  after the `=` sign)

```
[contact]
unit = name_of_your_unit_here
phone = phone_number
fax = fax_number
email = your_unit@your_company
logo = filename_without_extension
```

The logo file in pdf format must be present in the same directory. In other words,
if we want to use the file `~/.minvar/company_logo_bw.pdf`, then in the INI file we
will write `logo = company_logo_bw`.

#### Citation

MinVar (version 1, HIV-1 support only) has been introduced and validated in  
Huber, Metzner _et al._, (2017) MinVar: A rapid and versatile tool for HIV-1
drug resistance genotyping by deep sequencing _Journal of virological methods_
240:7-13, [doi:10.1016/j.jviromet.2016.11.008](http://dx.doi.org/10.1016/j.jviromet.2016.11.008)


## Output files

### Created by `prepare.py`

- `subtype_evidence.csv` percent of reads best aligned to each subtype (or
  genotype),
- `subtype_ref.fasta` references of the subtype identified,
- `cns_final.fasta`: sample consensus created by iteratively aligning reads and
  writing variants into the sequence,

### Created by `callvar.py`

- `hq_2_cns_final_recal.bam` sorted bam alignment of reads to the consensus
  sequence, recalibrated with either GATK or lofreq (indels only),
- `hq_2_cns_final_recal.vcf` VCF file of mutations found on reads with respect
  to consensus in `cns_final.fasta`.

### Created by `annotate.py`

- `merged_mutations_nt.csv` a list of all variants observed at single positions,
- `max_freq_muts_aa.csv` the amminoacid found at maximum frequency at each codon,
- `final.csv` mutations at amminoacid level with indication of
  the gene, the position on the gene, wild type and frequency

### Created by `reportdrm.py`

- `merged_muts_drm_annotated.csv` is the join of `final.csv` with the annotation
  of DRM/RAS,
- `report.md` and `report.pdf` final report with subtye estimate based on
  alignment of reads to different references and tables with mutations. The
  pdf report is created from the template `minvar/db/template.tex`.


## Add a new reference to the sequence database

MinVar looks for reference sequences in two files. Respectively, in

- `src/minvar/db/organism/subtype_references.fasta` (for non-recombinant forms) and
- `src/minvar/db/organism/recomb_references.fasta` (for recombinant forms),

where `organism` is HIV or HCV.

New reference sequences can be added there, provided that related data structures in
`src/minvar/common.py` are updated to reflect the reference names as outlined below.

#### If an HIV reference sequence was added

- Add a key:value pair to `hiv_map` where the key is the id as in the fasta header and the
  value is an abbreviation for it. Example `'CONSENSUS_12_BF':'CRF12_BF'`.
- Add a key:value pair to `org_dict` where the key is the id as in the fasta header and the
  value is `HIV`.

#### If an HCV reference sequence was added
- Edit (if the reference added is from a genotype already present) or add (if a new genotype is being added)
  the key:value pair in `acc_numbers` where the key is the genotype of the sequence being added (_e.g._ 1a)
  and the value is the list of sequence ids as in the fasta headers for that genotype.

#### Citation

MinVar has been introduced and validated in  
Huber, Metzner _et al._, (2017) MinVar: A rapid and versatile tool for HIV-1
drug resistance genotyping by deep sequencing _Journal of virological methods_
240:7-13, [doi:10.1016/j.jviromet.2016.11.008](http://dx.doi.org/10.1016/j.jviromet.2016.11.008)
