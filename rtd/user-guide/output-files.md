## Description of output

The most important outputs are `annotated_DRM.csv`, a comma-separated-value file
that can be conveniently opened in any spreadsheet program, and `report.md`,
a mostly comprehensive report of the sample written in markdown format and that
includes the columns from `annotated_DRM.csv`. Other files are discussed below.
The columns are, in order,

- `gene`, can be RT (reverse transcriptase), protease, or integrase (GagPolTF
  is included in `annotated_DRM.csv`, but not in the final report),
- `pos`, a number indicating the position on that gene where a
  mutation with respect to
  [consensus B](https://hivdb.stanford.edu/pages/documentPage/consensus_amino_acid_sequences.html)
  sequence was observed,
- `mut`, the amino acid found at that position, if different from consensus B,
- `freq`, the frequency at which it was found,
- `category`, the classification of the mutation according to HIVDB.

The report also includes a rough estimate of the subtype for the analyzed sample.
This is done by aligning a subset of reads to sequences representative of
different subtypes and choosing the best match. The distribution of best matches
gives an idea of the most likely subtype for the sample. A better approach is
to take the sample consensus and run HIV BLAST (see below).

### Other output

#### Sample consensus

A consensus sequence for the sample is found in `cns_final.fasta`. This can be
used, for example, to run
[HIV BLAST](https://www.hiv.lanl.gov/content/sequence/BASIC_BLAST/basic_blast.html)
from Los Alamos HIV Database.

#### Main alignment file

The alignment of (at most) 200,000 reads to the sample consensus is
in file `hq_2_cons_sorted.bam`. One can explore this alignment, for example,
with the command

    samtools tview hq_2_cons_sorted.bam cns_final.fasta
