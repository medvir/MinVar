## Description of files with reference sequences


- `HCV_ref.fasta` (15 sequences) contains reference sequences for HCV as in the
  list maintained by ICTV; the file is created with `../../scripts/fetch_hcv_ref.py`

- `HCV_single_cons.fasta` is the consensus created from the above file with
  muscle and cons (see command below)

- `HIV1_CON_2002_genome_DNA[_aligned].fasta` (20 sequences) was downloaded from
  HIV LANL and contains the Consensus/Ancestral sequences, full genome, from all
  subtypes including recombinant form (alignment ID: 102CG1), aligned and raw

- `HIV1_CON_pol_DNA.fasta` (24 sequences) pol region of consensus sequences
  from all subtypes

The HCV consensus was obtained with the commands

```
muscle -in HCV_ref.fasta -out msa.fasta
cons -sequence msa.fasta -outseq c.fasta -plurality 2 -setcase 1 -name HCV_CONSENSUS
```
