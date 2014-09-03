##OptimAssembly
Viruses are characterized by, among other things, much smaller genome with
respect to bacteria. As a consequence, the coverage obtained in massively
parallel sequencing tends to be very high. This is not always good for
genome assembly. In fact, too many reads can decrease the quality of the
assembly (and increase the time/memory requirements), due to the abundance of
sequencing errrors.

This simple tool filters the raw reads and takes a sample of them in order to
reach a coverage that should give optimal results in genome assembly.

Sampling reads multiple times increases the chances that different contigs
spanning the whole genome will be reconstructed. These are aligned and then
the consensus is generated.

### Dependencies

- Python;
- Biopython for the sequence manipulation;
- [`seqtk`](https://github.com/lh3/seqtk) for the trimming;
- `velvet` for the assembly;
- `needle` and `cons` from [EMBOSS](http://emboss.sourceforge.net), for alignment and consensus building;
- [`muscle`](http://drive5.com/muscle/) for MSA.


### Usage
	usage: optimassembly.py [-h] [-f FASTQ] [-r REFERENCE] [-l EXP_LENGTH]

	Optimise de novo assembly for short, viral genomes

	optional arguments:
	  -h, --help            show this help message and exit
	  -f FASTQ, --fastq FASTQ input file in fastq format <>
	  -r REFERENCE, --reference REFERENCE closest known genome reference
	  -l EXP_LENGTH, --length EXP_LENGTH expected length <10000>
