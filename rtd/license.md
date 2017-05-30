# License

MinVar is licensed for free for non-commercial purposes according to
[this agreement](agreement.md). Its commercial use requires a separately
executed written license agreement, please
[contact us](https://ozagordi.github.io/MinVar/#contact).

MinVar relies on the following libraries with their own specific licenses:

- pandas (library for data analysis) [simplified BSD licence](http://pandas.pydata.org/pandas-docs/stable/overview.html#license)
- Biopython (library for biological computation), distributed under [Biopython licence](https://github.com/biopython/biopython/blob/master/LICENSE.rst)

The following external tools are called by the program:

- BLAST (sequence alignment), code distributed in lesser GPL
- samtools (sequence manipulation), [MIT licence](https://github.com/samtools/samtools/blob/develop/LICENSE)
- GATK (base quality recalibration), commercial licence (see below)
- picard-tools (sequence manipulation), only used together with GATK, [MIT licence](https://github.com/broadinstitute/picard/blob/master/LICENSE.txt)
- bwa (sequence alignment), [GPL v3 licence](https://github.com/lh3/bwa/blob/master/COPYING)
- seqtk (sequence filtering), [MIT licence](https://github.com/lh3/seqtk/blob/master/LICENSE)
- lofreq (variant calling), [MIT licence](https://github.com/CSB5/lofreq/blob/master/LICENSE)

### Using recalibration in a commercial setting

It is worth noting that recalibration is done by GATK, which is itself
subject to a dual purpose [licensing](https://software.broadinstitute.org/gatk/download/licensing.php):
free for academic non-commercial research activities, subject to a fee for commercial use.

Starting from v1.2, MinVar by defaults does not recalibrate the base quality
scores (you can evaluate the difference in this
[Note on recalibration](user-guide/recal.md)). If you decide to use MinVar
with recalibration for commercial purposes, you will have to obtain a licence
for GATK from Broad Institute.
