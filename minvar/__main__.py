#!/usr/bin/env python3
"""Runs everything"""
files_to_remove = [
    'calls_1.vcf.gz', 'cnsref.amb', 'cnsref.ann', 'cnsref.bwt', 'hq_smp.fastq',
    'hq_smp.fasta', 'hq_2_cns_final.bam', 'hq_2_cns_final.bam.bai',
    'cnsref.pac', 'cnsref.sa', 'high_quality.fastq', 'loc_res.tsv',
    'outblast.xml', 'refcon_sorted.bam', 'refcon_sorted.bam.bai',
    'sample_hq.fasta', 'phased.csv']

def main():
    """What does the main do?"""
    import os
    import sys
    import argparse
    # from setuptools_scm import get_version
    from pkg_resources import get_distribution, DistributionNotFound, \
        resource_filename
    try:
        __version__ = get_distribution('minvar').version

    except DistributionNotFound:
       # package is not installed
        pass

    # parse command line
    parser = argparse.ArgumentParser()
    # First define all option groups
    group1 = parser.add_argument_group('Input files', 'Required input')
    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")
    group1.add_argument(
        "-r", "--recal", action="store_true",
        help="turn on recalibration with GATK <default: %(default)s>",
        default=False)
    group1.add_argument("-k", "--keep", action="store_true",
                        help="keep intermediate files <default: %(default)s>",
                        default=False)
    group1.add_argument('-v', '--version', action='version',
                        version=__version__)

    args = parser.parse_args()

    # exit so that log file is not written
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    import logging
    import logging.handlers
    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() \
%(lineno)d: \t%(message)s'
    logging.basicConfig(filename='minvar.log', level=logging.INFO,
                        format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from minvar import prepare
    cns_file, prepared_bam, org_found = prepare.main(args.f)

    from minvar import callvar
    called_file, called_bam = callvar.main(ref_file=cns_file,
                                           bamfile=prepared_bam,
                                           caller='lofreq',
                                           recalibrate=args.recal)

    from minvar import annotate
    annotate.main(vcf_file=called_file, ref_file=cns_file, bam_file=called_bam,
                  organism=org_found)

    from minvar import reportdrm
    reportdrm.main(org_found)

    if not args.keep:
        for f in files_to_remove:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    if org_found == 'HIV':
        from minvar import stats
        bed_file = resource_filename(__name__, 'db/HIV/consensus_B.bed')
        stats.coverage_stats_per_base(prepared_bam, bed_file)
        cum_cov = stats.gene_coverage(prepared_bam, bed_file)
        for k, v in cum_cov.items():
            print(k, v)

if __name__ == "__main__": #  and __package__ is None:
    main()
