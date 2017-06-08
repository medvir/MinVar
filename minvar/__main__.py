#!/usr/bin/env python3

def main():

    import sys
    import argparse
    from setuptools_scm import get_version

    from pkg_resources import get_distribution, DistributionNotFound
    try:
        __version__ = get_distribution('minvar').version
    except DistributionNotFound:
       # package is not installed
       pass

    #__version__ = get_version(root='..', relative_to=__file__)

    # parse command line
    parser = argparse.ArgumentParser()
    # First define all option groups
    group1 = parser.add_argument_group('Input files', 'Required input')
    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")
    group1.add_argument("-r", "--recal", action="store_true",
                        help="turn on recalibration with GATK")
    group1.add_argument('-v', '--version', action='version',
                         version=__version__)

    args = parser.parse_args()

    # exit so that log file is not written
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    import logging
    import logging.handlers

    logging.basicConfig(filename='minvar.log', level=logging.DEBUG,
                        format='%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))



    from minvar import prepare
    cns_file, prepared_bam = prepare.main(args.f)

    from minvar import callvar
    called_file, called_bam = callvar.main(ref_file=cns_file,
                                           bamfile=prepared_bam,
                                           caller='lofreq',
                                           recalibrate=args.recal)

    from minvar import annotate
    annotate.main(vcf_file=called_file, ref_file=cns_file, bam_file=called_bam)

    from minvar import reportdrm
    reportdrm.main()

    from minvar import stats
    stats.coverage_stats_per_base(prepared_bam)
    stats.coverage_above_threshold(prepared_bam)

if __name__ == "__main__": #  and __package__ is None:
    main()
