#!/usr/local/bin/python3.4

import os
import sys
import argparse

def main():

    # parse command line
    parser = argparse.ArgumentParser()
    # First define all option groups
    group1 = parser.add_argument_group('Input files', 'Required input')
    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")
    args = parser.parse_args()
    # manipulate path to import functions
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(1, parent_dir)
    mod = __import__('minvar')
    sys.modules["minvar"] = mod
    import logging
    import logging.handlers
    # Make a global logging object.
    mylog = logging.getLogger(__name__)
    dn = os.path.dirname(__file__)
    # set logging level
    mylog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './minvar.log'
    h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                             maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s%(lineno)d %(message)s")
    h.setFormatter(f)
    mylog.addHandler(h)
    mylog.info(' '.join(sys.argv))

    from minvar import prepare
    prepare.main(args.f)

    from minvar import callvar
    recal_file = callvar.main(ref_file='cns_final.fasta',
    	                      bamfile='hq_2_cons_sorted.bam', caller='lofreq')

    from minvar import annotate
    annotate.main(recal_file)

    from minvar import reportdrm
    reportdrm.main()

    from minvar import stats
    stats.coverage_stats('hq_2_cons_sorted_recal.bam')

if __name__ == "__main__": #  and __package__ is None:
    main()
