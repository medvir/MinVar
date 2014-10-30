#!/usr/bin/env python
''' Provides some methods to work with needle (EMBOSS) alignments'''
__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__copyright__ = ""
__license__ = ""

import sys
import os

dna_code = {
    'A': set(['A']),
    'C': set(['C']),
    'G': set(['G']),
    'T': set(['T']),

    'R': set(['G', 'A']),
    'Y': set(['T', 'C']),
    'M': set(['A', 'C']),
    'K': set(['G', 'T']),
    'S': set(['G', 'C']),
    'W': set(['A', 'T']),

    'H': set(['A', 'C', 'T']),
    'B': set(['C', 'G', 'T']),
    'V': set(['A', 'C', 'G']),
    'D': set(['A', 'G', 'T']),
    'N': set(['A', 'C', 'G', 'T']),
    '-': set(['A', 'C', 'G', 'T'])
}
GO = 6.0
GE = 3.0

class AlignDict(dict):
    '''
    A set of alignment instances. Subclass dict
    in order to have a 2D dictionary.
    '''

    def __init__(self, set_name, gap_open, gap_extend):
        '''
        Initialization needs a name, gap open, gap extension
        '''
        self.name = set_name
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.default = {}

    def __getitem__(self, key):
        '''
            subclassing dictionary getitem
            '''
        if key not in self:
            self[key] = {}
        return dict.__getitem__(self, key)


class AlignInstance:

    '''
    Represents a single pairwise alignment
    '''
    def __init__(self, id_a, id_b, seq_a, seq_b, score,
                 descr_a=None, descr_b=None):
        '''
        Initialization needs id of both sequences, zipped sequence
        pair, and score
        '''
        self.id_a = id_a
        self.id_b = id_b

        self.score = score
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.descr_a = descr_a
        self.descr_b = descr_b

    def summary(self):
        '''
        summary must be called to calculate start, stop,
        internal gaps, and identity
        '''
        import itertools
        import warnings

        start = None
        stop = None
        i = 0

        it_pair = zip(self.seq_a, self.seq_b)
        while True:
            i += 1
            try:
                p = next(it_pair)
            except StopIteration:
                break
            if p is None:
                break
            if start is None:
                if '-' not in p:
                    start = i
            else:
                if '-' not in p:
                    stop = i

        self.start = start
        self.stop = stop
        if start == None and stop == None:
            warnings.warn('The two sequences do not align')
            return
        self.insertions = 0  # gap in the first sequence
        self.deletions = 0  # gap in the second sequence
        self.ident = 0
        self.mismatches = 0

        it_pair = zip(self.seq_a[start - 1:stop],
                                 self.seq_b[start - 1:stop])

        # determine which sequence type to count mismatches correctly
        both = self.seq_a + self.seq_b
        both = both.replace('-', '').replace('*', '')
        nt_count = sum([both.upper().count(nt) for nt in ['A', 'C', 'G', 'T']])
        if nt_count > 0.75 * len(both):
            seq_type = 'dna'
        else:
            seq_type = 'aa'

        while True:
            i += 1
            try:
                p = next(it_pair)
            except StopIteration:
                break
            if p is None:
                break
            if p == ('-', '-'):
                print(' double gap in %s %s' % \
                    (self.id_a, self.id_b), file=sys.stderr)
            if p[0] == '-':
                self.insertions += 1
            elif p[1] == '-':
                self.deletions += 1
            if p[0].upper() == p[1].upper():
                self.ident += 1
            # count mismatches with appropriate method for dna and aa
            if seq_type == 'dna' and \
               (dna_code[p[0].upper()] & dna_code[p[1].upper()] == set([])):
                self.mismatches += 1
            elif seq_type == 'aa' and (p[0].upper() != p[1].upper()):
                self.mismatches += 1

        return

    def print_info(self):
        '''
            Print all information to stderr
            '''

        print('id_a:', self.id_a, file=sys.stderr)
        print('id_b:', self.id_b, file=sys.stderr)
        print('start:', self.start, file=sys.stderr)
        print('stop:', self.stop, file=sys.stderr)
        print('score:', self.score, file=sys.stderr)
        print('mismatches:', self.mismatches, file=sys.stderr)
        print('ident:', self.ident, file=sys.stderr)


def alstart(seq_a, seq_b):
    '''Returns starting position of an alignment, given two sequences'''
    import tempfile

    f = tempfile.NamedTemporaryFile()  # delete=False)
    fn = f.name
    f.close()
    needle_align(seq_a, seq_b, fn)
    af = alignfile2dict([fn], 'noname', GO, GE)
    ak1 = list(af.keys())[0]
    ak2 = list(af[ak1].keys())[0]
    a = af[ak1][ak2]
    a.summary()
    os.remove(fn)
    return a.start


def alignfile2dict(al_files, name='noname', gap_open=GO, gap_extend=GE):
    """
        Takes a list of files and returns an alignment 2D dictionary
        """
    from Bio import AlignIO

    al_set = AlignDict(name, gap_open, gap_extend)
    for f_file in al_files:
        alignment = AlignIO.read(f_file, 'fasta')
        assert len(alignment) == 2
        id_a, id_b = alignment[0].id.split('#')[0], \
            alignment[1].id.split('#')[0]
        descr_a = alignment[0].description
        descr_b = alignment[1].description
        a, b = str(alignment[0].seq).upper(), \
               str(alignment[1].seq).upper()

        al_set[id_a][id_b] = AlignInstance(id_a, id_b, a, b, descr_a, descr_b)

    return al_set


def needle_align(a_seq, b_seq, out_file, go=GO, ge=GE,
                 reverse1=False, reverse2=False):
    """ Does not require cmline, neither from Biopython,
        nor from pythonlib
    """
    import subprocess
    Verbose = False

    line1 = 'needle %s %s -outfile %s -auto' % (a_seq, b_seq, out_file)
    line1 += ' -gapopen %f -gapextend %f -aformat fasta' % (go, ge)
    if reverse1:
        line1 += ' -sreverse1'
    if reverse2:
        line1 += ' -sreverse2'
    if Verbose:
        print(line1, 'running', file=sys.stderr)
    if os.path.exists(out_file):
        return out_file, 'exists'
    try:
        retcode = subprocess.call(line1, shell=True)
        if retcode < 0:
            if Verbose:
                print("'%s'" % line1, file=sys.stderr)
                print("Child diri_sampler was terminated by signal", -retcode, file=sys.stderr)
        else:
            if Verbose:
                print("Child diri_sampler returned %i" % retcode)
    except OSError as ee:
        if Verbose:
            print("Execution of diri_sampler failed:", ee, file=sys.stderr)

    return a_seq, b_seq, 'aligned'


def main():
    '''
    Quick test
    '''
    print(' running a quick test '.center(60, '#'))
    print(' some information is printed '.center(60, '#'))
    a = 'AAAAACACACGCAT-------'
    b = '---AACAC-C-CATTTTTTTT'
    print(a)
    print(b)
    ai = AlignInstance('a', 'b', a, b, 120.0)
    ai.summary()
    ai.print_info()
    print(ai.ident)

if __name__ == '__main__':
    main()
