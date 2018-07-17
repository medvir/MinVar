#!/usr/bin/env python3
"""Simple test facility"""
import os
import pytest

from src.minvar.reportdrm import parse_merged


test_mer = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'small_merged_muts.csv')
assert os.path.exists(test_mer)


def test_parse_merged():
    out = parse_merged(test_mer)
    assert len(out) == 8, out
    
    # assert out_0['wt'] == 'T', out
    # assert out_0['pos'] == 4, out
    # assert out_0['mut'] == 'A', out
    # assert out_1['wt'] == 'G', out
    # assert out_1['pos'] == 11, out
    # assert out_1['mut'] == 'C', out
    # assert out.iloc[2]['wt'] == 'G', out
    # assert out.iloc[2]['mut'] == 'CAAA', out
    # assert out.iloc[3]['wt'] == 'TAGG', out
    # assert out.iloc[3]['mut'] == 'C', out


@pytest.mark.skip(reason='Still an empty box')
def test_parse_merged_2():
    import pandas as pd
    '''
    A, 1011, G
    A, 1050, T
    '''
    assert 1 > 0