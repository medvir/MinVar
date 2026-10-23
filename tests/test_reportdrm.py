#!/usr/bin/env python3
"""Simple test facility"""
import pandas as pd
import pytest

from src.minvar.reportdrm import annotate_hiv_sierra


def test_annotate_hiv_sierra():
    df = pd.DataFrame([
        {'gene': 'GagPolTF', 'wt': 'D', 'pos': 5, 'mut': 'N', 'freq': 1.0},
        {'gene': 'protease', 'wt': 'L', 'pos': 10, 'mut': 'I', 'freq': 0.0809},
        {'gene': 'protease', 'wt': 'K', 'pos': 20, 'mut': 'R', 'freq': 0.983},
        {'gene': 'protease', 'wt': 'L', 'pos': 89, 'mut': 'M', 'freq': 1.0},
        {'gene': 'RT', 'wt': 'D', 'pos': 123, 'mut': '*', 'freq': 0.0179},
        {'gene': 'RT', 'wt': 'V', 'pos': 179, 'mut': 'I', 'freq': 1.0},
        {'gene': 'RT', 'wt': 'D', 'pos': 192, 'mut': 'N', 'freq': 0.0153},
        {'gene': 'capsid', 'wt': 'M', 'pos': 66, 'mut': 'I', 'freq': 0.25},
    ])
    out, version, pubdate = annotate_hiv_sierra(df, sample_id='test')
    assert len(out) == 8, out
    assert version
    assert pubdate
    assert list(out['category']) == [
        'unannotated', 'Other', 'Other', 'unannotated',
        'unannotated', 'Other', 'unannotated', 'Major'
    ]


@pytest.mark.skip(reason='Still an empty box')
def test_parse_merged_2():
    '''
    A, 1011, G
    A, 1050, T
    '''
    assert 1 > 0
