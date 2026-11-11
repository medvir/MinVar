#!/usr/bin/env python3
"""Simple test facility"""
import io
import os
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.minvar.reportdrm import annotate_hiv_sierra, calc_APD_HCV, write_ambig_score_HCV


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


def test_calc_apd_hcv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sample_id = 'test_hcv'

    # Load H77 CDS as synthetic sample sequence (3011 codons = 9033 nt)
    h77_cds_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'src/minvar/db/HCV/H77_cds.fna',
    )
    h77_seq = str(list(SeqIO.parse(h77_cds_path, 'fasta'))[0].seq)

    SeqIO.write(
        SeqRecord(Seq(h77_seq), id=f'max_freq_cons_{sample_id}', description=''),
        f'cns_max_freq_{sample_id}.fasta',
        'fasta',
    )
    # Put one IUPAC wobble 'R' at nt 1507 (first nt of codon 503)
    ambi_list = list(h77_seq)
    ambi_list[1506] = 'R'
    SeqIO.write(
        SeqRecord(Seq(''.join(ambi_list)), id=f'ambiguous_consensus_{sample_id}', description=''),
        f'cns_ambiguous_{sample_id}.fasta',
        'fasta',
    )

    # 1-to-1 codon mapping for codons 1..3011
    aa_df = pd.DataFrame({
        'pos': list(range(1, 3012)),
        'wt': ['A'] * 3011,
        'mut': ['A'] * 3011,
        'in_pos': list(range(1, 3012)),
    })
    aa_df.to_csv(f'max_freq_muts_aa_{sample_id}.csv', index=False)

    # Create merged_mutations_nt with a 80%/20% split at nt 1507 (inside both 503-1004 and 384-1026)
    nt_rows = []
    for p in range(1, len(h77_seq) + 1):
        if p == 1507:
            nt_rows.append({'pos': p, 'mut': 'A', 'freq': 0.80, 'coverage': 500})
            nt_rows.append({'pos': p, 'mut': 'G', 'freq': 0.20, 'coverage': 500})
        else:
            nt_rows.append({'pos': p, 'mut': h77_seq[p - 1], 'freq': 1.0, 'coverage': 500})
    pd.DataFrame(nt_rows).to_csv(f'merged_mutations_nt_{sample_id}.csv', index=False)

    apd_sum_503, rlen_503, apd_pos_503 = calc_APD_HCV(
        sample_id=sample_id, start_codon=503, end_codon=1004
    )
    assert rlen_503 == 1506
    assert apd_pos_503 == 1
    # Diversity contribution at 0.8 / 0.2 is 0.8*0.2 + 0.2*0.8 = 0.32
    assert pytest.approx(apd_sum_503, rel=1e-6) == 0.32

    apd_sum_e2ns2, rlen_e2ns2, apd_pos_e2ns2 = calc_APD_HCV(
        sample_id=sample_id, start_codon=384, end_codon=1026
    )
    assert rlen_e2ns2 == 1929
    assert apd_pos_e2ns2 == 1
    assert pytest.approx(apd_sum_e2ns2, rel=1e-6) == 0.32

    buf = io.StringIO()
    write_ambig_score_HCV(buf, sample_id=sample_id)
    report_text = buf.getvalue()
    assert 'ORF codons 503-1004' in report_text
    assert 'Genes E2, p7 and NS2 (codons 384-1026)' in report_text
    assert 'Estimated Time Since Infection (TSI):' in report_text
    assert 'Ambiguity score:' not in report_text
