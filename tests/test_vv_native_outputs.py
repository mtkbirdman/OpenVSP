"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path

import pytest

from openvsppy.vv.VVNative import check_stability_condition
from openvsppy.vv.VVNative import read_mode_derivatives
from openvsppy.vv.VVNative import required_native_files
from openvsppy.vv.VVPostprocess import condition_match
from openvsppy.vv.VVPostprocess import read_history
from openvsppy.vv.VVPostprocess import read_polar

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/fixtures/SampleGlider_native_3.50.5.stab'
HEADER = FIXTURE.read_text().split('Case')[0]
MAPPING = {'ConGrp_1':'AILERON_GROUP','ConGrp_2':'ELEVATOR_GROUP','ConGrp_3':'RUDDER_GROUP'}
COEFS = ('CFx','CFy','CFz','CMx','CMy','CMz','CL','CD','CS','CMl','CMm','CMn')


def test_nonfinite_optional_polar_value_is_not_turned_into_zero(tmp_path):
    from openvsppy.vv.VVPostprocess import POLAR_CONDITIONS
    reference = ROOT / 'examples/models/CRM-HS/CRM-HS.polar'
    content = reference.read_text()
    lines = content.splitlines()
    header = next(line.split() for line in lines if line.split()[:3] == ['Beta', 'Mach', 'AoA'])
    data_index = next(i for i, line in enumerate(lines) if len(line.split()) == len(header)
                      and line.split()[0] == '0.000000000000')
    data = lines[data_index].split()
    data[header.index('L/D')] = 'nan'
    lines[data_index] = ' '.join(data)
    path = tmp_path / 'nonfinite_ratio.polar'
    path.write_text('\n'.join(lines), encoding='utf-8')
    row = read_polar(path)[0]
    assert row['L/D'] is None
    assert row['CLi'] is not None
    assert condition_match(row, {'mach': 0.6, 're_cref': 1e7, 'alpha_deg': 3., 'beta_deg': 0.})
    assert all(row[key] is not None for key in POLAR_CONDITIONS if key in row)

    data[header.index('Mach')] = 'nan'
    lines[data_index] = ' '.join(data)
    path.write_text('\n'.join(lines), encoding='utf-8')
    with pytest.raises(ValueError, match='Nonfinite polar condition'):
        read_polar(path)


def test_offline_direct_evidence_parsers_and_wake_count():
    from openvsppy.vv.VVPostprocess import read_polar,read_history,read_vspaero_settings
    evidence=Path(__file__).resolve().parents[1]/'vv_study/02_reference/20261031_api_vs_direct'
    polar=read_polar(evidence/'direct_success_outputs/model_effective.polar')
    history=read_history(evidence/'direct_success_outputs/model_effective.history')
    setup=read_vspaero_settings(evidence/'api_failed_inputs/model_effective.vspaero')
    assert len(polar)==1 and len(history)==20
    assert int(history[-1]['Iter'])==20
    assert polar[0]['CLtot']==pytest.approx(.509706605733,abs=1e-10)
    assert polar[0]['CDtot']==pytest.approx(.018312180664,abs=1e-10)
    assert setup['WakeIters']==20
    assert setup['NumWakeNodes']==8


def test_file_contracts():
    assert required_native_files('STEADY_SWEEP')==('.polar','.history','.lod','.adb')
    assert required_native_files('STEADY')==('.stab',)
    assert required_native_files('ADJOINT')==('.adjoint.stab',)
    assert required_native_files('P')==('.pstab',)
    assert required_native_files('Q')==('.qstab',)
    assert required_native_files('R')==('.rstab',)


def test_steady_native_and_reference_checks():
    header, rows, groups=read_mode_derivatives(FIXTURE,'STEADY')
    assert groups==MAPPING and len(rows)==132
    check_stability_condition(header,{'mach':0,'alpha_deg':0,'beta_deg':0,
        'sref':18,'bref':27,'cref':.625})
    with pytest.raises(ValueError,match='condition mismatch'):
        check_stability_condition(header,{'mach':.1,'alpha_deg':0,'beta_deg':0,
            'sref':18,'bref':27,'cref':.625})


@pytest.mark.parametrize('mode,suffix', [('P','_p'),('Q','_(q + alpha_dot)'),('R','_(r - beta_dot)')])
def test_pqr_native_names(tmp_path,mode,suffix):
    path=tmp_path/'native'
    path.write_text(HEADER+'Case header\n# Name Value\n'+
        ''.join(f'{c}{suffix}  {idx}.125\n' for idx,c in enumerate(COEFS)),encoding='utf-8')
    _,rows,_=read_mode_derivatives(path,mode)
    assert len(rows)==12
    assert {key for key,_,_ in rows}=={c+suffix for c in COEFS}
    path.write_text(path.read_text().replace('CMn'+suffix, 'BAD'+suffix),encoding='utf-8')
    with pytest.raises(ValueError,match='Expected 12'):
        read_mode_derivatives(path,mode)


def test_adjoint_multitable_and_incomplete_rejection(tmp_path):
    path=tmp_path/'native.adjoint.stab'
    column_names='Alpha Beta p q r Mach U ConGrp_1 ConGrp_2 ConGrp_3'
    text=HEADER+'Case header\n'
    for flavor in ('Total','Inviscid','Viscous'):
        text+=f'Coef {flavor} {column_names}\n#\n#\n'
        text+=''.join(f'{coef} '+(' 0.123'*11)+'\n' for coef in COEFS)
    path.write_text(text,encoding='utf-8')
    header,rows,groups=read_mode_derivatives(path,'ADJOINT')
    assert len(rows)==396 and groups=={}
    assert ('CFx','Total_Alpha',.123) in rows
    path.write_text(text.replace('CMn '+' 0.123'*11+'\n','',1),encoding='utf-8')
    with pytest.raises(ValueError,match='Incomplete adjoint'):
        read_mode_derivatives(path,'ADJOINT')


def test_named_native_history_and_column_reorder(tmp_path):
    f=tmp_path/'sample.history'
    f.write_text('Solver Case: 1\nIter CDiw CLi CLiw CDi CMiy L2 Residual Max Residual Wall_Time\n'
                 '1 0.012 0.4 0.5 0.01 -0.02 0.001 0.002 1.0\n'
                 'Solver Case: 2\nIter CLi CDiw CLiw CDi CMiy L2 Residual Max Residual Wall_Time\n'
                 '1 0.8 0.022 0.9 0.02 -0.03 0.003 0.004 1.5\n')
    rec=read_history(f)
    assert [r['case'] for r in rec]==[1,2]
    assert rec[0]['CDiw']==pytest.approx(.012)
    assert rec[1]['CDiw']==pytest.approx(.022)
    assert rec[0]['l2_residual']==pytest.approx(.001)


def test_history_broken_header_blocks_false_evaluation(tmp_path):
    path=tmp_path/'bad.history'
    path.write_text('Solver Case: 1\nIter CLi\n1 .2\n')
    with pytest.raises(ValueError,match='no valid data'):
        read_history(path)


def test_polar_and_stab_samples():
    rows=read_polar(ROOT/'examples/models/CRM-HS/CRM-HS.polar')
    assert len(rows)==1
    assert rows[0]['Mach']==pytest.approx(.6)
    assert {'CLi','CDi','CLiw','CDiw','CDo'}.issubset(rows[0])
    assert condition_match(rows[0],{'mach':.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.})
    assert not condition_match(rows[0],{'mach':.6,'re_cref':1e7,'alpha_deg':2.,'beta_deg':0.})
    from openvsppy.vspaero.VSPAEROStab import read_vspaero_stab
    stab=read_vspaero_stab(ROOT/'tests/fixtures/SampleGlider_native_3.50.5.stab')
    assert stab.derivatives.shape[0]==12
    assert stab.control_groups['ConGrp_2']=='ELEVATOR_GROUP'


def test_native_history_parser():
    from openvsppy.vv.VVPostprocess import read_history
    rows=read_history(ROOT/'tests/fixtures/G103A_native_history_3.50.5.history')
    assert [row['iteration'] for row in rows]==[1,2,3,4,5]
    assert rows[-1]['CLi']==pytest.approx(.509540463264)
    assert rows[-1]['CDiw']==pytest.approx(.005234956441)
    assert rows[-1]['CMiy']==pytest.approx(-.058901465952)
    assert rows[-1]['wall_time']==pytest.approx(75.88)


def test_native_mesh_summary():
    from openvsppy.vv.VVPostprocess import read_mesh_stats
    sample=ROOT/'tests/fixtures/G103A_native_stdout_3.50.5.log'
    result=read_mesh_stats(sample)
    assert result['triangles_total']==10316
    assert result['kutta_nodes']==254
    assert result['wake_trailing_nodes']==8
