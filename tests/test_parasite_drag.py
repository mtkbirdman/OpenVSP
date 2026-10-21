from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import src.AnalysisParasiteDrag as parasite_drag_module
from src.AnalysisParasiteDrag import (
    apply_parasite_drag_to_polar,
    vsp_parasite_drag,
)


class FakeParasiteDragVSP:
    V_UNIT_M_S = 1
    PD_UNITS_METRIC = 1
    TEMP_UNIT_K = 0
    ATMOS_TYPE_US_STANDARD_1976 = 0
    MANUAL_REF = 0

    def __init__(self):
        self.inputs = {
            'VelocityUnit': [0],
            'AltLengthUnit': [0],
            'TempUnit': [0],
            'FreestreamPropChoice': [0],
            'Vinf': [0.0],
            'Altitude': [0.0],
            'DeltaTemp': [0.0],
            'RecomputeGeom': [1],
            'RefFlag': [1],
            'Sref': [100.0],
            'GeomSet': [0],
            'FileName': ['parasite.csv'],
        }
        self.recompute_history = []
        self.result_inputs = {}
        self.exec_count = 0

    def SetAnalysisInputDefaults(self, analysis_name):
        assert analysis_name == 'ParasiteDrag'

    def GetAnalysisInputNames(self, analysis_name):
        assert analysis_name == 'ParasiteDrag'
        return list(self.inputs)

    def SetIntAnalysisInput(self, analysis_name, name, values, index=0):
        assert analysis_name == 'ParasiteDrag'
        self.inputs[name] = list(values)
        if name == 'RecomputeGeom':
            self.recompute_history.append(int(values[0]))

    def SetDoubleAnalysisInput(self, analysis_name, name, values, index=0):
        assert analysis_name == 'ParasiteDrag'
        self.inputs[name] = list(values)

    def SetStringAnalysisInput(self, analysis_name, name, values, index=0):
        assert analysis_name == 'ParasiteDrag'
        self.inputs[name] = list(values)

    def PrintAnalysisInputs(self, analysis_name):
        assert analysis_name == 'ParasiteDrag'

    def ExecAnalysis(self, analysis_name):
        assert analysis_name == 'ParasiteDrag'
        self.exec_count += 1
        result_id = f'pd-{self.exec_count}'
        self.result_inputs[result_id] = {
            key: values[0] for key, values in self.inputs.items()
        }
        return result_id

    def GetDoubleResults(self, result_id, result_name, index=0):
        inputs = self.result_inputs[result_id]
        velocity = float(inputs['Vinf'])
        values = {
            'FC_Vinf': velocity,
            'FC_Alt': float(inputs['Altitude']),
            'FC_Mach': velocity / 340.0,
            'Re_L': velocity * 1.0e6,
            'FC_Sref': float(inputs['Sref']),
            'Total_CD_Total': 0.01 + velocity * 1.0e-5,
            'Analysis_Duration_Sec': 0.02,
        }
        return [values[result_name]] if result_name in values else []


def test_vsp_parasite_drag_runs_requested_speeds_and_reuses_geometry():
    fake_vsp = FakeParasiteDragVSP()

    result = vsp_parasite_drag(
        fake_vsp,
        [30.0, 40.0],
        altitude_m=1200.0,
        dT=5.0,
        sref=18.0,
        geom_set=3,
        verbose=0,
    )

    assert list(result.columns) == [
        'Velocity_mps',
        'Altitude_m',
        'Mach',
        'Re_L',
        'Sref',
        'CD0',
        'AnalysisDuration_s',
    ]
    np.testing.assert_allclose(result['Velocity_mps'], [30.0, 40.0])
    np.testing.assert_allclose(result['Altitude_m'], [1200.0, 1200.0])
    np.testing.assert_allclose(result['Sref'], [18.0, 18.0])
    np.testing.assert_allclose(result['CD0'], [0.0103, 0.0104])
    assert fake_vsp.recompute_history == [1, 0]
    assert fake_vsp.inputs['VelocityUnit'] == [fake_vsp.V_UNIT_M_S]
    assert fake_vsp.inputs['AltLengthUnit'] == [fake_vsp.PD_UNITS_METRIC]
    assert fake_vsp.inputs['TempUnit'] == [fake_vsp.TEMP_UNIT_K]
    assert fake_vsp.inputs['FreestreamPropChoice'] == [
        fake_vsp.ATMOS_TYPE_US_STANDARD_1976
    ]
    assert fake_vsp.inputs['RefFlag'] == [fake_vsp.MANUAL_REF]
    assert fake_vsp.inputs['GeomSet'] == [3]


def test_apply_parasite_drag_to_polar_uses_cli_and_cd0_plus_cdi(monkeypatch):
    source = pd.DataFrame(
        {
            'InputMach': [0.10, 0.20],
            'CLi': [0.50, 0.80],
            'CL': [0.55, 0.85],
            'CLtot': [0.55, 0.85],
            'CDi': [0.010, 0.020],
            'CD0': [0.003, 0.004],
            'CDtot': [0.013, 0.024],
        }
    )
    original = source.copy(deep=True)

    def fake_parasite_drag(vsp, velocity_mps, **kwargs):
        assert kwargs['sref'] == 18.0
        assert len(velocity_mps) == 2
        return pd.DataFrame({'CD0': [0.020, 0.030]})

    monkeypatch.setattr(
        parasite_drag_module,
        'vsp_parasite_drag',
        fake_parasite_drag,
    )

    corrected = apply_parasite_drag_to_polar(
        object(),
        source,
        altitude_m=0.0,
        dT=0.0,
        sref=18.0,
        verbose=0,
    )

    pd.testing.assert_frame_equal(source, original)
    np.testing.assert_allclose(corrected['CL'], corrected['CLi'])
    np.testing.assert_allclose(corrected['CD0_VSPAERO'], original['CD0'])
    np.testing.assert_allclose(corrected['CD0'], [0.020, 0.030])
    np.testing.assert_allclose(corrected['CD'], [0.030, 0.050])
    np.testing.assert_allclose(corrected['L_D'], [0.50 / 0.030, 0.80 / 0.050])
    np.testing.assert_allclose(corrected['CLtot'], original['CLtot'])
    np.testing.assert_allclose(corrected['CDi'], original['CDi'])
    np.testing.assert_allclose(corrected['CDtot'], original['CDtot'])


def test_apply_parasite_drag_to_polar_reads_vspaero_sref(monkeypatch):
    captured = {}

    monkeypatch.setattr(
        parasite_drag_module,
        'find_vspaero_settings_container',
        lambda vsp: 'settings-id',
    )
    monkeypatch.setattr(
        parasite_drag_module,
        'get_container_parm_value',
        lambda vsp, container_id, parm_name: (17.8819, 'sref-parm'),
    )

    def fake_parasite_drag(vsp, velocity_mps, **kwargs):
        captured['sref'] = kwargs['sref']
        return pd.DataFrame({'CD0': [0.015]})

    monkeypatch.setattr(
        parasite_drag_module,
        'vsp_parasite_drag',
        fake_parasite_drag,
    )

    corrected = apply_parasite_drag_to_polar(
        object(),
        pd.DataFrame({'Mach': [0.1], 'CLi': [0.6], 'CDi': [0.012]}),
        verbose=0,
    )

    assert captured['sref'] == pytest.approx(17.8819)
    assert corrected['CD'].iloc[0] == pytest.approx(0.027)


def test_parasite_drag_g103a_integration():
    vsp = pytest.importorskip('openvsp')
    model_path = (
        REPO_ROOT
        / 'examples'
        / 'notebooks'
        / 'trimmed_polar'
        / 'G103A.vsp3'
    )

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(model_path))
    vsp.Update()

    result = vsp_parasite_drag(
        vsp,
        [30.0],
        altitude_m=0.0,
        sref=17.88185312625759,
        verbose=0,
    )

    assert len(result) == 1
    assert np.isfinite(result[['Mach', 'Re_L', 'Sref', 'CD0']].to_numpy()).all()
    assert result['CD0'].iloc[0] > 0.0
