from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.VSPAEROMeshDataset import (
    _build_experiment_cases,
    _case_id,
    _latin_hypercube,
    _parameter_levels,
    discover_vspaero_mesh_parameters,
)


class FakeDatasetVSP:
    def __init__(self):
        self.geom_parms = {
            "wing": ["wing_tw", "wing_le", "wing_te", "wing_cap", "wing_chord"],
            "body": ["body_tw", "body_length"],
        }
        self.names = {
            "wing_tw": "Tess_W",
            "wing_le": "LECluster",
            "wing_te": "TECluster",
            "wing_cap": "CapUMinTess",
            "wing_chord": "Chord",
            "body_tw": "Tess_W",
            "body_length": "Length",
            "wing_x1_u": "SectTess_U",
            "wing_x1_in": "InCluster",
            "wing_x1_out": "OutCluster",
            "body_x1_u": "SectTess_U",
            "body_x1_fwd": "FwdCluster",
            "body_x1_aft": "AftCluster",
        }
        self.values = {
            "wing_tw": 33.0,
            "wing_le": 0.5,
            "wing_te": 1.0,
            "wing_cap": 5.0,
            "wing_chord": 8.0,
            "body_tw": 17.0,
            "body_length": 20.0,
            "wing_x1_u": 12.0,
            "wing_x1_in": 1.0,
            "wing_x1_out": 0.5,
            "body_x1_u": 10.0,
            "body_x1_fwd": 1.0,
            "body_x1_aft": 1.0,
        }
        self.xsec_parms = {
            ("wing_x0", "SectTess_U"): "",
            ("wing_x1", "SectTess_U"): "wing_x1_u",
            ("wing_x1", "InCluster"): "wing_x1_in",
            ("wing_x1", "OutCluster"): "wing_x1_out",
            ("body_x0", "SectTess_U"): "",
            ("body_x1", "SectTess_U"): "body_x1_u",
            ("body_x1", "FwdCluster"): "body_x1_fwd",
            ("body_x1", "AftCluster"): "body_x1_aft",
        }

    def GetGeomName(self, geom_id):
        return {"wing": "Wing", "body": "Body"}[geom_id]

    def GetGeomTypeName(self, geom_id):
        return {"wing": "WingGeom", "body": "FuselageGeom"}[geom_id]

    def FindContainerParmIDs(self, geom_id):
        return self.geom_parms.get(geom_id, [])

    def GetParmName(self, parm_id):
        return self.names[parm_id]

    def GetParmGroupName(self, parm_id):
        return "Shape"

    def GetParmVal(self, parm_id):
        return self.values[parm_id]

    def GetParmLowerLimit(self, parm_id):
        return 2.0 if self.names[parm_id] in {"Tess_W", "SectTess_U", "CapUMinTess"} else 0.0

    def GetParmUpperLimit(self, parm_id):
        return 101.0 if self.names[parm_id] in {"Tess_W", "SectTess_U", "CapUMinTess"} else 2.0

    def GetNumXSecSurfs(self, geom_id):
        return 1

    def GetXSecSurf(self, geom_id, surf_index):
        return f"{geom_id}_surf"

    def GetNumXSec(self, surf_id):
        return 2

    def GetXSec(self, surf_id, xsec_index):
        geom_id = surf_id.removesuffix("_surf")
        return f"{geom_id}_x{xsec_index}"

    def GetXSecParm(self, xsec_id, name):
        return self.xsec_parms.get((xsec_id, name), "")


REPRESENTATION = {
    "active_geom_ids": ["wing", "body"],
    "thin_geom_ids": ["wing"],
    "thick_geom_ids": ["body"],
}


def test_parameter_discovery_uses_mesh_whitelist_and_keeps_per_section_controls():
    rows = discover_vspaero_mesh_parameters(FakeDatasetVSP(), REPRESENTATION)
    by_key = {row["parameter_key"]: row for row in rows}

    assert "geom:wing:Tess_W" in by_key
    assert "geom:wing:LECluster" in by_key
    assert "geom:wing:TECluster" in by_key
    assert "geom:wing:CapUMinTess" in by_key
    assert "xsec:wing:0:1:SectTess_U" in by_key
    assert "xsec:wing:0:1:InCluster" in by_key
    assert "xsec:wing:0:1:OutCluster" in by_key
    assert "xsec:body:0:1:FwdCluster" in by_key
    assert "xsec:body:0:1:AftCluster" in by_key

    # Geometry design parameters must never enter the mesh DOE automatically.
    assert all(row["parameter_name"] not in {"Chord", "Length"} for row in rows)


def parameter(key, kind, baseline, lower, upper):
    return {
        "parameter_key": key,
        "parameter_kind": kind,
        "baseline_value": baseline,
        "lower_limit": lower,
        "upper_limit": upper,
    }


def test_count_levels_scale_intervals_and_respect_parameter_limits():
    row = parameter("w", "count", 33.0, 9.0, 49.0)
    assert _parameter_levels(row, (0.25, 0.5, 1.0, 1.5, 2.0), (1.0,)) == [9.0, 17.0, 33.0, 49.0]


def test_zero_continuous_baseline_still_gets_nonzero_levels_from_finite_limits():
    row = parameter("cluster", "continuous", 0.0, 0.0, 2.0)
    levels = _parameter_levels(row, (1.0,), (0.5, 1.0, 1.5))
    assert levels == [0.0, 0.5, 1.0, 1.5]


def test_latin_hypercube_is_deterministic_and_stratified_per_dimension():
    first = _latin_hypercube(8, 3, 42)
    second = _latin_hypercube(8, 3, 42)
    assert np.allclose(first, second)
    assert np.all((first >= 0.0) & (first < 1.0))
    for dimension in range(first.shape[1]):
        bins = np.floor(first[:, dimension] * 8).astype(int)
        assert sorted(bins.tolist()) == list(range(8))


def test_experiment_plan_sweeps_every_parameter_and_lhs_varies_all_dimensions():
    parameters = [
        parameter("wing_w", "count", 33.0, 2.0, 101.0),
        parameter("body_w", "count", 17.0, 2.0, 101.0),
        parameter("wing_le", "continuous", 0.5, 0.0, 2.0),
    ]
    cases = _build_experiment_cases(
        parameters,
        count_scales=(0.5, 1.0, 1.5),
        continuous_scales=(0.5, 1.0, 1.5),
        lhs_samples=12,
        lhs_seed=7,
        include_single_parameter_sweeps=True,
        include_pairwise_extremes=False,
        explicit_cases=None,
    )

    baseline = next(case for case in cases if case["experiment_type"] == "baseline")
    for key in ("wing_w", "body_w", "wing_le"):
        assert any(
            case["experiment_type"] == "single" and case["values"][key] != baseline["values"][key]
            for case in cases
        )

    lhs = [case for case in cases if str(case["experiment_type"]).startswith("lhs_")]
    assert len(lhs) == 12
    for key in baseline["values"]:
        assert len({case["values"][key] for case in lhs}) > 1


def test_pairwise_extreme_block_is_available_without_being_mandatory():
    parameters = [
        parameter("a", "count", 17.0, 2.0, 101.0),
        parameter("b", "continuous", 1.0, 0.0, 2.0),
        parameter("c", "continuous", 0.5, 0.0, 2.0),
    ]
    without_pairwise = _build_experiment_cases(
        parameters,
        count_scales=(0.5, 1.0, 1.5),
        continuous_scales=(0.5, 1.0, 1.5),
        lhs_samples=0,
        lhs_seed=0,
        include_single_parameter_sweeps=False,
        include_pairwise_extremes=False,
        explicit_cases=None,
    )
    with_pairwise = _build_experiment_cases(
        parameters,
        count_scales=(0.5, 1.0, 1.5),
        continuous_scales=(0.5, 1.0, 1.5),
        lhs_samples=0,
        lhs_seed=0,
        include_single_parameter_sweeps=False,
        include_pairwise_extremes=True,
        explicit_cases=None,
    )
    assert len(without_pairwise) == 1
    assert len(with_pairwise) == 1 + 4 * 3  # baseline + four corners for each of three pairs


def test_case_id_is_order_independent_but_changes_with_solver_context():
    common = dict(
        source_sha256="abc",
        openvsp_version="OpenVSP 3.52.1",
        representation="ThickAll",
        flight_condition={"alpha_deg": 2.0, "mach": 0.1, "reynolds_number": 4.4e6},
    )
    first = _case_id(
        **common,
        values={"b": 2.0, "a": 1.0},
        solver_settings={"wake_num_iter": 12},
    )
    second = _case_id(
        **common,
        values={"a": 1.0, "b": 2.0},
        solver_settings={"wake_num_iter": 12},
    )
    changed = _case_id(
        **common,
        values={"a": 1.0, "b": 2.0},
        solver_settings={"wake_num_iter": 16},
    )
    assert first == second
    assert first != changed


def test_case_state_uses_case_json_as_commit_marker(tmp_path):
    from src.VSPAEROMeshDataset import _load_case_states

    case_dir = tmp_path / "cases" / "abc"
    first = case_dir / "attempt_01"
    first.mkdir(parents=True)
    (first / "case.json").write_text(
        '{"case_id":"abc","attempt":1,"status":"completed","elapsed_s":12.5}',
        encoding="utf-8",
    )

    # An interrupted later attempt has artifacts but no final case.json marker.
    second = case_dir / "attempt_02"
    second.mkdir()
    (second / "request.json").write_text('{"case_id":"abc","attempt":2}', encoding="utf-8")

    state = _load_case_states(tmp_path)["abc"]
    assert state["max_attempt"] == 2
    assert state["latest_committed_attempt"] == 1
    assert state["latest_status"] == "completed"
    assert state["latest_elapsed_s"] == 12.5


def test_rebuild_dataset_tables_reads_only_committed_attempts(tmp_path):
    import json
    import pandas as pd

    from src.VSPAEROMeshDataset import rebuild_mesh_dataset_tables

    committed = tmp_path / "cases" / "abc" / "attempt_01"
    committed.mkdir(parents=True)
    (committed / "case.json").write_text(
        json.dumps({"case_id": "abc", "attempt": 1, "status": "completed"}),
        encoding="utf-8",
    )
    pd.DataFrame([
        {"case_id": "abc", "attempt": 1, "parameter_key": "geom:wing:Tess_W", "requested_value": 17.0}
    ]).to_csv(committed / "parameters.csv", index=False)

    incomplete = tmp_path / "cases" / "def" / "attempt_01"
    incomplete.mkdir(parents=True)
    (incomplete / "request.json").write_text('{"case_id":"def","attempt":1}', encoding="utf-8")
    pd.DataFrame([
        {"case_id": "def", "attempt": 1, "parameter_key": "geom:wing:Tess_W", "requested_value": 33.0}
    ]).to_csv(incomplete / "parameters.csv", index=False)

    rebuilt = rebuild_mesh_dataset_tables(tmp_path)
    cases = pd.read_csv(rebuilt["cases_path"])
    parameters = pd.read_csv(rebuilt["parameters_path"])

    assert rebuilt["committed_attempt_count"] == 1
    assert cases[["case_id", "attempt", "status"]].to_dict(orient="records") == [
        {"case_id": "abc", "attempt": 1, "status": "completed"}
    ]
    assert parameters["case_id"].tolist() == ["abc"]
