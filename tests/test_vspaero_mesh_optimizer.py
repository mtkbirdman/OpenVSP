import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.VSPAEROMesh import (
    _evaluate_candidate_pair,
    _qoi_pair_convergence,
    _same_topology_family,
    _scaled_tessellation_count,
    _select_converged_case,
    _set_saved_vspaero_representation,
    resolve_vspaero_representation,
)


class FakeVSP:
    SET_NONE = -1

    def __init__(self):
        self.sets = {
            "ThinGeom": (3, ["wing", "tail"]),
            "ThickGeom": (4, ["body"]),
            "ThickAll": (5, ["wing", "tail", "body"]),
        }
        self.names = {"wing": "Wing A", "tail": "Tail A", "body": "Body A"}
        self.parm_values = {}

    def GetSetIndex(self, name):
        return self.sets.get(name, (self.SET_NONE, []))[0]

    def GetGeomSetAtIndex(self, set_index):
        for _, (index, members) in self.sets.items():
            if index == set_index:
                return members
        return []

    def GetGeomName(self, geom_id):
        return self.names[geom_id]

    def FindContainer(self, name, index):
        return "settings" if name == "VSPAEROSettings" and index == 0 else ""

    def FindContainerParmIDs(self, container_id):
        return ["geom_set_parm", "thin_geom_set_parm"] if container_id == "settings" else []

    def GetParmName(self, parm_id):
        return {"geom_set_parm": "GeomSet", "thin_geom_set_parm": "ThinGeomSet"}[parm_id]

    def SetParmVal(self, parm_id, value):
        self.parm_values[parm_id] = value

    def Update(self):
        pass


def test_representation_resolution_is_model_name_agnostic():
    vsp = FakeVSP()

    hybrid = resolve_vspaero_representation(vsp, "hybrid")
    assert hybrid["name"] == "Hybrid"
    assert hybrid["thin_set_index"] == 3
    assert hybrid["thick_set_index"] == 4
    assert hybrid["active_geom_names"] == ["Wing A", "Tail A", "Body A"]

    thin = resolve_vspaero_representation(vsp, "ThinWing")
    assert thin["thin_set_index"] == 3
    assert thin["thick_set_index"] == vsp.SET_NONE
    assert thin["active_geom_names"] == ["Wing A", "Tail A"]

    thick = resolve_vspaero_representation(vsp, "ThickWing")
    assert thick["thin_set_index"] == vsp.SET_NONE
    assert thick["thick_set_index"] == 3

    thick_all = resolve_vspaero_representation(vsp, "ThickAll")
    assert thick_all["thick_set_index"] == 5
    assert thick_all["active_geom_names"] == ["Wing A", "Tail A", "Body A"]




def test_selected_representation_is_persisted_to_vspaero_settings():
    vsp = FakeVSP()
    representation = resolve_vspaero_representation(vsp, "Hybrid")
    _set_saved_vspaero_representation(vsp, representation)
    assert vsp.parm_values["geom_set_parm"] == pytest.approx(4.0)
    assert vsp.parm_values["thin_geom_set_parm"] == pytest.approx(3.0)


def test_representation_can_use_nondefault_set_names():
    vsp = FakeVSP()
    vsp.sets = {
        "Lifting": (7, ["wing", "tail"]),
        "Bodies": (8, ["body"]),
        "Everything": (9, ["wing", "tail", "body"]),
    }
    result = resolve_vspaero_representation(
        vsp,
        "Hybrid",
        lifting_set_name="Lifting",
        body_set_name="Bodies",
        thick_all_set_name="Everything",
    )
    assert result["thin_set_index"] == 7
    assert result["thick_set_index"] == 8


def test_missing_representation_set_fails_loudly():
    vsp = FakeVSP()
    del vsp.sets["ThickGeom"]
    with pytest.raises(ValueError, match="ThickGeom"):
        resolve_vspaero_representation(vsp, "Hybrid")


def test_tessellation_scaling_scales_intervals_not_raw_count():
    assert _scaled_tessellation_count(33, 0.5) == 17
    assert _scaled_tessellation_count(33, 1.0) == 33
    assert _scaled_tessellation_count(33, 1.5) == 49
    assert _scaled_tessellation_count(2, 0.1) == 2


def test_qoi_convergence_uses_user_scale_for_near_zero_quantity():
    coarse = {"CLiw": 0.400, "CMytot": 0.0010}
    fine = {"CLiw": 0.402, "CMytot": 0.0015}
    passed, errors = _qoi_pair_convergence(
        coarse,
        fine,
        {"CLiw": 0.01, "CMytot": 0.02},
        {"CMytot": 0.05},
    )
    assert errors["CLiw"] < 0.01
    assert errors["CMytot"] == pytest.approx(0.01)
    assert passed is True


def topology_signature(junction_pairs=((1, 1), (2, 2))):
    return {
        "surfaces": [[1, "Wing"], [2, "Body"]],
        "junction_pairs": [[list(pair) for pair in junction_pairs]],
        "kutta_surfaces": [[1, "Wing", True]],
        "same_surface_non_manifold_edge_count": 0,
    }


def candidate(scale, cl, *, spikes=0, feasible=True, topology=None):
    return {
        "status": "completed",
        "feasible": feasible,
        "scale": scale,
        "CLiw": cl,
        "strong_mesh_advisory_count": 0,
        "n_local_cp_spikes": spikes,
        "lod_outlier_count": 0,
        "topology_signature_json": json.dumps(topology or topology_signature()),
        "case_name": f"case_{scale}",
        "vsp3_path": f"case_{scale}.vsp3",
    }


def test_topology_family_ignores_refinement_count_but_detects_structural_change():
    first = topology_signature()
    second = topology_signature()
    assert _same_topology_family(first, second) is True

    changed = topology_signature(junction_pairs=((1, 1), (3, 4)))
    assert _same_topology_family(first, changed) is False


def test_candidate_pair_rejects_coarse_mesh_when_refinement_removes_cp_spike():
    coarse = candidate(0.6, 0.400, spikes=1)
    fine = candidate(1.0, 0.401, spikes=0)
    result = _evaluate_candidate_pair(coarse, fine, {"CLiw": 0.01}, None)
    assert result["qoi_converged"] is True
    assert result["local_health_not_worse"] is False
    assert result["passed"] is False


def test_selector_returns_coarsest_case_confirmed_by_two_finer_levels():
    rows = [
        candidate(0.6, 0.4000),
        candidate(1.0, 0.4010),
        candidate(1.5, 0.4015),
    ]
    result = _select_converged_case(rows, {"CLiw": 0.01}, None)
    assert result["converged"] is True
    assert result["selected"]["scale"] == pytest.approx(0.6)


def test_selector_returns_finest_feasible_case_when_convergence_not_demonstrated():
    rows = [
        candidate(0.6, 0.30),
        candidate(1.0, 0.38),
        candidate(1.5, 0.40),
    ]
    result = _select_converged_case(rows, {"CLiw": 0.01}, None)
    assert result["converged"] is False
    assert result["selected"]["scale"] == pytest.approx(1.5)


def test_selector_does_not_skip_failed_intermediate_level_when_claiming_convergence():
    failed = candidate(1.0, 0.4005, feasible=False)
    failed["status"] = "failed"
    rows = [
        candidate(0.6, 0.4000),
        failed,
        candidate(1.5, 0.4008),
    ]
    result = _select_converged_case(rows, {"CLiw": 0.01}, None)
    assert result["converged"] is False
    assert result["selected"]["scale"] == pytest.approx(1.5)


def test_selector_does_not_ignore_a_later_finer_level_that_diverges():
    rows = [
        candidate(0.5, 0.4000),
        candidate(0.8, 0.4005),
        candidate(1.0, 0.4008),
        candidate(1.5, 0.4300),
    ]
    result = _select_converged_case(rows, {"CLiw": 0.01}, None)
    assert result["converged"] is False
    assert result["selected"]["scale"] == pytest.approx(1.5)
