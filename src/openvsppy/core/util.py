"""Low-level OpenVSP helper functions shared by the analysis modules.

This module deliberately contains only thin OpenVSP API operations: import,
stdout/workdir handling, Geom/Parm lookup and access, XSec parameter access, simple
AnalysisInput setting, control-surface setup, and Result-to-DataFrame reading.
Higher-level workflows such as VSPAERO stability runs, Vv-Gamma sweeps, .stab
parsing, and turn trim solvers remain in their own modules.
"""
from __future__ import annotations

import importlib
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

# Common OpenVSP wing-section parameter names.  The first name is preferred;
# aliases keep the code usable across small API/XML naming differences.
SPAN_NAMES = ("Span",)
ROOT_CHORD_NAMES = ("Root_Chord", "RootC")
TIP_CHORD_NAMES = ("Tip_Chord", "TipC")
SWEEP_NAMES = ("Sweep",)
SWEEP_LOCATION_NAMES = ("Sweep_Location", "SweepLoc")
DIHEDRAL_NAMES = ("Dihedral",)
X_LOCATION_NAMES = ("X_Rel_Location", "X_Location")
XFORM_GROUP_NAMES = ("XForm",)

def import_openvsp():
    try:
        return importlib.import_module("openvsp")
    except ImportError as exc:
        raise ImportError(
            "openvsp is required for .vsp3 geometry editing and VSPAERO runs. "
            "Use this code inside an OpenVSP Python environment."
        ) from exc

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

@contextmanager
def workdir(path: str | os.PathLike):
    previous = Path.cwd()
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)

def find_one_geom(vsp, geom_name: str) -> str:
    geom_ids = list(vsp.FindGeomsWithName(geom_name))
    if len(geom_ids) != 1:
        raise ValueError(f"Geom '{geom_name}' must exist exactly once. found={geom_ids}")
    return geom_ids[0]

def find_container_parm(vsp, container_id: str, parm_name: str) -> str:
    for parm_id in vsp.FindContainerParmIDs(container_id):
        if parm_name == vsp.GetParmName(parm_id):
            return parm_id
    return ""

def find_vspaero_settings_container(vsp) -> str:
    """Return the VSPAERO settings ParmContainer ID independent of its saved name.

    Current OpenVSP models normally name this ParmContainer ``VSPAEROSettings``,
    but some valid legacy/external .vsp3 files save the same container as
    ``Default``.  ``FindContainer`` searches the ParmContainer *name*, so a
    name-only lookup incorrectly reports that those models have no VSPAERO
    settings.  Fall back to identifying the unique container by the core
    VSPAERO parameter set it owns.
    """

    settings_id = vsp.FindContainer("VSPAEROSettings", 0)
    if settings_id:
        return settings_id

    required_parm_names = {
        "GeomSet",
        "ThinGeomSet",
        "Sref",
        "bref",
        "cref",
        "Xcg",
        "Ycg",
        "Zcg",
    }
    candidates = []
    for container_id in vsp.FindContainers():
        try:
            parm_names = {
                str(vsp.GetParmName(parm_id))
                for parm_id in vsp.FindContainerParmIDs(container_id)
            }
        except Exception:
            continue
        if required_parm_names.issubset(parm_names):
            candidates.append(str(container_id))

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise RuntimeError(
            "Multiple ParmContainers expose the core VSPAERO settings; "
            f"candidates={candidates}"
        )
    return ""

def get_container_parm_value(vsp, container_id: str, parm_name: str) -> tuple[float | None, str]:
    parm_id = find_container_parm(vsp, container_id, parm_name)
    if not parm_id:
        return None, ""
    return vsp.GetParmVal(parm_id), parm_id

def find_geom_parm(
    vsp,
    geom_id: str,
    parm_names: Sequence[str],
    group_names: Sequence[str],
) -> str:
    """Return the first matching Geom Parm ID, or an empty string."""
    for group_name in group_names:
        for parm_name in parm_names:
            parm_id = vsp.FindParm(geom_id, parm_name, group_name)
            if parm_id and str(parm_id).upper() != "NONE":
                return parm_id
    return ""

def get_geom_parm_value(
    vsp,
    geom_id: str,
    parm_names: Sequence[str],
    group_names: Sequence[str],
    default: float | None = None,
) -> float:
    """Read a Geom Parm using ordered Parm-name and group-name candidates."""
    parm_id = find_geom_parm(vsp, geom_id, parm_names, group_names)
    if parm_id:
        return float(vsp.GetParmVal(parm_id))
    if default is not None:
        return float(default)
    raise KeyError(
        "None of these Geom parameters were found: "
        f"parms={list(parm_names)}, groups={list(group_names)}"
    )

def set_geom_parm_value(
    vsp,
    geom_id: str,
    parm_names: Sequence[str],
    group_names: Sequence[str],
    value: float,
) -> None:
    """Set a Geom Parm selected from ordered Parm-name and group-name candidates."""
    parm_id = find_geom_parm(vsp, geom_id, parm_names, group_names)
    if not parm_id:
        raise KeyError(
            "None of these Geom parameters were found: "
            f"parms={list(parm_names)}, groups={list(group_names)}"
        )
    vsp.SetParmVal(parm_id, float(value))

def get_xsec_parm_id(vsp, xsec_id: str, parm_names: Sequence[str]) -> tuple[str, str]:
    for name in parm_names:
        try:
            parm_id = vsp.GetXSecParm(xsec_id, name)
        except Exception:
            parm_id = ""
        if parm_id:
            return parm_id, name
    raise KeyError(f"None of these XSec parameters were found: {list(parm_names)}")

def get_xsec_value(vsp, xsec_id: str, parm_names: Sequence[str], default: float | None = None) -> float:
    try:
        parm_id, _ = get_xsec_parm_id(vsp, xsec_id, parm_names)
    except KeyError:
        if default is None:
            raise
        return float(default)
    return float(vsp.GetParmVal(parm_id))

def set_xsec_value(vsp, xsec_id: str, parm_names: Sequence[str], value: float) -> None:
    parm_id, _ = get_xsec_parm_id(vsp, xsec_id, parm_names)
    vsp.SetParmVal(parm_id, float(value))

def set_wing_section_driver_for_scaling(vsp, geom_id: str, section_index: int) -> None:
    try:
        vsp.SetDriverGroup(
            geom_id,
            section_index,
            vsp.SPAN_WSECT_DRIVER,
            vsp.ROOTC_WSECT_DRIVER,
            vsp.TIPC_WSECT_DRIVER,
        )
    except Exception:
        # Some OpenVSP builds may already have the desired driver group and may
        # object to setting it before all parameters are valid.  The subsequent
        # SetParmVal calls will still fail loudly if the section cannot be edited.
        pass

def set_analysis_input_if_available(vsp, analysis_name: str, input_names: set[str], setter, input_name: str, values) -> bool:
    if input_name in input_names:
        setter(analysis_name, input_name, values, 0)
        return True
    return False

def analysis_duration_seconds(vsp, result_id: str) -> float | None:
    if not result_id:
        return None
    try:
        values = vsp.GetDoubleResults(result_id, "Analysis_Duration_Sec", 0)
    except Exception:
        return None
    if not values:
        return None
    try:
        return float(values[0])
    except Exception:
        return None

def results_dataframe(vsp, result_id: str, result_name: str) -> pd.DataFrame:
    data, columns = [], []
    for child_result_id in vsp.GetStringResults(result_id, "ResultsVec"):
        if vsp.GetResultsName(child_result_id) != result_name:
            continue
        for data_name in vsp.GetAllDataNames(child_result_id):
            values = vsp.GetDoubleResults(child_result_id, data_name, 0)
            if values:
                data.append(values)
                columns.append(data_name)
        break
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(np.array(data).T, columns=columns)

def set_control_surface(vsp, geom_name, deflection, cs_group_name, sub_id=0, gains=(1, 1), verbose=1):
    """Enable and set an existing/missing VSPAERO control group, checking the readback.

    A saved group can be present but disabled (ActiveFlag=0). Setting its
    DeflectionAngle alone then has *no aerodynamic effect*. Do not try to add
    surfaces that already belong to the group: available-surface indices are
    one-based and refer only to GetAvailableCSNameVec, not GetActiveCSNameVec.
    """
    geom_id = find_one_geom(vsp, geom_name)
    subsurf_id = vsp.GetSubSurf(geom_id, sub_id)
    if not subsurf_id or str(subsurf_id).upper() == "NONE":
        raise RuntimeError(f"Missing control subsurface {geom_name}[{sub_id}]")

    group_names = [vsp.GetVSPAEROControlGroupName(i)
                   for i in range(vsp.GetNumControlSurfaceGroups())]
    if group_names.count(cs_group_name) > 1:
        raise RuntimeError(f"Duplicate control group: {cs_group_name}")
    if cs_group_name in group_names:
        group_index = group_names.index(cs_group_name)
    else:
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)
        vsp.Update()

    def belongs_to_target(name):
        return name.startswith(f"{geom_name}_") and name.endswith(
            f"_{vsp.GetSubSurfName(subsurf_id)}"
        ) if hasattr(vsp, "GetSubSurfName") else f"{geom_name}_" in name

    active = list(vsp.GetActiveCSNameVec(group_index))
    available = list(vsp.GetAvailableCSNameVec(group_index))
    selected = [i + 1 for i, name in enumerate(available)
                if belongs_to_target(name)]
    if selected:
        vsp.AddSelectedToCSGroup(selected, group_index)
        vsp.Update()
        active = list(vsp.GetActiveCSNameVec(group_index))
    expected = [name for name in active if belongs_to_target(name)]
    if len(expected) != len(gains):
        raise RuntimeError(
            f"{cs_group_name}: expected {len(gains)} active {geom_name} surfaces, "
            f"found {expected}; available={available}. "
            "Check saved SubSurfaces and control-group membership."
        )
    if len(active) != len(expected):
        raise RuntimeError(
            f"{cs_group_name}: group also contains unrelated control surfaces: {active}"
        )

    settings_id = find_vspaero_settings_container(vsp)
    if not settings_id:
        raise RuntimeError("VSPAERO settings container not found")
    group_name = f"ControlSurfaceGroup_{group_index}"
    active_id = vsp.FindParm(settings_id, "ActiveFlag", group_name)
    deflection_id = vsp.FindParm(settings_id, "DeflectionAngle", group_name)
    if not active_id or str(active_id).upper() == "NONE" or not deflection_id or str(deflection_id).upper() == "NONE":
        raise RuntimeError(f"{cs_group_name}: missing ActiveFlag/DeflectionAngle")

    vsp.SetParmVal(active_id, 1.0)
    for i, gain in enumerate(gains):
        parm_name = f"Surf_{subsurf_id}_{i}_Gain"
        gain_id = vsp.FindParm(settings_id, parm_name, group_name)
        if not gain_id or str(gain_id).upper() == "NONE":
            raise RuntimeError(f"{cs_group_name}: missing {parm_name}")
        vsp.SetParmVal(gain_id, float(gain))
    vsp.SetParmVal(deflection_id, float(deflection))
    vsp.Update()
    if not bool(vsp.GetParmVal(active_id)):
        raise RuntimeError(f"{cs_group_name}: ActiveFlag did not become enabled")
    if abs(vsp.GetParmVal(deflection_id) - float(deflection)) > 1e-8:
        raise RuntimeError(f"{cs_group_name}: DeflectionAngle readback mismatch")
    for i, gain in enumerate(gains):
        gain_id = vsp.FindParm(settings_id, f"Surf_{subsurf_id}_{i}_Gain", group_name)
        if abs(vsp.GetParmVal(gain_id) - float(gain)) > 1e-8:
            raise RuntimeError(f"{cs_group_name}: gain {i} readback mismatch")
    if verbose:
        print(f"{cs_group_name}: active={expected}, gains={tuple(gains)}, "
              f"deflection={deflection} deg")
    return deflection_id
