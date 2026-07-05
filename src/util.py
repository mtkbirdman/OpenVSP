"""Low-level OpenVSP helper functions shared by the analysis modules.

This module deliberately contains only thin OpenVSP API operations: import,
stdout/workdir handling, Geom/Parm lookup, XSec parameter access, simple
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


def get_container_parm_value(vsp, container_id: str, parm_name: str) -> tuple[float | None, str]:
    parm_id = find_container_parm(vsp, container_id, parm_name)
    if not parm_id:
        return None, ""
    return vsp.GetParmVal(parm_id), parm_id


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


def geom_x_location(vsp, geom_id: str) -> float:
    parm_id = vsp.FindParm(geom_id, "X_Rel_Location", "XForm")
    if not parm_id:
        parm_id = vsp.FindParm(geom_id, "X_Location", "XForm")
    if not parm_id:
        raise KeyError("Could not find X_Rel_Location/X_Location on Geom XForm.")
    return float(vsp.GetParmVal(parm_id))


def set_geom_x_location(vsp, geom_id: str, value: float) -> None:
    parm_id = vsp.FindParm(geom_id, "X_Rel_Location", "XForm")
    if not parm_id:
        parm_id = vsp.FindParm(geom_id, "X_Location", "XForm")
    if not parm_id:
        raise KeyError("Could not find X_Rel_Location/X_Location on Geom XForm.")
    vsp.SetParmVal(parm_id, float(value))


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
    current_group_names = [
        vsp.GetVSPAEROControlGroupName(group_index)
        for group_index in range(vsp.GetNumControlSurfaceGroups())
    ]
    if cs_group_name in current_group_names:
        group_index = current_group_names.index(cs_group_name)
    else:
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)

    cs_name_vec = vsp.GetAvailableCSNameVec(group_index)
    selected = [1 + i for i, cs_name in enumerate(cs_name_vec) if geom_name in cs_name]
    vsp.AddSelectedToCSGroup(selected, group_index)

    if verbose:
        print("\n", cs_group_name, ":", vsp.GetActiveCSNameVec(group_index))
        print("\t", f"{'container_id':14s}", f"{'parm_name':24s}", f"{'group_name':22s}", f"{'parm_id':14s}", f"{'value':8s}")

    geom_id = vsp.FindGeomsWithName(geom_name)[0]
    group_name = "ControlSurfaceGroup_" + str(group_index)
    cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)
    for i, gain in enumerate(gains):
        parm_name = "Surf_" + vsp.GetSubSurf(geom_id, sub_id) + "_" + str(i) + "_Gain"
        parm_id = find_container_parm(vsp, cs_group_container_id, parm_name)
        if parm_id != "":
            vsp.SetParmVal(parm_id, gain)
            if verbose:
                print("\t", f"{cs_group_container_id:14s}", f"{parm_name:24s}", f"{group_name:22s}", f"{parm_id:14s}", f"{gain:8.3f}")
        elif verbose:
            print("\t", "Failed to find " + parm_name)

    parm_name = "DeflectionAngle"
    parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name)
    vsp.SetParmVal(parm_id, deflection)
    vsp.Update()

    if deflection == vsp.GetParmVal(parm_id):
        if verbose:
            print("\t", f"{cs_group_container_id:14s}", f"{parm_name:24s}", f"{group_name:22s}", f"{parm_id:14s}", f"{deflection:8.3f}")
    else:
        if verbose:
            print("Failed to set deflection angle")
        exit()

    return parm_id
