"""Read VSPAERO .stab files and resolve their control-derivative columns.

This module owns the .stab file format.  Flight-mechanics solvers and chart
post-processing import the parsed data instead of depending on one another.
Angles exposed by :class:`VSPAEROStab` are radians; scalar values otherwise
remain in the unit system written by VSPAERO.
"""
from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from typing import Mapping

import pandas as pd

BASE_AERO_COLUMNS = [
    "CFx", "CFy", "CFz", "CMx", "CMy", "CMz",
    "CL", "CD", "CS", "CMl", "CMm", "CMn",
]
STABILITY_CASE_NAMES = {
    "Base_Aero", "Alpha", "Beta", "Roll__Rate", "Pitch_Rate", "Yaw___Rate", "Mach",
}
FALLBACK_DERIVATIVE_COLUMNS = [
    "Alpha", "Beta", "p", "q", "r", "Mach", "U",
    "ConGrp_1", "ConGrp_2", "ConGrp_3",
]
CONTROL_NAME_HINTS = {
    "delta_a": "AILERON",
    "delta_e": "ELEVATOR",
    "delta_r": "RUDDER",
}
CONTROL_COLUMNS_FALLBACK = {
    "delta_a": "ConGrp_1",
    "delta_e": "ConGrp_2",
    "delta_r": "ConGrp_3",
}
LINEAR_AERO_COEFFICIENTS = ("CL", "CD", "CS", "CMl", "CMm", "CMn")
RAW_FORCE_COEFFICIENTS = ("CFx", "CFy", "CFz")

@dataclass
class VSPAEROStab:
    """Parsed values from one VSPAERO .stab file."""

    path: str
    references: dict
    base_condition: dict
    base_aero: dict
    derivatives: pd.DataFrame
    control_groups: dict

    @property
    def Sref(self) -> float:
        return float(self.references["Sref"])

    @property
    def Cref(self) -> float:
        return float(self.references["Cref"])

    @property
    def Bref(self) -> float:
        return float(self.references["Bref"])

    @property
    def alpha0(self) -> float:
        return math.radians(float(self.base_condition.get("AoA", 0.0)))

    @property
    def beta0(self) -> float:
        return math.radians(float(self.base_condition.get("Beta", 0.0)))

    @property
    def V0(self) -> float:
        return float(self.base_condition.get("Vinf", 0.0) or 0.0)

    @property
    def rho0(self) -> float | None:
        value = self.base_condition.get("Rho")
        return None if value is None else float(value)

def read_vspaero_stab(stab_path: str | os.PathLike) -> VSPAEROStab:
    """Read the scalar header, base case, control groups, and derivative table."""

    stab_path = os.fspath(stab_path)
    with open(stab_path, "r", encoding="utf-8", errors="ignore") as file:
        lines = file.readlines()

    scalar_values: dict[str, float] = {}
    for line in lines:
        if line.lstrip().startswith("Case"):
            break
        match = re.match(
            r"^\s*([A-Za-z][A-Za-z0-9_]*_?)\s+([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\b",
            line,
        )
        if match:
            name = match.group(1).strip()
            if name.endswith("_"):
                name = name[:-1]
            scalar_values[name] = float(match.group(2))

    references = {
        "Sref": scalar_values["Sref"],
        "Cref": scalar_values["Cref"],
        "Bref": scalar_values["Bref"],
        "Xcg": scalar_values.get("Xcg"),
        "Ycg": scalar_values.get("Ycg"),
        "Zcg": scalar_values.get("Zcg"),
    }
    base_condition = {
        "Mach": scalar_values.get("Mach"),
        "AoA": scalar_values.get("AoA", 0.0),
        "Beta": scalar_values.get("Beta", 0.0),
        "Rho": scalar_values.get("Rho"),
        "Vinf": scalar_values.get("Vinf", 0.0),
    }

    base_aero = None
    derivative_header_index = None
    control_group_names: list[str] = []
    for index, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith("Base_Aero"):
            parts = line.split()
            base_aero = {
                name: float(value)
                for name, value in zip(BASE_AERO_COLUMNS, parts[3:15])
            }

        parts = line.split()
        if len(parts) >= 15 and parts[0] not in STABILITY_CASE_NAMES and not parts[0].startswith("#"):
            try:
                float(parts[1])
            except ValueError:
                pass
            else:
                control_group_names.append(parts[0])

        if stripped.startswith("Coef") and "Alpha" in line and "ConGrp_1" in line:
            derivative_header_index = index
            break

    if base_aero is None:
        raise ValueError("Base_Aero row was not found in the .stab file.")
    if derivative_header_index is None:
        raise ValueError("Derivative table was not found in the .stab file.")

    header_parts = lines[derivative_header_index].split()
    derivative_column_names = header_parts[2:] or FALLBACK_DERIVATIVE_COLUMNS
    derivative_rows = []
    derivative_start = derivative_header_index + 4
    for line in lines[derivative_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < len(derivative_column_names) + 2:
            break
        coefficient = parts[0]
        if coefficient not in BASE_AERO_COLUMNS:
            break
        values = [float(value) for value in parts[1:len(derivative_column_names) + 2]]
        derivative_rows.append([coefficient] + values)

    columns = ["Coef", "Base"] + derivative_column_names
    derivatives = pd.DataFrame(derivative_rows, columns=columns).set_index("Coef")
    control_groups = {
        f"ConGrp_{index + 1}": name
        for index, name in enumerate(control_group_names)
    }
    return VSPAEROStab(
        path=stab_path,
        references=references,
        base_condition=base_condition,
        base_aero=base_aero,
        derivatives=derivatives,
        control_groups=control_groups,
    )

def resolve_control_columns_from_stab(
    stab: VSPAEROStab,
    control_map: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Map delta_a/e/r to .stab ConGrp_* derivative columns."""

    available_columns = set(stab.derivatives.columns)
    group_to_column = {
        group_name: column
        for column, group_name in stab.control_groups.items()
    }

    if control_map:
        result: dict[str, str] = {}
        for delta_name, target in control_map.items():
            if target in available_columns:
                result[delta_name] = target
            elif target in group_to_column:
                result[delta_name] = group_to_column[target]
        return result

    result = {}
    for delta_name, hint in CONTROL_NAME_HINTS.items():
        for column, group_name in stab.control_groups.items():
            if hint in group_name.upper():
                result[delta_name] = column
                break

    for delta_name, fallback_column in CONTROL_COLUMNS_FALLBACK.items():
        if delta_name not in result and fallback_column in available_columns:
            result[delta_name] = fallback_column
    return result

def stab_coefficient_value(stab: VSPAEROStab, coefficient: str, column: str) -> float:
    """Return one raw value from the parsed .stab coefficient table."""

    if coefficient not in stab.derivatives.index:
        raise KeyError(f"Coefficient row {coefficient!r} was not found in {stab.path}.")
    if column not in stab.derivatives.columns:
        raise KeyError(f"Derivative column {column!r} was not found in {stab.path}.")
    return float(stab.derivatives.loc[coefficient, column])

def evaluate_stab_linear_aero(
    stab: VSPAEROStab,
    increments: Mapping[str, float],
    *,
    alpha: float,
    beta: float,
    include_raw_force_coefficients: bool = False,
) -> dict[str, float]:
    """Evaluate .stab linear coefficients and convert forces to solver axes."""

    coefficients: dict[str, float] = {}
    for coefficient in LINEAR_AERO_COEFFICIENTS:
        row = stab.derivatives.loc[coefficient]
        value = float(row["Base"])
        for derivative_name, increment in increments.items():
            if derivative_name in row.index:
                value += float(row[derivative_name]) * float(increment)
        coefficients[coefficient] = value

    if include_raw_force_coefficients:
        for coefficient in RAW_FORCE_COEFFICIENTS:
            if coefficient not in stab.derivatives.index:
                continue
            row = stab.derivatives.loc[coefficient]
            value = float(row["Base"])
            for derivative_name, increment in increments.items():
                if derivative_name in row.index:
                    value += float(row[derivative_name]) * float(increment)
            coefficients[f"{coefficient}_raw"] = value

    CL = coefficients["CL"]
    CD = coefficients["CD"]
    CS = coefficients["CS"]
    axial_opposite_forward = CD * math.cos(beta) + CS * math.sin(beta)
    coefficients["CX"] = -axial_opposite_forward * math.cos(alpha) + CL * math.sin(alpha)
    coefficients["CY"] = -CD * math.sin(beta) + CS * math.cos(beta)
    coefficients["CZ"] = -axial_opposite_forward * math.sin(alpha) - CL * math.cos(alpha)
    coefficients["Cl"] = coefficients["CMl"]
    coefficients["Cm"] = coefficients["CMm"]
    coefficients["Cn"] = coefficients["CMn"]
    return coefficients

def solver_axis_derivative_value(
    stab: VSPAEROStab,
    coefficient: str,
    column: str,
) -> float:
    """Return an exact first derivative in the solver-axis convention."""

    moment_aliases = {"Cl": "CMl", "Cm": "CMm", "Cn": "CMn"}
    if coefficient in moment_aliases:
        return stab_coefficient_value(stab, moment_aliases[coefficient], column)
    if coefficient in LINEAR_AERO_COEFFICIENTS + RAW_FORCE_COEFFICIENTS:
        return stab_coefficient_value(stab, coefficient, column)
    if coefficient not in {"CX", "CY", "CZ"}:
        raise KeyError(f"Unsupported solver-axis coefficient {coefficient!r}.")

    alpha = stab.alpha0
    beta = stab.beta0
    CL = stab_coefficient_value(stab, "CL", "Base")
    CD = stab_coefficient_value(stab, "CD", "Base")
    CS = stab_coefficient_value(stab, "CS", "Base")
    dCL = stab_coefficient_value(stab, "CL", column)
    dCD = stab_coefficient_value(stab, "CD", column)
    dCS = stab_coefficient_value(stab, "CS", column)
    d_alpha = 1.0 if column == "Alpha" else 0.0
    d_beta = 1.0 if column == "Beta" else 0.0

    axial_opposite_forward = CD * math.cos(beta) + CS * math.sin(beta)
    d_axial_opposite_forward = (
        dCD * math.cos(beta)
        - CD * math.sin(beta) * d_beta
        + dCS * math.sin(beta)
        + CS * math.cos(beta) * d_beta
    )

    if coefficient == "CX":
        return float(
            -d_axial_opposite_forward * math.cos(alpha)
            + axial_opposite_forward * math.sin(alpha) * d_alpha
            + dCL * math.sin(alpha)
            + CL * math.cos(alpha) * d_alpha
        )
    if coefficient == "CY":
        return float(
            -dCD * math.sin(beta)
            - CD * math.cos(beta) * d_beta
            + dCS * math.cos(beta)
            - CS * math.sin(beta) * d_beta
        )
    return float(
        -d_axial_opposite_forward * math.sin(alpha)
        - axial_opposite_forward * math.cos(alpha) * d_alpha
        - dCL * math.cos(alpha)
        + CL * math.sin(alpha) * d_alpha
    )
