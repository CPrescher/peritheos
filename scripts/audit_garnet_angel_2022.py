"""Audit Angel (2022) source q-compromise parameterizations, without refitting.

Run ``python -m scripts.audit_garnet_angel_2022``. Source files are the unchanged
CC-BY-4.0 official supplements. Independent SciPy quadrature implements the
printed law; it does not call Peritheos. The complete original mixed V/Ks fit
inputs, selections, scale-factor protocol and weights are not distributed by
the source supplements. This script does not invent a joint-refit result.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "docs/data/garnet-angel-2022"
NAMES = ("pyrope", "grossular")


def source(name):
    """Read all 59 parameter tokens and the complete source covariance."""
    path = FOLDER / "sources" / f"{name}-official.eos"
    text = path.read_text(encoding="utf-8")
    parameters = {
        int(k): float(v) for k, v in re.findall(r"^Param =\s*(\d+)\s+(\S+)", text, re.M)
    }
    rows = [
        re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)
        for line in text.split("Variance-Covariance matrix=")[1].splitlines()
    ]
    covariance = np.array([[float(x) for x in r] for r in rows if len(r) == 59])
    if len(parameters) != 59 or covariance.shape != (59, 59):
        raise ValueError("Incomplete source parameter/covariance block")
    return parameters, covariance


def independent_pressure(parameters, volume_cm3_mol, temperature, gas_constant=R):
    """Independent SI/molar evaluation of BM3 + constant-theta Debye pressure."""
    p = parameters
    v0, k0, kp, theta, n, gamma = (p[i] for i in (1, 2, 3, 11, 13, 18))
    eta = (v0 / volume_cm3_mol) ** (1 / 3)
    cold = 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))

    def energy(t):
        y = theta / t
        integral = quad(
            lambda x: x**3 / np.expm1(x) if x else 0.0,
            0.0,
            y,
            epsabs=1e-11,
            epsrel=1e-12,
        )[0]
        return 9 * n * gas_constant * t * integral / y**3

    return cold + gamma * (energy(temperature) - energy(298.0)) / v0 / 1000


def reproduce():
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import QCompromiseDebye

    report = {
        "scope": "source_parameterization_only",
        "author_workflow_reproduced": False,
        "source_doi": "10.1007/s00410-022-01918-5",
        "records": {},
    }
    for name in NAMES:
        p, covariance = source(name)
        model = QCompromiseDebye(BM3(p[1] / 10, p[2], p[3]), 298.0, p[11], p[18], p[13])
        checks = []
        for ratio in (0.85, 0.95, 1.0, 1.03):
            for t in (90.0, 298.0, 700.0, 1000.0, 1700.0):
                v = p[1] * ratio
                expected = independent_pressure(p, v, t)
                observed = float(model.pressure(v / 10, t))
                checks.append(
                    {
                        "volume_cm3_mol": v,
                        "temperature_k": t,
                        "independent_pressure_gpa": expected,
                        "native_pressure_gpa": observed,
                        "difference_gpa": observed - expected,
                    }
                )
        indices = [1, 2, 10, 17]
        c = covariance[np.ix_(indices, indices)]
        np.testing.assert_allclose(c, c.T, atol=0, rtol=0)
        if np.linalg.eigvalsh(c)[0] <= 0:
            raise ValueError("Nonpositive EOS covariance")
        h = 0.001
        alpha = (
            independent_pressure(p, p[1], 298 + h)
            - independent_pressure(p, p[1], 298 - h)
        ) / (2 * h * p[2])
        report["records"][name] = {
            "source_file": f"sources/{name}-official.eos",
            "source_sha256": hashlib.sha256(
                (FOLDER / "sources" / f"{name}-official.eos").read_bytes()
            ).hexdigest(),
            "source_parameters": {str(i): p[i] for i in (1, 2, 3, 11, 13, 14, 18)},
            "volume0_cell_a3": p[1] * 8e24 / Avogadro,
            "covariance_parameter_order": [
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "theta0",
                "gamma0",
            ],
            "covariance": c.tolist(),
            "covariance_standard_errors": np.sqrt(c.diagonal()).tolist(),
            "alpha_298k_inverse_k": alpha,
            "table1_alpha_inverse_k": {"pyrope": 2.256e-5, "grossular": 2.047e-5}[name],
            "maximum_independent_difference_gpa": max(
                abs(x["difference_gpa"]) for x in checks
            ),
            "calculated_checks_not_observations": checks,
            "unresolved_fit_inputs": "Complete joint V/Ks input project, exact point masks, rescalings, weights (including increased Hartwig weights), and fit covariance confidence convention.",
        }
    return report


def ledger_outcome(record):
    """A source-equation check is not coefficient refitting or author workflow."""
    name = record["identifier"].split("_")[0]
    result = reproduce()["records"][name]
    return {
        "status": "not_refittable",
        "dataset_identifiers": [],
        "fit_kind": "source_curve_validation_only",
        "reason": result["unresolved_fit_inputs"],
        "assessment_basis": "Published parameterization verified by independent quadrature and official EosFit checks; complete original joint fit is unavailable.",
        "qualification": "Covariance is source-reported and conditional; no original-fit coefficient or formal-uncertainty reproduction is claimed.",
        "source_equation_maximum_difference_gpa": result[
            "maximum_independent_difference_gpa"
        ],
    }


if __name__ == "__main__":
    payload = reproduce()
    (FOLDER / "audit.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    for name, record in payload["records"].items():
        print(
            name,
            record["maximum_independent_difference_gpa"],
            record["alpha_298k_inverse_k"],
        )
