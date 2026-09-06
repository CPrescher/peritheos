"""Reproduce Caracas et al. (2005) CaSiO3 Table 2 EOS checks."""

from __future__ import annotations

from peritheos.eos.rt import BM3, BM4

MOLAR_MASS_G_MOL = 116.161
AVOGADRO = 6.02214076e23

PARAMETERS = {
    "Pm-3m": ((44.579, 250.0, 4.098), (44.588, 248.0, 4.206, -0.002)),
    "I4/mcm": ((44.537, 249.0, 4.090), (44.547, 247.0, 4.213, -0.002)),
    "Imma": ((44.567, 249.0, 4.094), (44.576, 247.0, 4.218, -0.002)),
    "R-3c": ((44.821, 247.0, 4.100), (44.832, 244.0, 4.236, -0.002)),
    "P4/mbm": ((44.629, 247.0, 4.124), (44.641, 244.0, 4.261, -0.002)),
    "I4/mmm": ((44.599, 250.0, 4.103), (44.609, 248.0, 4.219, -0.002)),
    "Im-3": ((44.600, 250.0, 4.104), (44.610, 247.0, 4.229, -0.002)),
    "P42/nmc": ((44.576, 248.0, 4.092), (44.566, 251.0, 3.977, -0.001)),
    "Pnma": ((44.576, 249.0, 4.104), (44.588, 246.0, 4.248, -0.002)),
}


def density_g_cm3(volume_a3_per_formula_unit: float) -> float:
    """Convert a one-formula-unit volume to density."""
    return MOLAR_MASS_G_MOL / (AVOGADRO * volume_a3_per_formula_unit * 1.0e-24)


def reproduce() -> dict[str, dict[str, float]]:
    """Return source-coefficient curve and independent density diagnostics."""
    metrics = {}
    for phase, (bm3_parameters, bm4_parameters) in PARAMETERS.items():
        bm3 = BM3(*bm3_parameters)
        bm4 = BM4(*bm4_parameters)
        volume_130 = float(bm3.volume(130.0))
        metrics[phase] = {
            "bm3_density_0_g_cm3": density_g_cm3(bm3.V0),
            "bm3_density_130_g_cm3": density_g_cm3(volume_130),
            "bm3_bm4_delta_160_gpa": abs(
                float(bm4.pressure(bm3.volume(160.0))) - 160.0
            ),
        }
    return metrics


if __name__ == "__main__":
    for phase, values in reproduce().items():
        print(
            f"{phase}: rho0={values['bm3_density_0_g_cm3']:.6f}, "
            f"rho130={values['bm3_density_130_g_cm3']:.6f} g/cm3, "
            f"BM3-BM4@BM3(160)={values['bm3_bm4_delta_160_gpa']:.6f} GPa"
        )
