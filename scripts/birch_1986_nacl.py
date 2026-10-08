"""Source-transcribed Birch NaCl ratio scale; conditional Walker normalization.

Birch (1986), Tables 5-6 (p. 4951), Equation 8 (p. 4952).
No absolute reference cell or author reduction spreadsheet is inferred.
"""

from scipy.optimize import brentq

SOURCE = {
    "doi": "10.1029/JB091iB05p04949",
    "pdf_sha256": "3b5b75419a3ceb70fc46e7c676b77cca2ad43dddb0449a4a5e3f2f9bc28c8e3a",
    "locations": [
        "Table 5, p. 4951",
        "Table 6, 25 Celsius extended row, p. 4951",
        "Equation 8, p. 4952",
    ],
    "K0_gpa": 23.88,
    "a": 1.796,
    "b": -5.00,
    "Tr_k": 298.15,
    "thermal_coefficient_gpa_per_k": 0.00286,
    "qualification": "Use the printed polynomial a,b, rather than reconstructing them from rounded elastic derivatives. Ratio reference is the zero-pressure 25 Celsius volume. Source construction: 25-500 Celsius, 0-30 GPa.",
}

WALKER_B2_ANCHORS = (
    {
        "nacl_file": "r34439",
        "temperature_k": 296.15,
        "nacl_lattice_a_angstrom": 5.6414,
        "nacl_lattice_a_esd_angstrom": 0.0014,
    },
    {
        "nacl_file": "r35101",
        "temperature_k": 297.15,
        "nacl_lattice_a_angstrom": 5.6468,
        "nacl_lattice_a_esd_angstrom": 0.0004,
    },
)


def pressure_ratio(ratio, temperature=298.15):
    """GPa from V/V00; Equation 8 with the final extended 25 C polynomial."""
    f = (ratio ** (-2.0 / 3.0) - 1.0) / 2.0
    return 3.0 * 23.88 * f * (1.0 + 2.0 * f) ** 2.5 * (
        1.0 + 1.796 * f - 5.0 * f**2
    ) + 0.00286 * (temperature - 298.15)


def normalized_pressure(lattice, temperature, anchor):
    """Normalize to zero pressure at the supplied anchor temperature.

    Callers must distinguish measured temperatures from conditional assumptions.
    """
    ambient_ratio = brentq(
        lambda ratio: pressure_ratio(ratio, anchor["temperature_k"]), 0.99, 1.01
    )
    ratio = (lattice / anchor["nacl_lattice_a_angstrom"]) ** 3 * ambient_ratio
    return pressure_ratio(ratio, temperature)


def table5_check():
    """Independent printed 25 C benchmarks, kbar converted to GPa."""
    return [
        {
            "volume_ratio": ratio,
            "printed_pressure_gpa": printed * 0.1,
            "calculated_pressure_gpa": pressure_ratio(ratio),
            "difference_gpa": pressure_ratio(ratio) - printed * 0.1,
        }
        for ratio, printed in [
            (1.0, 0.0),
            (0.95, 13.98),
            (0.90, 32.89),
            (0.85, 58.38),
            (0.80, 92.66),
            (0.75, 138.71),
            (0.70, 200.54),
            (0.65, 283.51),
        ]
    ]
