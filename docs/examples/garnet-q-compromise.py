"""Installed-package example: select each published garnet thermal EOS."""

from peritheos import get_material

for name in ("pyrope", "grossular"):
    record = get_material(name).get_eos_record(
        f"{name}_angel_2022_bm3_q_compromise_debye"
    )
    volume = record.volume(10.0, 1000.0)
    print(name, "V(10 GPa, 1000 K) =", volume, "angstrom^3/cell")
    assert abs(record.pressure(volume, 1000.0) - 10.0) < 1e-7
