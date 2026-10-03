"""Compare published Miozzi coefficients with pinned CrysFML/EosFit routines.

Requires gfortran and a local clone of
https://code.ill.fr/scientific-software/crysfml.git (passed as --source-checkout).
Extracted upstream routines are compiled in a temporary directory; this is a
calculation harness, not the EosFit console program or its refinement driver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from reproduce_miozzi_2020_iron import RECORDS, ROOT, observations, source_pressure

REVISION = "5bb3630d684387230065fa7c5933af2d81d19645"
SOURCES = {"eos": "Src/CFML_EoS_Mod.f90", "math": "Src/CFML_Math_Gen.f90"}
HEADER = """module audit
implicit none
integer,parameter:: cp=kind(1.d0),dp=kind(1.d0),n_eospar=24
logical:: err_eos=.false.,err_mathgen=.false.
character(200)::err_eos_mess,err_mathgen_mess
type eos_type
real(cp)::params(24)=0,tref=300
integer::itherm=7
logical::linear=.false.
character(15)::pscale_name='GPA',vscale_name='CM3/MOL'
end type
contains
subroutine eos_to_vec(e,ev)
type(eos_type),intent(in)::e
real(cp),intent(out)::ev(24)
ev=e%params
end subroutine
function u_case(a) result(b)
character(*),intent(in)::a
character(len(a))::b
b=a
end function
logical function vscaleMGD(e)
type(eos_type),intent(in)::e
vscaleMGD=index(e%vscale_name,'CM')>0
end function
real(dp) function debye(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
debye=debye3(x)
end function
"""
# The wrapper supports only uppercase GPa / cm3/mol, volume data, no
# transitions, and the third Debye function. Those are the settings audited.
# The cold formula below is the BM3 branch of upstream Get_Pressure, expressed
# directly in finite strain. Pthermal and its called functions are unmodified.
DRIVER = """
end module
program main
use audit
implicit none
type(eos_type)::e
real(cp)::v,t,f,p
integer::io
e%params(1)=22.81d0*6.02214076d23/2d24
e%params(2)=129.d0
e%params(3)=6.24d0
e%params(11)=420.d0
e%params(13)=1.d0
e%params(18)=1.11d0
e%params(19)=0.3d0
do
read(*,*,iostat=io) v,t
if(io/=0)exit
f=((e%params(1)/v)**(2.d0/3.d0)-1.d0)/2.d0
p=3.d0*e%params(2)*f*(1.d0+2.d0*f)**2.5d0*(1.d0+1.5d0*(e%params(3)-4.d0)*f)
write(*,'(3ES25.16)')p,pthermal(v,t,e),p+pthermal(v,t,e)
end do
end program
"""


def extract_function(source: str, name: str) -> str:
    match = re.search(
        rf"^\s*(?:pure )?function {name}\(.*?^\s*end function {name}\b",
        source,
        re.S | re.M | re.I,
    )
    if match is None:
        raise ValueError(f"Missing upstream routine: {name}")
    return match.group()


def audit(checkout: Path, compiler: str) -> dict:
    payloads = {
        key: subprocess.run(
            ["git", "show", f"{REVISION}:{path}"],
            cwd=checkout,
            capture_output=True,
            check=True,
        ).stdout
        for key, path in SOURCES.items()
    }
    sources = {key: value.decode("latin-1") for key, value in payloads.items()}
    program = HEADER + "\n".join(
        extract_function(sources["eos"], name)
        for name in ("EthDebye", "Get_DebyeT", "Get_Grun_V", "Pthermal")
    )
    program += "\n" + "\n".join(
        extract_function(sources["math"], name) for name in ("Cheval", "Debye3")
    )
    program += DRIVER
    data = observations()
    volume, temperature, pressure = (
        data[key] for key in ("volume_a3", "temperature_k", "pressure_gpa")
    )
    stdin = "".join(
        f"{v * 6.02214076e23 / 2e24:.17g} {t:.17g}\n"
        for v, t in zip(volume, temperature)
    )
    with tempfile.TemporaryDirectory(prefix="miozzi-eosfit-") as scratch:
        path = Path(scratch)
        (path / "audit.f90").write_text(program)
        subprocess.run(
            [compiler, "-O2", "-ffree-line-length-none", "-o", "audit", "audit.f90"],
            cwd=path,
            capture_output=True,
            check=True,
        )
        result = subprocess.run(
            [str(path / "audit")],
            input=stdin,
            text=True,
            capture_output=True,
            check=True,
        )
    calculated = np.loadtxt(result.stdout.splitlines())
    independent = source_pressure(
        volume, temperature, RECORDS["iron_miozzi_2020_bm3_mgd"], n=1
    )
    # Upstream uses default-real literals for R and the unit factor. Account
    # for their float32 representation before comparing only the equations.
    upstream_r = float(np.float32(8.314))
    upstream_factor = float(np.float32(1e-9)) * 1e6
    harmonized = calculated[:, 0] + calculated[:, 1] * (
        8.31446261815324 / upstream_r * 1e-3 / upstream_factor
    )
    np.testing.assert_allclose(harmonized, independent, rtol=0, atol=1e-9)

    with (ROOT / "peritheos/data/datasets/iron-miozzi-2020-mgo-pvt.csv").open() as f:
        rows = list(csv.DictReader(f))
    mgo_volume = np.array([float(row["mgo_volume_a3"]) for row in rows])
    cold = temperature[15:] == 300
    printed = pressure[15:][cold]
    standards = {
        "Speziale_2001_BM3": [74.71, 160.2, 3.99],
        "Tange_2009_Fit3_BM3": [74.698, 160.64, 4.221],
    }
    predictions = {
        key: source_pressure(mgo_volume[cold], 300, values)
        for key, values in standards.items()
    }
    return {
        "source_revision": REVISION,
        "source_files": {
            key: {"path": SOURCES[key], "sha256": hashlib.sha256(value).hexdigest()}
            for key, value in payloads.items()
        },
        "scope": "Unmodified thermal routines in a restricted calculation harness, not an EosFit console fit; double precision. All 131 original observations and published coefficients remain unchanged.",
        "thermal_equation_comparison": {
            "observations": len(volume),
            "max_difference_gpa": float(np.max(np.abs(calculated[:, 2] - independent))),
            "max_difference_after_R_harmonization_gpa": float(
                np.max(np.abs(harmonized - independent))
            ),
            "compiled_source_vs_printed_pressure_rmse_gpa": float(
                np.sqrt(np.mean((calculated[:, 2] - pressure) ** 2))
            ),
            "independent_vs_printed_pressure_rmse_gpa": float(
                np.sqrt(np.mean((independent - pressure) ** 2))
            ),
        },
        "cold_mgo_pressure_comparison": {
            "observations": int(np.sum(cold)),
            "qualification": "Room-temperature equations only. Closeness does not prove which scale generated the supplied column; no pressures are replaced.",
            "standards": {
                key: {
                    "parameters": standards[key],
                    "rmse_gpa": float(np.sqrt(np.mean((value - printed) ** 2))),
                }
                for key, value in predictions.items()
            },
            "rows": [
                {
                    "source_row": int(row["source_row"]),
                    "mgo_volume_a3": float(row["mgo_volume_a3"]),
                    "printed_pressure_gpa": float(row["pressure_gpa"]),
                    **{key: float(value[i]) for key, value in predictions.items()},
                }
                for i, row in enumerate(
                    [r for r in rows if float(r["temperature_k"]) == 300]
                )
            ],
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-checkout", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(args.source_checkout, args.compiler), indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")
