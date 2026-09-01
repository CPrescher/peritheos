"""Peritheos: thermodynamic equations of state calculations."""

from peritheos.eosmat import (
    EOSMAT_FORMAT,
    EOSMAT_FORMAT_VERSION,
    eosmat_schema,
    get_material_document,
    list_material_documents,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
from peritheos.materials import (
    EOSRecord,
    Material,
    get_eos_record,
    get_material,
    list_eos_records,
    list_materials,
    material_from_dict,
)
from peritheos.uncertainty import (
    EOSUncertainty,
    ParameterUncertainty,
    PredictionUncertainty,
)

__version__ = "0.6.0"

__all__ = [
    "EOSMAT_FORMAT",
    "EOSMAT_FORMAT_VERSION",
    "EOSUncertainty",
    "ParameterUncertainty",
    "EOSRecord",
    "Material",
    "PredictionUncertainty",
    "__version__",
    "eosmat_schema",
    "get_eos_record",
    "get_material",
    "get_material_document",
    "list_eos_records",
    "list_materials",
    "list_material_documents",
    "load_eosmat",
    "material_from_dict",
    "save_eosmat",
    "validate_eosmat_document",
]
