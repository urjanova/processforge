"""Claims, Arguments and Evidence (CAE) reporting.

CAE data lives outside the flowsheet schema and is linked to flowsheets through
``cae/index.yaml``. Use :func:`generate_cae_reports` (or ``pf runs ... --cae``)
to render an assurance report for a concrete simulation run.
"""

from .evidence import ResolvedEvidence, resolve_case_evidence
from .file_repository import FileCAERepository
from .loader import CAEError, find_cae_root, load_case, load_index
from .repository import CAERepository
from .service import CAEReport, generate_cae_reports
from .types import (
    CAEArgument,
    CAECase,
    CAECaseLink,
    CAECaseMeta,
    CAECaseRef,
    CAEClaim,
    CAEEvidence,
    CAEEvidenceRef,
    CAEFlowsheetLink,
    CAEIndex,
    CAERequirement,
)

__all__ = [
    "CAEArgument",
    "CAECase",
    "CAECaseLink",
    "CAECaseMeta",
    "CAECaseRef",
    "CAEClaim",
    "CAEError",
    "CAEEvidence",
    "CAEEvidenceRef",
    "CAEFlowsheetLink",
    "CAEIndex",
    "CAEReport",
    "CAERequirement",
    "CAERepository",
    "FileCAERepository",
    "ResolvedEvidence",
    "find_cae_root",
    "generate_cae_reports",
    "load_case",
    "load_index",
    "resolve_case_evidence",
]
