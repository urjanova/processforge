"""Typed models for Claims, Arguments and Evidence (CAE) cases.

A CAE case is a versioned YAML document that records a regulatory assurance
argument independently of any flowsheet. Cases reference flowsheet elements
(units, streams, artifacts) only through :class:`CAEEvidenceRef` values, which
are resolved against a concrete run at report time.

The schema is intentionally small and additive: new optional fields can be
introduced without breaking existing cases (``schema_version`` gates any
breaking change).
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, model_validator

#: How a raw result value is reduced before it is shown in a report.
ExtractOp = Literal["value", "final", "std_dev", "mean", "max", "min"]

#: The kind of flowsheet element an evidence item points at.
EvidenceRefType = Literal["unit", "stream", "artifact"]


class CAECaseMeta(BaseModel):
    """Identity and revision metadata for a CAE case."""

    id: str
    title: str
    version: str = "1.0.0"
    description: str = ""
    jurisdictions: list[str] = Field(default_factory=list)


class CAERequirement(BaseModel):
    """An external requirement the claims are traced to."""

    id: str
    statement: str = ""
    source: str = ""


class CAEClaim(BaseModel):
    """A regulatory claim, optionally scoped to a jurisdiction."""

    id: str
    statement: str
    jurisdiction: str = ""
    requirement_refs: list[str] = Field(default_factory=list)
    argument_refs: list[str] = Field(default_factory=list)
    summary: str = ""


class CAEArgument(BaseModel):
    """Reasoning that links a claim to one or more evidence items."""

    id: str
    statement: str
    evidence_refs: list[str] = Field(default_factory=list)


class CAEEvidenceRef(BaseModel):
    """A pointer at a flowsheet element whose result backs an argument.

    ``type`` selects the namespace (unit / stream / artifact), ``name`` is the
    element name as written in the flowsheet, ``field`` is the engine output
    field (required for units and streams), and ``extract`` selects how a
    timeseries/vector value is reduced. ``unit`` is the owning unit for an
    artifact reference; it is ignored for other types.
    """

    type: EvidenceRefType
    name: str
    field: str | None = None
    unit: str | None = None
    extract: ExtractOp = "value"

    @model_validator(mode="after")
    def _validate_ref(self) -> "CAEEvidenceRef":
        if self.type in ("unit", "stream") and not self.field:
            raise ValueError(f"evidence ref of type '{self.type}' requires a field")
        if self.type == "artifact" and self.field:
            raise ValueError("artifact evidence refs must not declare a field")
        return self


class CAEEvidence(BaseModel):
    """A named piece of evidence backed by a flowsheet result."""

    id: str
    description: str
    ref: CAEEvidenceRef


class CAECase(BaseModel):
    """A complete, validated CAE case."""

    schema_version: str = "1.0"
    case: CAECaseMeta
    requirements: list[CAERequirement] = Field(default_factory=list)
    claims: list[CAEClaim] = Field(default_factory=list)
    arguments: list[CAEArgument] = Field(default_factory=list)
    evidence: list[CAEEvidence] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_references(self) -> "CAECase":
        _require_unique("requirement", [r.id for r in self.requirements])
        _require_unique("claim", [c.id for c in self.claims])
        _require_unique("argument", [a.id for a in self.arguments])
        _require_unique("evidence", [e.id for e in self.evidence])

        requirement_ids = {r.id for r in self.requirements}
        argument_ids = {a.id for a in self.arguments}
        evidence_ids = {e.id for e in self.evidence}

        for claim in self.claims:
            missing_req = set(claim.requirement_refs) - requirement_ids
            if missing_req:
                raise ValueError(
                    f"claim '{claim.id}' references unknown requirement(s): "
                    f"{sorted(missing_req)}"
                )
            missing_arg = set(claim.argument_refs) - argument_ids
            if missing_arg:
                raise ValueError(
                    f"claim '{claim.id}' references unknown argument(s): "
                    f"{sorted(missing_arg)}"
                )

        for argument in self.arguments:
            missing_evd = set(argument.evidence_refs) - evidence_ids
            if missing_evd:
                raise ValueError(
                    f"argument '{argument.id}' references unknown evidence: "
                    f"{sorted(missing_evd)}"
                )
        return self

    def jurisdictions(self) -> list[str]:
        """Return every jurisdiction declared by the case or its claims."""
        seen: dict[str, None] = {}
        for jurisdiction in self.case.jurisdictions:
            seen.setdefault(jurisdiction, None)
        for claim in self.claims:
            if claim.jurisdiction:
                seen.setdefault(claim.jurisdiction, None)
        return list(seen)

    def filter_jurisdiction(self, jurisdiction: str) -> "CAECase":
        """Return a copy containing only claims for *jurisdiction*.

        Arguments and evidence are pruned to those reachable from the retained
        claims, so the report never shows orphaned reasoning.
        """
        claims = [
            c
            for c in self.claims
            if not c.jurisdiction or c.jurisdiction == jurisdiction
        ]
        argument_ids = {ref for c in claims for ref in c.argument_refs}
        arguments = [a for a in self.arguments if a.id in argument_ids]
        evidence_ids = {ref for a in arguments for ref in a.evidence_refs}
        evidence = [e for e in self.evidence if e.id in evidence_ids]
        requirement_ids = {ref for c in claims for ref in c.requirement_refs}
        requirements = [r for r in self.requirements if r.id in requirement_ids]
        return self.model_copy(
            update={
                "claims": claims,
                "arguments": arguments,
                "evidence": evidence,
                "requirements": requirements,
            }
        )


class CAECaseLink(BaseModel):
    """A flowsheet-to-case association recorded in ``cae/index.yaml``."""

    case_id: str
    version: str | None = None


class CAEFlowsheetLink(BaseModel):
    """All CAE cases linked to one flowsheet path."""

    path: str
    cases: list[CAECaseLink] = Field(default_factory=list)


class CAEIndex(BaseModel):
    """The repository-wide flowsheet-to-case index."""

    schema_version: str = "1.0"
    flowsheets: list[CAEFlowsheetLink] = Field(default_factory=list)


class CAECaseRef(BaseModel):
    """Lightweight description of a case available in a repository."""

    case_id: str
    version: str
    title: str = ""
    jurisdictions: list[str] = Field(default_factory=list)
    path: str = ""


def _require_unique(kind: str, ids: list[str]) -> None:
    duplicates = sorted({i for i in ids if ids.count(i) > 1})
    if duplicates:
        raise ValueError(f"duplicate {kind} id(s): {duplicates}")
