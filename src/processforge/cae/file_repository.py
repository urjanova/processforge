"""Local filesystem-backed :class:`CAERepository` implementation."""

from __future__ import annotations

import os

from .loader import (
    CASES_DIRNAME,
    CAEError,
    iter_case_paths,
    load_case,
    load_index,
    select_version,
    version_key,
)
from .repository import CAERepository
from .types import CAECase, CAECaseLink, CAECaseRef


class FileCAERepository(CAERepository):
    """Read CAE cases from a local ``cae/`` directory.

    Layout::

        cae/
        ├── index.yaml
        └── cases/
            ├── closed-loop-cae.ca.yaml
            └── hydraulic-cae.ca.yaml
    """

    def __init__(self, root: str) -> None:
        self.root = os.path.abspath(os.path.expanduser(root))
        self._cases_dir = os.path.join(self.root, CASES_DIRNAME)
        self._cases: dict[tuple[str, str], CAECase] | None = None
        self._index = None

    # ------------------------------------------------------------------
    # Index / case loading
    # ------------------------------------------------------------------
    @property
    def index(self):
        if self._index is None:
            self._index = load_index(self.root)
        return self._index

    def _load_all(self) -> dict[tuple[str, str], CAECase]:
        if self._cases is not None:
            return self._cases

        cases: dict[tuple[str, str], CAECase] = {}
        paths: dict[tuple[str, str], str] = {}
        for path in iter_case_paths(self._cases_dir):
            case = load_case(path)
            key = (case.case.id, case.case.version)
            if key in cases:
                raise CAEError(
                    f"duplicate CAE case '{case.case.id}' version "
                    f"'{case.case.version}' in '{paths[key]}' and '{path}'"
                )
            cases[key] = case
            paths[key] = path
        self._cases = cases
        return cases

    def list_cases(self) -> list[CAECaseRef]:
        refs = [
            CAECaseRef(
                case_id=case.case.id,
                version=case.case.version,
                title=case.case.title,
                jurisdictions=case.jurisdictions(),
            )
            for case in self._load_all().values()
        ]
        return sorted(refs, key=lambda ref: (ref.case_id, version_key(ref.version)))

    def load_case(self, case_id: str, version: str | None = None) -> CAECase:
        cases = self._load_all()
        versions = [v for (cid, v) in cases if cid == case_id]
        chosen = select_version(versions, version)
        if chosen is None:
            available = ", ".join(sorted(versions, key=version_key)) or "none"
            if version:
                raise CAEError(
                    f"CAE case '{case_id}' has no version '{version}' "
                    f"(available: {available})"
                )
            raise CAEError(f"CAE case '{case_id}' not found (available: {available})")
        return cases[(case_id, chosen)]

    # ------------------------------------------------------------------
    # Flowsheet links
    # ------------------------------------------------------------------
    def links_for_flowsheet(self, flowsheet_path: str) -> list[CAECaseLink]:
        best_score = 0
        best_links: list[CAECaseLink] = []
        matches: list[str] = []

        for entry in self.index.flowsheets:
            score = self._match_score(entry.path, flowsheet_path)
            if score == 0:
                continue
            if score > best_score:
                best_score = score
                best_links = list(entry.cases)
                matches = [entry.path]
            elif score == best_score:
                matches.append(entry.path)

        if best_score == 0:
            raise CAEError(
                f"No CAE case linked to flowsheet '{flowsheet_path}'. "
                f"Add it to '{os.path.join(self.root, 'index.yaml')}'."
            )
        if len(matches) > 1:
            raise CAEError(
                f"Ambiguous CAE index entries for flowsheet '{flowsheet_path}': "
                f"{matches}. Use a more specific path."
            )
        return best_links

    def _match_score(self, indexed_path: str, flowsheet_path: str) -> int:
        repo_root = os.path.dirname(self.root)
        indexed_abs = _abspath(indexed_path, repo_root)
        given_abs = os.path.abspath(os.path.expanduser(flowsheet_path))

        if indexed_abs == given_abs:
            return 100

        indexed_rel = os.path.normpath(indexed_path)
        given_rel = os.path.normpath(flowsheet_path)
        if indexed_rel == given_rel:
            return 90

        indexed_name = os.path.splitext(os.path.basename(indexed_path))[0]
        given_name = os.path.splitext(os.path.basename(flowsheet_path))[0]
        if indexed_name != given_name:
            return 0

        # Same basename: prefer a match on the trailing path components.
        if indexed_rel.endswith(given_rel) or given_rel.endswith(indexed_rel):
            return 80
        indexed_parent = os.path.basename(os.path.dirname(indexed_rel))
        given_parent = os.path.basename(os.path.dirname(given_rel))
        if indexed_parent and indexed_parent == given_parent:
            return 60
        return 40


def _abspath(path: str, base: str) -> str:
    expanded = os.path.expanduser(path)
    if os.path.isabs(expanded):
        return os.path.abspath(expanded)
    return os.path.abspath(os.path.join(base, expanded))
