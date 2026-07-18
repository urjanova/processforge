from pydantic import BaseModel


class StreamSchema(BaseModel):
    variables: list[str]
    dtypes: dict[str, str]
    units: dict[str, str] = {}
    shape: list[int]
    has_time: bool = False
    has_phase: bool = False


class SolverUnitSchema(BaseModel):
    variables: list[str]
    dtypes: dict[str, str]


class ResultSchema(BaseModel):
    version: int = 1
    store_type: str = "simulation_results"
    created: str
    mode: str
    processforge_version: str = ""
    streams: dict[str, StreamSchema] = {}
    solver_units: dict[str, SolverUnitSchema] = {}
    provenance: dict = {}
