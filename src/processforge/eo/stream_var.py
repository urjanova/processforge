"""StreamVar: variable container representing one process stream in the EO system."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass


@dataclass
class StreamVar:
    """
    Represents a single process stream in the global EO variable vector.

    Variable layout within the global x vector starting at global_offset:
        x[offset + 0]     = T  (temperature, K)
        x[offset + 1]     = P  (pressure, Pa)
        x[offset + 2]     = F  (molar flowrate, mol/s)
        x[offset + 3 + i] = z[comp_i]  (mole fraction, alphabetical order)

    All streams in a flowsheet share the same ordered ``components`` list,
    determined once at EOFlowsheet build time (alphabetical union of all
    components found in the JSON config).
    """

    name: str
    components: list[str]
    global_offset: int

    # Pyomo mode variables (set by PyomoBackend)
    pyo_T: Any = field(default=None, repr=False)
    pyo_P: Any = field(default=None, repr=False)
    pyo_F: Any = field(default=None, repr=False)
    pyo_z: Any = field(default=None, repr=False)  # IndexedVar over components

    # CasADi mode symbols (set by CasADiBackend)
    ca_T: Any = field(default=None, repr=False)
    ca_P: Any = field(default=None, repr=False)
    ca_F: Any = field(default=None, repr=False)
    ca_z: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def n_vars(self) -> int:
        """Number of scalar variables: T, P, F, z[0..N_c-1]."""
        return 3 + len(self.components)

    def to_vector_slice(self, x: np.ndarray) -> dict:
        """Extract this stream's values from the global x vector as a stream dict."""
        o = self.global_offset
        return {
            "T": float(x[o]),
            "P": float(x[o + 1]),
            "flowrate": float(x[o + 2]),
            "z": {c: float(x[o + 3 + i]) for i, c in enumerate(self.components)},
        }

    def from_stream_dict(self, stream: dict) -> np.ndarray:
        """Pack a stream dict into a flat numpy array segment (for initial guess)."""
        z = stream.get("z", {})
        return np.array(
            [
                stream.get("T", 298.15),
                stream.get("P", 101325.0),
                stream.get("flowrate", 1.0),
            ]
            + [z.get(c, 0.0) for c in self.components],
            dtype=float,
        )

    def write_to_vector(self, x: np.ndarray, stream: dict) -> None:
        """Write stream dict values into the global x vector at this stream's offset."""
        packed = self.from_stream_dict(stream)
        o = self.global_offset
        x[o: o + self.n_vars] = packed
