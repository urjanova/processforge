"""ModelicaJacobianMixin: Jacobian contributions via FMI 2.0 directional derivatives.

Implements the ``JacobianContributor`` protocol for ``ModelicaProvider`` by
subclassing ``BaseJacobianMixin`` and overriding ``compute_jacobian_block``.

FMI 2.0 directional derivatives
--------------------------------
FMI 2.0 Model Exchange (ME) FMUs may expose::

    fmu.getDirectionalDerivative(v_unknown_refs, v_known_refs, dv_known)

which computes the Jacobian-vector product J·dv for a single direction dv.
To extract column i of the Jacobian, call with ``dv = e_i`` (unit vector).

**Important constraint:** FMI 2.0 Co-Simulation (CS) slaves — which
``ModelicaProvider`` uses for time-stepping via ``FMU2Slave`` — do *not*
mandate ``getDirectionalDerivative``.  Only ME FMUs with the flag::

    providesDirectionalDerivatives = true

in ``modelDescription.xml`` support it.

This mixin creates a *second* FMU instance (``FMU2Model``, ME type) solely for
Jacobian evaluation, separate from the CS slave used by ``compute_unit()``.
If the compiled FMU does not include a ME interface or the flag is False,
``_jac_fmu`` is set to ``None`` and ``compute_jacobian_block`` returns ``[]``
for all unit types — the global FD path is used transparently.

Usage
-----
Added to ``ModelicaProvider`` by extending its MRO::

    class ModelicaProvider(AbstractProvider, ModelicaJacobianMixin): ...

And calling ``_initialize_jacobian_fmu`` at the end of ``initialize()``::

    self._initialize_jacobian_fmu(fmu_path, self._model_desc)
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from loguru import logger

from .base_jacobian_mixin import BaseJacobianMixin

if TYPE_CHECKING:
    from .jacobian_contributor import ReferenceState


class ModelicaJacobianMixin(BaseJacobianMixin):
    """Jacobian mixin for ModelicaProvider using FMI 2.0 directional derivatives.

    If the FMU supports directional derivatives (ME interface +
    ``providesDirectionalDerivatives=true``), uses ``getDirectionalDerivative``
    to compute exact Jacobian columns.

    Otherwise falls back to ``[]`` for all unit types (global FD used).
    """

    # ------------------------------------------------------------------
    # Initialization (called from ModelicaProvider.initialize())
    # ------------------------------------------------------------------

    def _initialize_jacobian_fmu(self, fmu_path: str, model_desc) -> None:
        """Attempt to set up a Model Exchange FMU instance for directional derivatives.

        Must be called at the end of ``ModelicaProvider.initialize()`` after
        ``self._model_desc``, ``self._fmu_dir``, and ``self._var_refs`` are set.

        Sets ``self._jac_fmu`` to an initialized ``FMU2Model`` instance if the
        FMU supports directional derivatives, otherwise sets it to ``None``.

        Args:
            fmu_path:    Path to the ``.fmu`` archive.
            model_desc:  Parsed model description from ``fmpy.read_model_description``.
        """
        self._jac_fmu = None
        self._jac_var_refs: dict[str, int] = {}

        # Check if the FMU exposes a Model Exchange interface at all.
        me = getattr(model_desc, "modelExchange", None)
        if me is None:
            logger.debug(
                "ModelicaJacobianMixin: FMU has no Model Exchange interface; "
                "using global FD for Jacobian."
            )
            return

        provides_dd = getattr(me, "providesDirectionalDerivatives", False)
        if not provides_dd:
            logger.debug(
                "ModelicaJacobianMixin: FMU does not set "
                "providesDirectionalDerivatives=true; using global FD."
            )
            return

        try:
            from fmpy.fmi2 import FMU2Model
        except ImportError:
            logger.debug(
                "ModelicaJacobianMixin: fmpy not available; using global FD."
            )
            return

        try:
            jac_fmu = FMU2Model(
                guid=model_desc.guid,
                unzipDirectory=self._fmu_dir,
                modelIdentifier=me.modelIdentifier,
                instanceName="processforge_jac",
            )
            jac_fmu.instantiate()
            jac_fmu.setupExperiment(startTime=0.0)
            jac_fmu.enterInitializationMode()
            jac_fmu.exitInitializationMode()
            jac_fmu.enterEventMode()
            jac_fmu.newDiscreteStates()
            jac_fmu.enterContinuousTimeMode()

            self._jac_fmu = jac_fmu
            self._jac_var_refs = {
                v.name: v.valueReference
                for v in model_desc.modelVariables
            }
            logger.info(
                "ModelicaJacobianMixin: ME FMU ready for directional derivatives."
            )
        except Exception as exc:
            logger.warning(
                f"ModelicaJacobianMixin: failed to initialize ME FMU "
                f"({exc}); falling back to global FD."
            )
            self._jac_fmu = None

    # ------------------------------------------------------------------
    # JacobianContributor protocol
    # ------------------------------------------------------------------

    def compute_jacobian_block(
        self,
        unit_type: str,
        config: dict,
        ref_state: "ReferenceState",
        row_offset: int,
        col_offset: int,
        n_vars: int,
    ) -> list[tuple[int, int, float]]:
        """Return Jacobian triplets via FMI 2.0 directional derivatives.

        Returns ``[]`` immediately if the ME FMU instance is unavailable,
        causing the caller to fall back to global finite differences.

        Args:
            unit_type:   Unit class name (not used — FMU handles all units).
            config:      Unit configuration dict.
            ref_state:   Inlet stream snapshot at the current Newton iterate.
            row_offset:  Global row offset for the outlet residual block.
            col_offset:  Global column offset for the inlet variable block.
            n_vars:      Variables per stream (3 + n_components).

        Returns:
            List of ``(global_row, global_col, value)`` triplets, or ``[]``.
        """
        if not getattr(self, "_jac_fmu", None):
            return []

        jac_fmu = self._jac_fmu
        var_refs = self._jac_var_refs

        outlet_stream = config.get("out", "")
        comps = sorted(ref_state.z.keys())

        # Build value reference lists for unknown (output) variables.
        # Variable naming convention from ModelicaProvider: {stream}_{T|P|flowrate}
        unknown_names = [
            f"{outlet_stream}_T",
            f"{outlet_stream}_P",
            f"{outlet_stream}_flowrate",
        ]
        unknown_refs = [var_refs[n] for n in unknown_names if n in var_refs]
        if not unknown_refs:
            return []

        # Build value reference lists for known (input) variables.
        # Inlet stream name is derived from config.
        inlet_stream = config.get("in", "")
        known_names = [
            f"{inlet_stream}_T",
            f"{inlet_stream}_P",
            f"{inlet_stream}_flowrate",
        ]
        known_refs = [var_refs[n] for n in known_names if n in var_refs]
        if not known_refs:
            return []

        # Set FMU state to the current reference state.
        try:
            state_refs_vals: list[tuple[int, float]] = []
            if f"{inlet_stream}_T" in var_refs:
                state_refs_vals.append((var_refs[f"{inlet_stream}_T"], ref_state.T))
            if f"{inlet_stream}_P" in var_refs:
                state_refs_vals.append((var_refs[f"{inlet_stream}_P"], ref_state.P))
            if f"{inlet_stream}_flowrate" in var_refs:
                state_refs_vals.append(
                    (var_refs[f"{inlet_stream}_flowrate"], ref_state.flowrate)
                )
            if state_refs_vals:
                refs, vals = zip(*state_refs_vals)
                jac_fmu.setReal(list(refs), list(vals))
        except Exception as exc:
            logger.debug(
                f"ModelicaJacobianMixin: setReal failed ({exc}); falling back to FD."
            )
            return []

        # For each known input variable, compute one directional derivative column.
        triplets: list[tuple[int, int, float]] = []
        for i, kref in enumerate(known_refs):
            try:
                dd = jac_fmu.getDirectionalDerivative(
                    unknown_refs, [kref], [1.0]
                )
            except Exception as exc:
                logger.debug(
                    f"ModelicaJacobianMixin: getDirectionalDerivative failed "
                    f"({exc}); aborting analytic Jacobian for this block."
                )
                return []

            # dd[j] = d(outlet_unknown[j]) / d(inlet_known[i])
            for j, val in enumerate(dd):
                if abs(val) > 1e-20:
                    triplets.append((row_offset + j, col_offset + i, float(val)))

        return triplets

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    def _teardown_jacobian_fmu(self) -> None:
        """Terminate and free the ME FMU instance.

        Called by ``ModelicaProvider.teardown()`` after adding this mixin.
        """
        jac_fmu = getattr(self, "_jac_fmu", None)
        if jac_fmu is not None:
            try:
                jac_fmu.terminate()
                jac_fmu.freeInstance()
            except Exception as exc:
                logger.debug(f"ModelicaJacobianMixin: ME FMU teardown error ({exc})")
            finally:
                self._jac_fmu = None
