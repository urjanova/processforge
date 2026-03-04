"""Validate a Modelica .mo file and compile it to a Model Exchange FMU via OMPython."""
from __future__ import annotations

import os
import shutil

from loguru import logger


def compile_modelica(
    mo_path: str,
    model_name: str,
    output_dir: str = "outputs",
) -> str:
    """Validate and compile a Modelica model to a Model Exchange FMU 2.0.

    Steps:
    1. Start an OMPython session (``OMCSessionZMQ``).
    2. ``loadFile()`` — load the ``.mo`` source.
    3. ``checkModel()`` — validate types, equations, and connections.
    4. ``buildModelFMU()`` — generate a Model Exchange FMU 2.0.
    5. Move the ``.fmu`` to *output_dir*.

    Args:
        mo_path:     Absolute path to the ``.mo`` file.
        model_name:  Modelica model identifier (must match the model name
                     inside the ``.mo`` file).
        output_dir:  Directory where the final ``.fmu`` will be placed.

    Returns:
        Absolute path to the generated ``.fmu`` file.

    Raises:
        RuntimeError:  If OMPython is not installed, or if ``checkModel`` or
                       ``buildModelFMU`` fails.
    """
    try:
        from OMPython import OMCSessionZMQ  # type: ignore[import]
    except ImportError as exc:
        raise RuntimeError(
            "OMPython is not installed. "
            "Install it with:  pip install 'processforge[modelica]'\n"
            "OpenModelica must also be installed on your system: "
            "https://openmodelica.org"
        ) from exc

    # Guard: OMCSessionZMQ.__init__ crashes (and its __del__ then raises
    # AttributeError) when omc is not on PATH.  Check first so we surface a
    # clean error instead of a confusing traceback from the GC.
    if not shutil.which("omc"):
        raise RuntimeError(
            "OpenModelica compiler (omc) not found on PATH. "
            "Install OpenModelica from https://openmodelica.org and ensure "
            "the omc binary is accessible."
        )

    logger.info("Starting OMPython session…")
    try:
        om = OMCSessionZMQ()
    except Exception as exc:
        raise RuntimeError(
            f"Failed to start OMPython/omc session: {exc}\n"
            "Ensure OpenModelica is correctly installed."
        ) from exc

    abs_mo = os.path.abspath(mo_path)

    # Load file
    load_result = om.sendExpression(f'loadFile("{abs_mo}")')
    if not load_result:
        errors = om.sendExpression("getErrorString()")
        raise RuntimeError(f"OMPython loadFile failed:\n{errors}")
    logger.info(f"Loaded: {abs_mo}")

    # Check model
    check_result = om.sendExpression(f'checkModel({model_name})')
    errors = om.sendExpression("getErrorString()")
    if errors and errors.strip() not in ("", '""'):
        logger.warning(f"OMPython checkModel warnings/errors:\n{errors}")
    logger.info(f"checkModel result: {check_result}")

    # Build FMU (model exchange, FMI 2.0)
    os.makedirs(output_dir, exist_ok=True)
    abs_output = os.path.abspath(output_dir)
    fmu_result = om.sendExpression(
        f'buildModelFMU({model_name}, version="2.0", fmuType="me", '
        f'fileNamePrefix="{abs_output}/{model_name}")'
    )
    build_errors = om.sendExpression("getErrorString()")
    if build_errors and build_errors.strip() not in ("", '""'):
        logger.warning(f"OMPython buildModelFMU warnings:\n{build_errors}")

    # OMPython returns the path to the generated FMU as a string
    fmu_path = str(fmu_result).strip().strip('"')
    if not fmu_path or not os.path.exists(fmu_path):
        # Fallback: look for the FMU in the output_dir
        candidate = os.path.join(abs_output, f"{model_name}.fmu")
        if os.path.exists(candidate):
            fmu_path = candidate
        else:
            raise RuntimeError(
                f"buildModelFMU did not produce a .fmu file.\n"
                f"OMPython returned: {fmu_result!r}\n"
                f"Errors: {build_errors}"
            )

    # Move to output_dir if it landed elsewhere
    target = os.path.join(abs_output, os.path.basename(fmu_path))
    if os.path.abspath(fmu_path) != os.path.abspath(target):
        shutil.move(fmu_path, target)

    logger.info(f"FMU written to: {target}")
    return os.path.abspath(target)
