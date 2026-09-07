"""``pf init`` — initialise the .processforge/ project directory."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess

import typer
from loguru import logger

from .common import (
    extract_providers,
    is_local_provider_url,
)
from ..lock import flowsheet_env_dir, read_lock, write_lock
from ..providers.registry import (
    get_provider_docker_image,
    get_provider_default_port,
    is_containerized,
    _PROVIDER_CATALOG,
)
from ..utils.validate_flowsheet import validate_flowsheet


def _migrate_legacy_env(pf_dir: str) -> None:
    """Move a legacy root-level lock.json / docker-compose.yml into a per-flowsheet dir.

    Older processforge stored a single environment at ``.processforge/lock.json``
    and ``.processforge/docker-compose.yml``. This relocates those into the hashed
    env dir derived from the flowsheet recorded in the legacy lock, so an existing
    repo keeps its provider environment after upgrading. The target always uses the
    same hashed env dir that ``pf init``/``read_lock`` locate, so the migrated env
    is never orphaned.
    """
    legacy_lock = os.path.join(pf_dir, "lock.json")
    legacy_compose = os.path.join(pf_dir, "docker-compose.yml")
    if not (os.path.exists(legacy_lock) or os.path.exists(legacy_compose)):
        return

    recorded = None
    if os.path.exists(legacy_lock):
        try:
            recorded = read_lock(pf_dir).get("flowsheet")
        except Exception:
            recorded = None

    # Always target the hashed env dir keyed by the recorded flowsheet path so it
    # matches what `pf init <flowsheet>` / `read_lock` expect on subsequent runs.
    env_dir = flowsheet_env_dir(pf_dir, recorded or "legacy")
    if recorded:
        logger.info(
            f"Migrating legacy .processforge/lock.json + docker-compose.yml into "
            f"'{os.path.relpath(env_dir, pf_dir)}/'."
        )
    else:
        logger.warning(
            "Migrating legacy .processforge/lock.json + docker-compose.yml into "
            f"'{os.path.relpath(env_dir, pf_dir)}/' (no recorded flowsheet found — "
            "re-run `pf init <flowsheet.json>` to restore the correct env dir)."
        )

    if os.path.exists(env_dir):
        logger.warning(
            f"Migration target {os.path.relpath(env_dir, pf_dir)}/ already exists — "
            "skipping legacy migration."
        )
        return

    os.makedirs(env_dir, exist_ok=True)
    for src in (legacy_lock, legacy_compose):
        if os.path.exists(src):
            shutil.move(src, os.path.join(env_dir, os.path.basename(src)))
    logger.info(
        "Migrated legacy .processforge/ environment into per-flowsheet dir "
        f"'{os.path.relpath(env_dir, pf_dir)}/'."
    )


def init(
    flowsheet: str | None = typer.Argument(
        default=None,
        help="Flowsheet JSON to initialise environment for (omit for scaffold only)",
    ),
    path: str = typer.Option(
        ".",
        "--path",
        help="Root directory to initialise in (default: current directory)",
    ),
    no_pull: bool = typer.Option(
        False,
        "--no-pull",
        help="Generate docker-compose.yml but skip pulling provider images",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Re-initialise even if the environment already exists",
    ),
) -> None:
    """Initialise the .processforge/ project directory."""
    from ..compose import generate_compose

    root = path or "."
    pf_dir = os.path.join(root, ".processforge")
    outputs_dir = os.path.join(root, "outputs")

    os.makedirs(pf_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    # Migrate a legacy single-environment layout (.processforge/lock.json and
    # .processforge/docker-compose.yml at the root) into a per-flowsheet env
    # dir so existing repos aren't silently broken by the new structure.
    _migrate_legacy_env(pf_dir)

    # Write config.json (always)
    config_path = os.path.join(pf_dir, "config.json")
    if not os.path.exists(config_path):
        default_config = {
            "version": 1,
            "default_backend": "scipy",
            "outputs_dir": "outputs",
        }
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(default_config, f, indent=2)
        logger.info(f"Created {config_path}")
    else:
        logger.info(f"{config_path} already exists — skipped.")

    # Honour a configured outputs_dir for the rest of init.
    outputs_dir_name = "outputs"
    try:
        with open(config_path, encoding="utf-8") as f:
            outputs_dir_name = json.load(f).get("outputs_dir", "outputs")
    except Exception:
        pass
    outputs_dir = os.path.join(root, outputs_dir_name)

    # Remove stale .pfstate snapshot directories from outputs/
    stale_count = 0
    for entry in os.listdir(outputs_dir):
        if entry.endswith(".pfstate"):
            stale = os.path.join(outputs_dir, entry)
            if os.path.isdir(stale):
                shutil.rmtree(stale)
                stale_count += 1
    if stale_count:
        logger.info(f"Removed {stale_count} stale snapshot(s) from {outputs_dir}/.")

    # No flowsheet → scaffold only
    if not flowsheet:
        logger.info(".processforge/ initialised successfully.")
        logger.info("To set up providers: pf init <flowsheet.json>")
        return

    # Read providers from flowsheet
    flowsheet_path = flowsheet
    if not os.path.exists(flowsheet_path):
        logger.error(f"Flowsheet '{flowsheet_path}' not found.")
        raise typer.Exit(code=1)

    try:
        validate_flowsheet(flowsheet_path)
    except Exception as e:
        logger.error(f"Failed to validate flowsheet '{flowsheet_path}': {e}")
        raise typer.Exit(code=1)

    providers = extract_providers(flowsheet_path)
    logger.info(f"Reading providers from {flowsheet_path}...")

    # Categorize providers, reusing the shared local/remote/pip classification.
    local_docker_providers: dict[str, dict] = {}
    remote_docker_providers: dict[str, dict] = {}
    pip_providers: dict[str, dict] = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = cfg.get("url")
            if is_local_provider_url(url):
                port = get_provider_default_port(ptype) or 9000
                if not url:
                    url = f"http://localhost:{port}"
                local_docker_providers[name] = {
                    "type": ptype,
                    "url": url,
                    "docker_image": cfg.get("docker_image")
                    or get_provider_docker_image(ptype),
                    "port": port,
                }
                logger.info(f"  {name}: type={ptype}, url={url} (Docker, local)")
            else:
                remote_docker_providers[name] = {
                    "type": ptype,
                    "url": url,
                    "docker_image": cfg.get("docker_image")
                    or get_provider_docker_image(ptype),
                }
                logger.info(
                    f"  {name}: type={ptype}, url={url} (Docker, remote — skipping compose)"
                )
        else:
            pip_providers[name] = {"type": ptype}
            logger.info(f"  {name}: type={ptype} (pip)")

    # Validate pip providers are importable
    for name, info in pip_providers.items():
        ptype = info["type"]
        catalog = _PROVIDER_CATALOG.get(ptype)
        module = catalog.module if catalog else ""
        try:
            importlib.util.find_spec(module)
            logger.info(f"  {name} — importable")
        except (ModuleNotFoundError, ValueError):
            dep = catalog.optional_dep if catalog else None
            hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
            logger.warning(f"  {name} — not installed. Install with: {hint}")

    # Generate compose + pull images only for locally-managed Docker providers.
    # Providers with an explicit remote URL are assumed to be running elsewhere
    # (e.g. a cloud deployment of the ghcr.io image) and are not touched here.
    if local_docker_providers:
        env_dir = flowsheet_env_dir(pf_dir, flowsheet_path)
        compose_path = os.path.join(env_dir, "docker-compose.yml")
        if os.path.exists(compose_path) and not force:
            logger.warning(
                f"Environment already initialized — reinitializing from {flowsheet_path}"
            )

        generate_compose(
            pf_dir, local_docker_providers, outputs_dir, flowsheet=flowsheet_path
        )
        logger.info(f"Generated {compose_path}")

        if no_pull:
            logger.info("Skipping Docker image pull (--no-pull).")
        else:
            # Attempt docker compose pull. Capture output and stream it via the
            # logger; on timeout the child is reaped (subprocess.run kills it)
            # so no orphaned/blocked process is left behind.
            try:
                logger.info("Pulling Docker images (this may take a while)...")
                result = subprocess.run(
                    ["docker", "compose", "-f", compose_path, "pull"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=600,
                )
                for line in result.stdout.splitlines():
                    if line.strip():
                        logger.info(line)
                if result.returncode == 0:
                    logger.info("Pulled Docker images.")
                else:
                    logger.warning(
                        f"docker compose pull failed with exit code {result.returncode}."
                    )
            except FileNotFoundError:
                logger.warning(
                    "Docker not found. Install Docker to use containerized providers."
                )
            except subprocess.TimeoutExpired:
                logger.warning("docker compose pull timed out after 600s.")

        logger.info(
            "To start the containerized provider(s), run:\n"
            f"  docker compose -f {compose_path} up -d"
        )
        if any(info["type"] == "openmc" for info in local_docker_providers.values()):
            logger.info(
                "Running OpenMC — set OPENMC_DATA_ROOT to mount your OpenMC "
                "data/workspace into the container (defaults to 'outputs'):\n"
                f"  OPENMC_DATA_ROOT=/path/to/data docker compose -f {compose_path} up -d"
            )
    elif remote_docker_providers:
        logger.info(
            "All containerized providers use remote URLs — skipping Docker setup."
        )
    else:
        logger.info("No containerized providers — skipping Docker setup.")

    # Write lock file, reusing the shared URL resolution for the recorded url.
    lock_providers: dict[str, dict] = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = cfg.get("url")
            if not url:
                port = get_provider_default_port(ptype) or 9000
                url = f"http://localhost:{port}"
            lock_providers[name] = {
                "docker_image": cfg.get("docker_image")
                or get_provider_docker_image(ptype),
                "url": url,
            }
        else:
            lock_providers[name] = {
                "docker_image": None,
                "url": None,
            }

    from .. import __version__ as pf_version

    write_lock(pf_dir, flowsheet_path, lock_providers, pf_version)
    logger.info(
        f"Wrote {os.path.join(flowsheet_env_dir(pf_dir, flowsheet_path), 'lock.json')}"
    )
    logger.info(".processforge/ initialised successfully.")
