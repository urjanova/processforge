"""``pf init`` — initialise the .processforge/ project directory."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import shutil
import subprocess

from loguru import logger

from .common import extract_providers


def add_init_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``init`` subcommand."""
    parser.add_argument(
        "flowsheet", nargs="?", default=None,
        help="Flowsheet JSON to initialise environment for (omit for scaffold only)",
    )
    parser.add_argument(
        "--path", default=".",
        help="Root directory to initialise in (default: current directory)",
    )


def cmd_init(args: argparse.Namespace) -> None:
    """Initialise the .processforge/ project directory."""
    from ..providers.registry import (
        is_containerized,
        get_provider_docker_image,
        get_provider_default_port,
        _PROVIDER_CATALOG,
    )
    from ..lock import write_lock
    from ..compose import generate_compose

    root = args.path or "."
    pf_dir = os.path.join(root, ".processforge")
    outputs_dir = os.path.join(root, "outputs")

    os.makedirs(pf_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

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

    # No flowsheet → scaffold only
    if not args.flowsheet:
        logger.info(".processforge/ initialised successfully.")
        logger.info("To set up providers: pf init <flowsheet.json>")
        return

    # Read providers from flowsheet
    flowsheet_path = args.flowsheet
    if not os.path.exists(flowsheet_path):
        logger.error(f"Flowsheet '{flowsheet_path}' not found.")
        raise SystemExit(1)

    providers = extract_providers(flowsheet_path)
    logger.info(f"Reading providers from {flowsheet_path}...")

    # Categorize providers
    pip_providers: dict[str, dict] = {}
    docker_providers: dict[str, dict] = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = cfg.get("url")
            if not url:
                port = get_provider_default_port(ptype) or 9000
                url = f"http://localhost:{port}"
            docker_providers[name] = {
                "type": ptype,
                "url": url,
                "docker_image": cfg.get("docker_image") or get_provider_docker_image(ptype),
                "port": get_provider_default_port(ptype),
            }
            logger.info(f"  {name}: type={ptype}, url={url} (Docker)")
        else:
            pip_providers[name] = {"type": ptype}
            logger.info(f"  {name}: type={ptype} (pip)")

    # Validate pip providers are importable
    for name, info in pip_providers.items():
        ptype = info["type"]
        catalog = _PROVIDER_CATALOG.get(ptype, {})
        module = catalog.get("module", "")
        try:
            importlib.util.find_spec(module)
            logger.info(f"  {name} — importable")
        except (ModuleNotFoundError, ValueError):
            dep = catalog.get("optional_dep")
            hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
            logger.warning(f"  {name} — not installed. Install with: {hint}")

    # Generate compose for Docker providers
    if docker_providers:
        compose_path = os.path.join(pf_dir, "docker-compose.yml")
        if os.path.exists(compose_path):
            logger.warning(
                f"Environment already initialized — reinitializing from {flowsheet_path}"
            )

        generate_compose(pf_dir, docker_providers, outputs_dir)
        logger.info(f"Generated {compose_path}")

        # Attempt docker compose pull
        try:
            result = subprocess.run(
                ["docker", "compose", "-f", compose_path, "pull"],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode == 0:
                logger.info("Pulled Docker images.")
            else:
                logger.warning(f"docker compose pull failed: {result.stderr}")
        except FileNotFoundError:
            logger.warning(
                "Docker not found. Install Docker to use containerized providers."
            )
        except subprocess.TimeoutExpired:
            logger.warning("docker compose pull timed out after 120s.")
    else:
        logger.info("No containerized providers — skipping Docker setup.")

    # Write lock file
    lock_providers: dict[str, dict] = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            lock_providers[name] = {
                "docker_image": cfg.get("docker_image") or get_provider_docker_image(ptype),
                "url": cfg.get("url")
                or f"http://localhost:{get_provider_default_port(ptype) or 9000}",
            }
        else:
            lock_providers[name] = {
                "docker_image": None,
                "url": None,
            }

    from .. import __version__ as pf_version

    write_lock(pf_dir, flowsheet_path, lock_providers, pf_version)
    logger.info(f"Wrote {os.path.join(pf_dir, 'lock.json')}")
    logger.info(".processforge/ initialised successfully.")
