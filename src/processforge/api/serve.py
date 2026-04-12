"""Standalone entry point for the ProcessForge API server.

Usage::

    # Via registered script
    pf-serve [--host 0.0.0.0] [--port 9000]

    # Or directly
    python -m processforge.api.serve [--host 0.0.0.0] [--port 9000]

    # Or via Docker (override entrypoint)
    docker run -p 9000:9000 ghcr.io/urjanova/processforge:latest pf-serve

AWS credentials are read from the standard environment variables
(AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_DEFAULT_REGION) or from an
attached IAM role — no code changes required.
"""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Start the ProcessForge FastAPI server",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--host", default="0.0.0.0", help="Bind host")
    parser.add_argument("--port", type=int, default=9000, help="Bind port")
    args = parser.parse_args()

    import uvicorn

    from processforge.api.app import app

    uvicorn.run(app, host=args.host, port=args.port, workers=1)


if __name__ == "__main__":
    main()
