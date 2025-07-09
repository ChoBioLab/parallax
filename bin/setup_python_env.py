#!/usr/bin/env python3
"""
Shared environment setup for nf-core/resolvinf pipeline
Handles container-safe numba configuration
"""

import os
import sys
import tempfile
import warnings
from pathlib import Path


def setup_container_environment():
    """Configure Python environment for container execution"""

    # Detect container environment
    is_container = any(
        [
            os.path.exists("/.singularity.d"),
            os.path.exists("/.dockerenv"),
            os.environ.get("SINGULARITY_CONTAINER"),
            os.environ.get("APPTAINER_CONTAINER"),
            os.environ.get("CONTAINER_ENGINE"),
        ]
    )

    if is_container:
        # Container-safe numba configuration
        numba_cache = os.path.join("/tmp", f"numba_cache_{os.getpid()}")
        os.makedirs(numba_cache, exist_ok=True)
        os.environ["NUMBA_CACHE_DIR"] = numba_cache

        # Disable problematic features, keep JIT
        os.environ["NUMBA_DISABLE_PARALLEL"] = "1"
        os.environ["NUMBA_DISABLE_TBB"] = "1"

        # Matplotlib backend for headless containers
        os.environ["MPLBACKEND"] = "Agg"

        print(
            f"Container environment detected. Numba cache: {numba_cache}",
            file=sys.stderr,
        )
    else:
        # Native environment - full performance
        os.environ["NUMBA_DISABLE_JIT"] = "0"
        print(
            "Native environment detected. Full numba performance enabled.",
            file=sys.stderr,
        )

    # Common settings
    os.environ["PYTHONNOUSERSITE"] = "1"
    warnings.filterwarnings("ignore")

    return is_container


if __name__ == "__main__":
    setup_container_environment()
