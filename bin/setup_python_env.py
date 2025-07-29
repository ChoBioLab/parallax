#!/usr/bin/env python3
"""
Comprehensive Python environment setup for nf-core/parallax pipeline
Handles container environments, GPU configuration, and performance optimization
"""

import os
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Dict, Optional


class EnvironmentSetup:
    def __init__(self):
        self.config = {
            "is_container": self._detect_container(),
            "container_type": self._detect_container_type(),
            "has_gpu": self._detect_gpu(),
            "temp_dir": None,
            "applied_settings": [],
        }

    def _detect_container(self) -> bool:
        """Detect if running in any container environment"""
        indicators = [
            "/.singularity.d",
            "/.dockerenv",
            os.environ.get("SINGULARITY_CONTAINER"),
            os.environ.get("APPTAINER_CONTAINER"),
            os.environ.get("CONTAINER_ENGINE"),
        ]
        return any(
            os.path.exists(indicator)
            if isinstance(indicator, str) and indicator.startswith("/")
            else bool(indicator)
            for indicator in indicators
        )

    def _detect_container_type(self) -> Optional[str]:
        """Detect specific container type"""
        if os.path.exists("/.singularity.d") or os.environ.get("SINGULARITY_CONTAINER"):
            return "singularity"
        elif os.path.exists("/.dockerenv"):
            return "docker"
        elif os.environ.get("APPTAINER_CONTAINER"):
            return "apptainer"
        return None

    def _detect_gpu(self) -> bool:
        """Detect GPU availability"""
        try:
            import torch

            return torch.cuda.is_available()
        except ImportError:
            return False

    def setup_numba_config(self):
        """Configure Numba for container and performance"""
        if self.config["is_container"]:
            # Container-safe configuration
            temp_dir = tempfile.mkdtemp(prefix="numba_cache_")
            os.environ["NUMBA_CACHE_DIR"] = temp_dir
            self.config["temp_dir"] = temp_dir

            # Disable problematic features in containers
            os.environ["NUMBA_DISABLE_PARALLEL"] = "1"
            os.environ["NUMBA_DISABLE_TBB"] = "1"
            os.environ["NUMBA_DISABLE_OPENMP"] = "1"

            # Keep JIT compilation enabled for performance
            os.environ["NUMBA_DISABLE_JIT"] = "0"

            self.config["applied_settings"].append(f"Numba cache: {temp_dir}")
            self.config["applied_settings"].append(
                "Numba parallel features disabled (container-safe)"
            )
        else:
            # Native environment - optimize for performance
            os.environ["NUMBA_DISABLE_JIT"] = "0"

            # Enable threading if available
            if "OMP_NUM_THREADS" not in os.environ:
                import multiprocessing

                os.environ["OMP_NUM_THREADS"] = str(multiprocessing.cpu_count())

            self.config["applied_settings"].append(
                "Numba full performance mode enabled"
            )

    def setup_matplotlib_config(self):
        """Configure matplotlib for headless operation"""
        # Always use Agg backend for reproducibility
        os.environ["MPLBACKEND"] = "Agg"

        # Configure matplotlib cache
        if self.config["is_container"]:
            mpl_cache = os.path.join(
                tempfile.gettempdir(), f"matplotlib_cache_{os.getpid()}"
            )
            os.makedirs(mpl_cache, exist_ok=True)
            os.environ["MPLCONFIGDIR"] = mpl_cache
            self.config["applied_settings"].append(f"Matplotlib cache: {mpl_cache}")

        self.config["applied_settings"].append("Matplotlib Agg backend enabled")

    def setup_pytorch_config(self):
        """Configure PyTorch for optimal performance"""
        if self.config["has_gpu"]:
            # GPU-specific settings
            os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:512"

            # Enable optimized attention if available
            os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"

            self.config["applied_settings"].append("PyTorch GPU optimizations enabled")
        else:
            # CPU-specific optimizations
            import multiprocessing

            cpu_count = multiprocessing.cpu_count()

            if "OMP_NUM_THREADS" not in os.environ:
                os.environ["OMP_NUM_THREADS"] = str(cpu_count)

            if "MKL_NUM_THREADS" not in os.environ:
                os.environ["MKL_NUM_THREADS"] = str(cpu_count)

            self.config["applied_settings"].append(
                f"PyTorch CPU optimizations enabled ({cpu_count} threads)"
            )

    def setup_memory_config(self):
        """Configure memory management"""
        # Disable memory mapping warnings in containers
        if self.config["is_container"]:
            os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
            self.config["applied_settings"].append(
                "HDF5 file locking disabled (container-safe)"
            )

        # Configure garbage collection for large datasets
        import gc

        gc.set_threshold(700, 10, 10)  # More aggressive GC for memory efficiency

    def setup_warning_filters(self):
        """Configure warning filters"""
        # Suppress common warnings that clutter output
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
        warnings.filterwarnings("ignore", category=UserWarning, module="anndata")

        # Keep important warnings
        warnings.filterwarnings("default", category=RuntimeWarning)
        warnings.filterwarnings("default", category=ImportWarning)

        self.config["applied_settings"].append("Warning filters configured")

    def setup_reproducibility(self):
        """Configure settings for reproducibility"""
        # Python hash randomization
        if "PYTHONHASHSEED" not in os.environ:
            os.environ["PYTHONHASHSEED"] = "0"

        # Disable user site packages
        os.environ["PYTHONNOUSERSITE"] = "1"

        # Configure numpy random seed (will be overridden by pipeline)
        os.environ["NUMPY_RANDOM_SEED"] = "42"

        self.config["applied_settings"].append("Reproducibility settings applied")

    def setup_container_specific(self):
        """Container-specific configurations"""
        if not self.config["is_container"]:
            return

        container_type = self.config["container_type"]

        if container_type == "singularity":
            # Singularity-specific settings
            os.environ["SINGULARITY_BIND"] = (
                os.environ.get("SINGULARITY_BIND", "") + ":/tmp"
            )
            self.config["applied_settings"].append(
                "Singularity /tmp binding configured"
            )

        elif container_type == "docker":
            # Docker-specific settings
            os.environ["DOCKER_TMPDIR"] = "/tmp"
            self.config["applied_settings"].append("Docker temp directory configured")

    def apply_all_configurations(self):
        """Apply all environment configurations"""
        print(
            "🔧 Setting up Python environment for nf-core/parallax...", file=sys.stderr
        )

        self.setup_numba_config()
        self.setup_matplotlib_config()
        self.setup_pytorch_config()
        self.setup_memory_config()
        self.setup_warning_filters()
        self.setup_reproducibility()
        self.setup_container_specific()

        return self.config

    def print_summary(self):
        """Print configuration summary"""
        print(f"\n🔧 Environment Configuration Summary:", file=sys.stderr)
        print(
            f"   Container: {self.config['container_type'] or 'Native'}",
            file=sys.stderr,
        )
        print(f"   GPU Available: {self.config['has_gpu']}", file=sys.stderr)
        print(f"   Applied Settings:", file=sys.stderr)
        for setting in self.config["applied_settings"]:
            print(f"     • {setting}", file=sys.stderr)

    def save_config(self, output_path: str = "environment_config.json"):
        """Save configuration to file"""
        import json

        with open(output_path, "w") as f:
            json.dump(self.config, f, indent=2, default=str)


def main():
    """Main execution function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Setup Python environment for nf-core/parallax"
    )
    parser.add_argument("--output-config", help="Save configuration to JSON file")
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    args = parser.parse_args()

    setup = EnvironmentSetup()

    try:
        config = setup.apply_all_configurations()

        if not args.quiet:
            setup.print_summary()

        if args.output_config:
            setup.save_config(args.output_config)

        print("✅ Environment setup complete!", file=sys.stderr)
        return 0

    except Exception as e:
        print(f"💥 Environment setup failed: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
