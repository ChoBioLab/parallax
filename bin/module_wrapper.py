#!/usr/bin/env python3
"""
Comprehensive module wrapper for nf-core/parallax pipeline
Provides standardized logging, error handling, and resource monitoring
"""

import json
import logging
import os
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import psutil
import yaml


class ModuleWrapper:
    def __init__(self, module_name: str, sample_id: str, output_dir: str = "."):
        self.module_name = module_name
        self.sample_id = sample_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.start_time = datetime.now()
        self.process = psutil.Process()
        self.initial_memory = self.process.memory_info().rss

        self.metrics = {
            "module": module_name,
            "sample_id": sample_id,
            "start_time": self.start_time.isoformat(),
            "end_time": None,
            "duration_seconds": None,
            "peak_memory_mb": 0,
            "cpu_percent": [],
            "status": "running",
            "error": None,
        }

        self.setup_logging()

    def setup_logging(self):
        """Setup comprehensive logging with multiple handlers"""
        # Create formatters
        detailed_formatter = logging.Formatter(
            "%(asctime)s | %(name)s | %(levelname)s | %(funcName)s:%(lineno)d | %(message)s"
        )
        simple_formatter = logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s"
        )

        # Setup logger
        self.logger = logging.getLogger(f"nf-core.parallax.{self.module_name}")
        self.logger.setLevel(logging.INFO)

        # Clear any existing handlers
        self.logger.handlers.clear()

        # File handler for detailed logs
        log_file = self.output_dir / f"{self.module_name}_{self.sample_id}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(detailed_formatter)
        self.logger.addHandler(file_handler)

        # Console handler for important messages
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(simple_formatter)
        self.logger.addHandler(console_handler)

        # Error file handler
        error_file = self.output_dir / f"{self.module_name}_{self.sample_id}_errors.log"
        error_handler = logging.FileHandler(error_file)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(detailed_formatter)
        self.logger.addHandler(error_handler)

    def log_start(self, additional_info: Optional[Dict] = None):
        """Log module start with system information"""
        self.logger.info(f"🚀 Starting {self.module_name} for sample {self.sample_id}")
        self.logger.info(f"📅 Timestamp: {self.start_time.isoformat()}")
        self.logger.info(
            f"🖥️  System: {psutil.cpu_count()} CPUs, {psutil.virtual_memory().total / (1024**3):.1f} GB RAM"
        )

        # Log GPU information if available
        try:
            import torch

            if torch.cuda.is_available():
                gpu_count = torch.cuda.device_count()
                gpu_memory = torch.cuda.get_device_properties(0).total_memory / (
                    1024**3
                )
                self.logger.info(
                    f"🎮 GPU: {gpu_count} device(s), {gpu_memory:.1f} GB memory"
                )
        except ImportError:
            pass

        # Log additional information
        if additional_info:
            for key, value in additional_info.items():
                self.logger.info(f"ℹ️  {key}: {value}")

        # Log environment variables
        relevant_env_vars = {
            k: v
            for k, v in os.environ.items()
            if any(
                keyword in k.upper()
                for keyword in ["CUDA", "OMP", "MKL", "NUMBA", "PYTORCH"]
            )
        }
        if relevant_env_vars:
            self.logger.debug(f"Environment variables: {relevant_env_vars}")

    def monitor_resources(self):
        """Monitor resource usage during execution"""
        try:
            # Memory usage
            memory_mb = self.process.memory_info().rss / (1024 * 1024)
            self.metrics["peak_memory_mb"] = max(
                self.metrics["peak_memory_mb"], memory_mb
            )

            # CPU usage
            cpu_percent = self.process.cpu_percent()
            self.metrics["cpu_percent"].append(cpu_percent)

            # Log if memory usage is high
            if memory_mb > 1000:  # > 1GB
                self.logger.debug(f"📊 Memory usage: {memory_mb:.1f} MB")

        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

    def log_completion(
        self, success: bool = True, additional_metrics: Optional[Dict] = None
    ):
        """Log completion with comprehensive metrics"""
        end_time = datetime.now()
        duration = end_time - self.start_time

        self.metrics.update(
            {
                "end_time": end_time.isoformat(),
                "duration_seconds": duration.total_seconds(),
                "status": "completed" if success else "failed",
            }
        )

        if additional_metrics:
            self.metrics.update(additional_metrics)

        # Calculate average CPU usage
        if self.metrics["cpu_percent"]:
            avg_cpu = sum(self.metrics["cpu_percent"]) / len(
                self.metrics["cpu_percent"]
            )
            self.metrics["avg_cpu_percent"] = avg_cpu

        if success:
            self.logger.info(
                f"✅ Completed {self.module_name} for sample {self.sample_id}"
            )
        else:
            self.logger.error(
                f"❌ Failed {self.module_name} for sample {self.sample_id}"
            )

        self.logger.info(f"⏱️  Duration: {duration.total_seconds():.2f} seconds")
        self.logger.info(f"📊 Peak memory: {self.metrics['peak_memory_mb']:.1f} MB")

        if "avg_cpu_percent" in self.metrics:
            self.logger.info(f"🖥️  Average CPU: {self.metrics['avg_cpu_percent']:.1f}%")

    def handle_error(self, error: Exception, create_placeholder: bool = True):
        """Handle errors with comprehensive logging and recovery"""
        error_info = {
            "error_type": type(error).__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
            "timestamp": datetime.now().isoformat(),
        }

        self.metrics["error"] = error_info
        self.metrics["status"] = "failed"

        # Log error
        self.logger.error(f"💥 Error in {self.module_name}: {error}")
        self.logger.debug(f"Full traceback:\n{traceback.format_exc()}")

        # Save error report
        error_file = self.output_dir / f"error_report_{self.sample_id}.json"
        with open(error_file, "w") as f:
            json.dump(error_info, f, indent=2)

        # Create placeholder outputs if requested
        if create_placeholder:
            self.create_placeholder_outputs()

        return error_info

    def create_placeholder_outputs(self):
        """Create standardized placeholder outputs for failed modules"""
        placeholder_info = {
            "status": "failed",
            "module": self.module_name,
            "sample_id": self.sample_id,
            "timestamp": datetime.now().isoformat(),
            "reason": "Module execution failed",
        }

        # Module-specific placeholders
        if self.module_name in ["RESOLVI_TRAIN", "SCVIVA_TRAIN"]:
            # Create empty model placeholder
            model_file = self.output_dir / "model_placeholder.pkl"
            model_file.touch()

            info_file = self.output_dir / "model_info.json"
            with open(info_file, "w") as f:
                json.dump(placeholder_info, f, indent=2)

        elif self.module_name in ["RESOLVI_ANALYZE", "SCVIVA_ANALYZE"]:
            # Create empty results placeholder
            results_file = self.output_dir / "results_placeholder.h5ad"
            results_file.touch()

            info_file = self.output_dir / "analysis_info.json"
            with open(info_file, "w") as f:
                json.dump(placeholder_info, f, indent=2)

        elif self.module_name in ["RESOLVI_VISUALIZE"]:
            # Create empty plot placeholders
            for plot_type in ["spatial", "umap", "metrics"]:
                plot_file = self.output_dir / f"{plot_type}_placeholder.png"
                plot_file.touch()

        # Always create versions file
        self.create_versions_file()

    def create_versions_file(self):
        """Create comprehensive versions.yml with error handling"""
        try:
            versions = {self.module_name: self._get_package_versions()}
        except Exception as e:
            versions = {
                self.module_name: {
                    "error": f"Could not determine versions: {e}",
                    "python": sys.version.split()[0],
                    "timestamp": datetime.now().isoformat(),
                }
            }

        versions_file = self.output_dir / "versions.yml"
        with open(versions_file, "w") as f:
            yaml.dump(versions, f, default_flow_style=False, sort_keys=True)

    def _get_package_versions(self) -> Dict[str, str]:
        """Get versions of all relevant packages with robust error handling"""
        import pkg_resources

        packages_to_check = [
            "scvi-tools",
            "spatialdata",
            "scanpy",
            "anndata",
            "torch",
            "numpy",
            "pandas",
            "scipy",
            "matplotlib",
            "seaborn",
            "zarr",
            "h5py",
            "numba",
            "packaging",
        ]

        versions = {
            "python": sys.version.split()[0],
            "platform": sys.platform,
            "timestamp": datetime.now().isoformat(),
        }

        for package in packages_to_check:
            try:
                # Try standard approach first
                version = pkg_resources.get_distribution(package).version
                versions[package] = version
            except pkg_resources.DistributionNotFound:
                try:
                    # Try importing directly for packages with non-standard names
                    if package == "scvi-tools":
                        import scvi

                        versions[package] = scvi.__version__
                    elif package == "torch":
                        import torch

                        versions[package] = torch.__version__
                    else:
                        # Try generic import
                        module = __import__(package.replace("-", "_"))
                        versions[package] = getattr(module, "__version__", "unknown")
                except (ImportError, AttributeError):
                    versions[package] = "not_installed"
            except Exception:
                versions[package] = "unknown"

        return versions

    def save_metrics(self):
        """Save execution metrics to JSON file"""
        metrics_file = self.output_dir / f"execution_metrics_{self.sample_id}.json"
        with open(metrics_file, "w") as f:
            json.dump(self.metrics, f, indent=2, default=str)

    def execute_with_monitoring(self, func: Callable, *args, **kwargs) -> Any:
        """Execute function with comprehensive monitoring"""
        self.log_start()

        try:
            # Start resource monitoring
            monitor_interval = 30  # seconds
            last_monitor = time.time()

            # Execute function
            result = func(*args, **kwargs)

            # Monitor resources periodically during execution
            current_time = time.time()
            if current_time - last_monitor > monitor_interval:
                self.monitor_resources()
                last_monitor = current_time

            self.log_completion(success=True)
            self.create_versions_file()
            self.save_metrics()

            return result

        except Exception as e:
            self.handle_error(e, create_placeholder=True)
            self.log_completion(success=False)
            self.save_metrics()
            raise


def safe_execution(
    func: Callable,
    module_name: str,
    sample_id: str,
    output_dir: str = ".",
    *args,
    **kwargs,
) -> Any:
    """Wrapper for safe execution with comprehensive monitoring"""
    wrapper = ModuleWrapper(module_name, sample_id, output_dir)
    return wrapper.execute_with_monitoring(func, *args, **kwargs)


def main():
    """Main execution for testing"""
    import argparse

    parser = argparse.ArgumentParser(description="Module wrapper for nf-core/parallax")
    parser.add_argument("--module", required=True, help="Module name")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    args = parser.parse_args()

    wrapper = ModuleWrapper(args.module, args.sample, args.output_dir)
    wrapper.log_start()
    wrapper.log_completion()
    wrapper.create_versions_file()
    wrapper.save_metrics()


if __name__ == "__main__":
    main()
