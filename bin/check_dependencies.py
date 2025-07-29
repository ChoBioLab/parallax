#!/usr/bin/env python3
"""
Robust dependency checker for nf-core/parallax pipeline
Handles container environments and provides detailed diagnostics
"""

import importlib
import json
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pkg_resources
from packaging import version

# Configurable dependency requirements
REQUIRED_VERSIONS = {
    "scvi": "1.1.0",
    "spatialdata": "0.1.0",
    "scanpy": "1.9.0",
    "anndata": "0.9.0",
    "torch": "2.0.0",
    "numpy": "1.20.0",
    "pandas": "1.5.0",
    "scipy": "1.9.0",
    "matplotlib": "3.5.0",
    "seaborn": "0.11.0",
    "zarr": "2.12.0",
    "h5py": "3.7.0",
}

OPTIONAL_PACKAGES = {"cupy": "11.0.0", "rapids-singlecell": "0.9.0", "squidpy": "1.2.0"}


class DependencyChecker:
    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {
            "status": "unknown",
            "missing": [],
            "incompatible": [],
            "satisfied": [],
            "optional_available": [],
            "environment_info": self._get_environment_info(),
            "recommendations": [],
        }

    def _get_environment_info(self) -> Dict:
        """Get comprehensive environment information"""
        return {
            "python_version": sys.version,
            "python_executable": sys.executable,
            "platform": sys.platform,
            "is_container": self._detect_container(),
            "cuda_available": self._check_cuda(),
            "environment_variables": {
                k: v
                for k, v in os.environ.items()
                if any(
                    keyword in k.upper()
                    for keyword in [
                        "CUDA",
                        "PYTHON",
                        "CONDA",
                        "MAMBA",
                        "SINGULARITY",
                        "DOCKER",
                    ]
                )
            },
        }

    def _detect_container(self) -> bool:
        """Detect if running in container environment"""
        container_indicators = [
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
            for indicator in container_indicators
        )

    def _check_cuda(self) -> Dict:
        """Check CUDA availability"""
        cuda_info = {"available": False, "version": None, "devices": 0}
        try:
            import torch

            cuda_info["available"] = torch.cuda.is_available()
            if cuda_info["available"]:
                cuda_info["version"] = torch.version.cuda
                cuda_info["devices"] = torch.cuda.device_count()
        except ImportError:
            pass
        return cuda_info

    def check_package(
        self, package: str, min_version: str
    ) -> Tuple[bool, str, Optional[str]]:
        """Check individual package with robust error handling"""
        try:
            # Special handling for packages with non-standard version access
            if package == "scvi":
                import scvi

                current_version = scvi.__version__
            elif package == "spatialdata":
                import spatialdata

                current_version = spatialdata.__version__
            elif package == "torch":
                import torch

                current_version = torch.__version__
            else:
                # Standard approach
                try:
                    current_version = pkg_resources.get_distribution(package).version
                except pkg_resources.DistributionNotFound:
                    # Try importing directly
                    module = importlib.import_module(package)
                    current_version = getattr(module, "__version__", "unknown")

            # Version comparison
            if current_version == "unknown":
                return False, "version_unknown", current_version

            if version.parse(current_version) >= version.parse(min_version):
                return True, "satisfied", current_version
            else:
                return False, "incompatible", current_version

        except ImportError:
            return False, "missing", None
        except Exception as e:
            return False, f"error: {str(e)}", None

    def check_all_dependencies(self) -> bool:
        """Check all required and optional dependencies"""
        print("🔍 Checking dependencies for nf-core/parallax...")

        # Check required packages
        all_satisfied = True
        for package, min_version in REQUIRED_VERSIONS.items():
            satisfied, status, current_version = self.check_package(
                package, min_version
            )

            if satisfied:
                self.results["satisfied"].append(
                    {
                        "package": package,
                        "required": min_version,
                        "installed": current_version,
                    }
                )
                print(f"✅ {package}: {current_version} (>= {min_version})")
            elif status == "missing":
                self.results["missing"].append(package)
                print(f"❌ {package}: MISSING (required >= {min_version})")
                all_satisfied = False
            elif status == "incompatible":
                self.results["incompatible"].append(
                    {
                        "package": package,
                        "required": min_version,
                        "installed": current_version,
                    }
                )
                print(f"⚠️  {package}: {current_version} < {min_version} (INCOMPATIBLE)")
                all_satisfied = False
            else:
                print(f"❓ {package}: {status}")
                all_satisfied = False

        # Check optional packages
        print("\n🔍 Checking optional dependencies...")
        for package, min_version in OPTIONAL_PACKAGES.items():
            satisfied, status, current_version = self.check_package(
                package, min_version
            )
            if satisfied:
                self.results["optional_available"].append(
                    {"package": package, "version": current_version}
                )
                print(f"✅ {package}: {current_version} (optional)")
            else:
                print(f"➖ {package}: not available (optional)")

        self.results["status"] = "satisfied" if all_satisfied else "failed"
        return all_satisfied

    def generate_recommendations(self):
        """Generate installation recommendations"""
        if self.results["missing"] or self.results["incompatible"]:
            self.results["recommendations"].append(
                "Install missing/incompatible packages using conda/mamba:"
            )

            packages_to_install = []
            packages_to_install.extend(self.results["missing"])
            packages_to_install.extend(
                [pkg["package"] for pkg in self.results["incompatible"]]
            )

            conda_cmd = f"conda install -c conda-forge -c bioconda {' '.join(packages_to_install)}"
            self.results["recommendations"].append(conda_cmd)

        if self.results["environment_info"]["is_container"]:
            self.results["recommendations"].append(
                "Running in container - dependencies should be pre-installed"
            )

    def save_results(self):
        """Save detailed results to JSON file"""
        output_file = self.output_dir / "dependency_check.json"
        with open(output_file, "w") as f:
            json.dump(self.results, f, indent=2, default=str)

        # Also create a simple status file for Nextflow
        status_file = self.output_dir / "dependency_status.txt"
        with open(status_file, "w") as f:
            f.write(self.results["status"])

    def print_summary(self):
        """Print comprehensive summary"""
        print(f"\n📊 Dependency Check Summary:")
        print(f"   ✅ Satisfied: {len(self.results['satisfied'])}")
        print(f"   ❌ Missing: {len(self.results['missing'])}")
        print(f"   ⚠️  Incompatible: {len(self.results['incompatible'])}")
        print(f"   ➕ Optional available: {len(self.results['optional_available'])}")

        if self.results["recommendations"]:
            print(f"\n💡 Recommendations:")
            for rec in self.results["recommendations"]:
                print(f"   {rec}")


def main():
    """Main execution function"""
    import argparse

    parser = argparse.ArgumentParser(description="Check nf-core/parallax dependencies")
    parser.add_argument(
        "--output-dir", default=".", help="Output directory for results"
    )
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    args = parser.parse_args()

    if not args.quiet:
        print("🧬 nf-core/parallax Dependency Checker")
        print("=" * 50)

    checker = DependencyChecker(args.output_dir)

    try:
        success = checker.check_all_dependencies()
        checker.generate_recommendations()
        checker.save_results()

        if not args.quiet:
            checker.print_summary()

        if success:
            print("\n🎉 All required dependencies satisfied!")
            return 0
        else:
            print("\n💥 Dependency check failed!")
            return 1

    except Exception as e:
        print(f"💥 Dependency check failed with error: {e}")
        checker.results["status"] = "error"
        checker.results["error"] = str(e)
        checker.save_results()
        return 1


if __name__ == "__main__":
    sys.exit(main())
