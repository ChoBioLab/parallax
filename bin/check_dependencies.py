# bin/check_dependencies.py
#!/usr/bin/env python3
import importlib
import sys

import pkg_resources
from packaging import version

REQUIRED_VERSIONS = {
    "scvi": "1.1.0",
    "spatialdata": "0.1.0",
    "scanpy": "1.9.0",
    "anndata": "0.9.0",
    "torch": "2.0.0",
    "numpy": "1.20.0",
    "pandas": "1.5.0",
}


def check_dependencies():
    """Check if all required dependencies are available and compatible"""
    missing = []
    incompatible = []

    for package, min_version in REQUIRED_VERSIONS.items():
        try:
            if package == "scvi":
                import scvi

                current_version = scvi.__version__
            elif package == "spatialdata":
                import spatialdata

                current_version = spatialdata.__version__
            else:
                current_version = pkg_resources.get_distribution(package).version

            if version.parse(current_version) < version.parse(min_version):
                incompatible.append(f"{package}: {current_version} < {min_version}")
            else:
                print(f"✓ {package}: {current_version}")

        except ImportError:
            missing.append(package)
        except Exception as e:
            print(f"Warning: Could not check {package}: {e}")

    if missing:
        print(f"Missing packages: {', '.join(missing)}")
        return False

    if incompatible:
        print(f"Incompatible versions: {', '.join(incompatible)}")
        return False

    return True


if __name__ == "__main__":
    if not check_dependencies():
        sys.exit(1)
    print("All dependencies satisfied!")
