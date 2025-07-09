# bin/module_wrapper.py
#!/usr/bin/env python3
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import yaml


class ModuleWrapper:
    def __init__(self, module_name, sample_id):
        self.module_name = module_name
        self.sample_id = sample_id
        self.start_time = datetime.now()
        self.setup_logging()

    def setup_logging(self):
        """Setup consistent logging"""
        log_format = "%(asctime)s | %(name)s | %(levelname)s | %(message)s"
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(f"{self.module_name}_{self.sample_id}.log"),
                logging.StreamHandler(sys.stderr),
            ],
        )
        self.logger = logging.getLogger(f"nf-core.resolvinf.{self.module_name}")

    def log_start(self):
        """Log module start with consistent format"""
        self.logger.info(f"Starting {self.module_name} for sample {self.sample_id}")
        self.logger.info(f"Timestamp: {self.start_time.isoformat()}")

    def log_completion(self):
        """Log successful completion"""
        end_time = datetime.now()
        duration = end_time - self.start_time
        self.logger.info(f"Completed {self.module_name} for sample {self.sample_id}")
        self.logger.info(f"Duration: {duration.total_seconds():.2f} seconds")

    def create_versions_file(self):
        """Create comprehensive versions.yml"""
        try:
            versions = {self.module_name: self._get_package_versions()}

            with open("versions.yml", "w") as f:
                yaml.dump(versions, f, default_flow_style=False, sort_keys=True)

        except Exception as e:
            self.logger.error(f"Failed to create versions file: {e}")

    def _get_package_versions(self):
        """Get versions of all relevant packages"""
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
        ]

        versions = {"python": sys.version.split()[0]}

        for package in packages_to_check:
            try:
                versions[package] = pkg_resources.get_distribution(package).version
            except pkg_resources.DistributionNotFound:
                versions[package] = "not_installed"
            except Exception:
                versions[package] = "unknown"

        return versions


# Usage in modules:
if __name__ == "__main__":
    wrapper = ModuleWrapper("${task.process}", "${meta.id}")
    wrapper.log_start()

    try:
        # Your module code here
        wrapper.log_completion()
        wrapper.create_versions_file()
    except Exception as e:
        wrapper.logger.error(f"Module failed: {e}")
        wrapper.create_versions_file()
        raise
