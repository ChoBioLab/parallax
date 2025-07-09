# bin/error_handler.py
#!/usr/bin/env python3
import json
import logging
import sys
import traceback
from datetime import datetime
from pathlib import Path


class PipelineErrorHandler:
    def __init__(self, module_name, sample_id):
        self.module_name = module_name
        self.sample_id = sample_id
        self.setup_logging()

    def setup_logging(self):
        """Setup consistent logging format"""
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler(f"{self.module_name}_{self.sample_id}.log"),
                logging.StreamHandler(sys.stderr),
            ],
        )
        self.logger = logging.getLogger(self.module_name)

    def handle_error(self, error, create_placeholder=True):
        """Standardized error handling"""
        error_info = {
            "module": self.module_name,
            "sample_id": self.sample_id,
            "timestamp": datetime.now().isoformat(),
            "error_type": type(error).__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }

        # Log error
        self.logger.error(
            f"Error in {self.module_name} for sample {self.sample_id}: {error}"
        )
        self.logger.debug(f"Full traceback: {traceback.format_exc()}")

        # Write error report
        with open(f"error_report_{self.sample_id}.json", "w") as f:
            json.dump(error_info, f, indent=2)

        # Create placeholder outputs if requested
        if create_placeholder:
            self.create_placeholder_outputs()

        return error_info

    def create_placeholder_outputs(self):
        """Create standardized placeholder outputs"""
        placeholder_info = {
            "status": "failed",
            "module": self.module_name,
            "sample_id": self.sample_id,
            "timestamp": datetime.now().isoformat(),
        }

        # Create placeholder files based on module type
        if self.module_name in ["RESOLVI_TRAIN", "SCVIVA_TRAIN"]:
            # Create empty model placeholder
            Path("model_placeholder.pkl").touch()
            with open("model_info.json", "w") as f:
                json.dump(placeholder_info, f, indent=2)

        elif self.module_name in ["RESOLVI_ANALYZE", "SCVIVA_ANALYZE"]:
            # Create empty results placeholder
            Path("results_placeholder.h5ad").touch()
            with open("analysis_info.json", "w") as f:
                json.dump(placeholder_info, f, indent=2)

        # Always create versions file
        self.create_versions_file()

    def create_versions_file(self):
        """Create standardized versions.yml"""
        try:
            import numpy as np
            import pandas as pd
            import scanpy as sc
            import scvi

            versions = {
                self.module_name: {
                    "scanpy": sc.__version__,
                    "scvi-tools": scvi.__version__,
                    "pandas": pd.__version__,
                    "numpy": np.__version__,
                    "python": sys.version.split()[0],
                }
            }
        except ImportError as e:
            versions = {
                self.module_name: {
                    "error": f"Could not determine versions: {e}",
                    "python": sys.version.split()[0],
                }
            }

        with open("versions.yml", "w") as f:
            import yaml

            yaml.dump(versions, f, default_flow_style=False)


# Usage in modules:
def safe_execution(func, module_name, sample_id, *args, **kwargs):
    """Wrapper for safe execution with standardized error handling"""
    handler = PipelineErrorHandler(module_name, sample_id)

    try:
        result = func(*args, **kwargs)
        handler.logger.info(f"Successfully completed {module_name} for {sample_id}")
        return result
    except Exception as e:
        error_info = handler.handle_error(e, create_placeholder=True)
        # Re-raise with additional context
        raise RuntimeError(
            f"Module {module_name} failed for sample {sample_id}: {e}"
        ) from e
