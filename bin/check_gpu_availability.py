#!/usr/bin/env python3
"""
Comprehensive GPU availability checker for nf-core/parallax
Handles multiple container environments and provides detailed diagnostics
"""

import json
import os
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional


class GPUChecker:
    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.gpu_info = {
            "available_gpus": 0,
            "gpus": [],
            "cuda_available": False,
            "nvidia_driver_version": None,
            "cuda_version": None,
            "environment_info": self._get_environment_info(),
            "adjusted_gpus": 0,
            "recommendations": [],
        }

    def _get_environment_info(self) -> Dict:
        """Get environment information relevant to GPU usage"""
        return {
            "container_type": self._detect_container_type(),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "nvidia_visible_devices": os.environ.get("NVIDIA_VISIBLE_DEVICES"),
            "singularity_nv": os.environ.get("SINGULARITY_NV"),
            "docker_runtime": self._detect_docker_runtime(),
        }

    def _detect_container_type(self) -> Optional[str]:
        """Detect container type"""
        if os.path.exists("/.singularity.d") or os.environ.get("SINGULARITY_CONTAINER"):
            return "singularity"
        elif os.path.exists("/.dockerenv"):
            return "docker"
        elif os.environ.get("APPTAINER_CONTAINER"):
            return "apptainer"
        return None

    def _detect_docker_runtime(self) -> Optional[str]:
        """Detect Docker runtime for GPU support"""
        try:
            result = subprocess.run(
                ["docker", "info", "--format", "{{.Runtimes}}"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if "nvidia" in result.stdout:
                return "nvidia"
        except (
            subprocess.TimeoutExpired,
            FileNotFoundError,
            subprocess.CalledProcessError,
        ):
            pass
        return None

    def check_nvidia_smi(self) -> bool:
        """Check nvidia-smi availability and parse output"""
        try:
            # Check nvidia-smi availability
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=index,name,memory.total,driver_version,cuda_version",
                    "--format=csv,noheader,nounits",
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )

            gpus = []
            for line in result.stdout.strip().split("\n"):
                if line.strip():
                    parts = [part.strip() for part in line.split(", ")]
                    if len(parts) >= 3:
                        gpu_info = {
                            "index": int(parts[0]),
                            "name": parts[1],
                            "memory_mb": int(parts[2]),
                        }

                        # Driver and CUDA version (may not always be available)
                        if len(parts) > 3 and parts[3] != "N/A":
                            self.gpu_info["nvidia_driver_version"] = parts[3]
                        if len(parts) > 4 and parts[4] != "N/A":
                            self.gpu_info["cuda_version"] = parts[4]

                        gpus.append(gpu_info)

            self.gpu_info["gpus"] = gpus
            self.gpu_info["available_gpus"] = len(gpus)

            print(f"✅ nvidia-smi detected {len(gpus)} GPU(s)", file=sys.stderr)
            return True

        except subprocess.TimeoutExpired:
            print("⚠️  nvidia-smi timeout - GPU may be busy", file=sys.stderr)
            return False
        except subprocess.CalledProcessError as e:
            print(f"⚠️  nvidia-smi failed: {e}", file=sys.stderr)
            return False
        except FileNotFoundError:
            print("❌ nvidia-smi not found", file=sys.stderr)
            return False
        except Exception as e:
            print(f"❌ nvidia-smi check failed: {e}", file=sys.stderr)
            return False

    def check_pytorch_cuda(self) -> bool:
        """Check PyTorch CUDA availability"""
        try:
            import torch

            cuda_available = torch.cuda.is_available()
            device_count = torch.cuda.device_count()

            if cuda_available:
                print(
                    f"✅ PyTorch CUDA available with {device_count} device(s)",
                    file=sys.stderr,
                )
                self.gpu_info["pytorch_cuda_available"] = True
                self.gpu_info["pytorch_device_count"] = device_count

                # Get PyTorch CUDA version
                if hasattr(torch.version, "cuda") and torch.version.cuda:
                    self.gpu_info["pytorch_cuda_version"] = torch.version.cuda

                return True
            else:
                print("❌ PyTorch CUDA not available", file=sys.stderr)
                self.gpu_info["pytorch_cuda_available"] = False
                return False

        except ImportError:
            print("⚠️  PyTorch not available for CUDA check", file=sys.stderr)
            self.gpu_info["pytorch_cuda_available"] = False
            return False
        except Exception as e:
            print(f"❌ PyTorch CUDA check failed: {e}", file=sys.stderr)
            self.gpu_info["pytorch_cuda_available"] = False
            return False

    def adjust_gpu_count(self) -> int:
        """Adjust GPU count based on environment constraints"""
        available = self.gpu_info["available_gpus"]

        # Check CUDA_VISIBLE_DEVICES constraint
        cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
        if cuda_visible is not None:
            if cuda_visible == "":
                adjusted = 0
                print(
                    "🚫 CUDA_VISIBLE_DEVICES is empty - forcing CPU mode",
                    file=sys.stderr,
                )
            else:
                try:
                    visible_devices = [
                        int(x.strip())
                        for x in cuda_visible.split(",")
                        if x.strip().isdigit()
                    ]
                    adjusted = min(available, len(visible_devices))
                    print(
                        f"🎯 CUDA_VISIBLE_DEVICES limits to {adjusted} GPU(s)",
                        file=sys.stderr,
                    )
                except ValueError:
                    adjusted = available
                    print(
                        f"⚠️  Invalid CUDA_VISIBLE_DEVICES format, using all {available} GPU(s)",
                        file=sys.stderr,
                    )
        else:
            adjusted = available

        # Check container GPU support
        container_type = self.gpu_info["environment_info"]["container_type"]
        if container_type and adjusted > 0:
            if container_type == "singularity" and not os.environ.get("SINGULARITY_NV"):
                print(
                    "⚠️  Singularity container without --nv flag may not support GPU",
                    file=sys.stderr,
                )
                self.gpu_info["recommendations"].append(
                    "Use --nv flag with Singularity for GPU support"
                )
            elif (
                container_type == "docker"
                and not self.gpu_info["environment_info"]["docker_runtime"]
            ):
                print(
                    "⚠️  Docker container may not have GPU runtime configured",
                    file=sys.stderr,
                )
                self.gpu_info["recommendations"].append(
                    "Use --gpus flag with Docker for GPU support"
                )

        self.gpu_info["adjusted_gpus"] = adjusted
        return adjusted

    def generate_recommendations(self):
        """Generate GPU setup recommendations"""
        if self.gpu_info["available_gpus"] == 0:
            self.gpu_info["recommendations"].extend(
                [
                    "No GPUs detected - pipeline will run in CPU mode",
                    "For GPU acceleration, ensure NVIDIA drivers and CUDA are installed",
                    "Check that nvidia-smi command works",
                ]
            )
        elif self.gpu_info["adjusted_gpus"] < self.gpu_info["available_gpus"]:
            self.gpu_info["recommendations"].append(
                f"Using {self.gpu_info['adjusted_gpus']}/{self.gpu_info['available_gpus']} available GPUs due to environment constraints"
            )

        # PyTorch-specific recommendations
        if self.gpu_info["available_gpus"] > 0 and not self.gpu_info.get(
            "pytorch_cuda_available", False
        ):
            self.gpu_info["recommendations"].extend(
                [
                    "GPUs detected but PyTorch CUDA not available",
                    "Install PyTorch with CUDA support for GPU acceleration",
                ]
            )

    def save_results(self):
        """Save GPU information to files"""
        # Detailed JSON output
        json_file = self.output_dir / "gpu_info.json"
        with open(json_file, "w") as f:
            json.dump(self.gpu_info, f, indent=2)

        # Simple text output for backward compatibility
        txt_file = self.output_dir / "gpu_info.txt"
        with open(txt_file, "w") as f:
            json.dump(
                {
                    "available_gpus": self.gpu_info["available_gpus"],
                    "adjusted_gpus": self.gpu_info["adjusted_gpus"],
                    "cuda_available": self.gpu_info["cuda_available"],
                },
                f,
            )

        # Status file for Nextflow
        status_file = self.output_dir / "gpu_status.txt"
        with open(status_file, "w") as f:
            f.write(f"{self.gpu_info['adjusted_gpus']}")

    def print_summary(self):
        """Print comprehensive GPU summary"""
        print("\n🎮 GPU Check Summary:", file=sys.stderr)
        print(f"   Available GPUs: {self.gpu_info['available_gpus']}", file=sys.stderr)
        print(f"   Usable GPUs: {self.gpu_info['adjusted_gpus']}", file=sys.stderr)
        print(f"   CUDA Available: {self.gpu_info['cuda_available']}", file=sys.stderr)

        if self.gpu_info["nvidia_driver_version"]:
            print(
                f"   NVIDIA Driver: {self.gpu_info['nvidia_driver_version']}",
                file=sys.stderr,
            )

        if self.gpu_info["cuda_version"]:
            print(f"   CUDA Version: {self.gpu_info['cuda_version']}", file=sys.stderr)

        if self.gpu_info["recommendations"]:
            print("\n💡 Recommendations:", file=sys.stderr)
            for rec in self.gpu_info["recommendations"]:
                print(f"   {rec}", file=sys.stderr)


def main():
    """Main execution function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Check GPU availability for nf-core/parallax"
    )
    parser.add_argument(
        "--output-dir", default=".", help="Output directory for results"
    )
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    args = parser.parse_args()

    if not args.quiet:
        print("🎮 nf-core/parallax GPU Checker", file=sys.stderr)
        print("=" * 40, file=sys.stderr)

    checker = GPUChecker(args.output_dir)

    try:
        # Check nvidia-smi
        nvidia_available = checker.check_nvidia_smi()

        # Check PyTorch CUDA
        pytorch_cuda = checker.check_pytorch_cuda()

        # Set overall CUDA availability
        checker.gpu_info["cuda_available"] = nvidia_available and pytorch_cuda

        # Adjust GPU count based on constraints
        adjusted_gpus = checker.adjust_gpu_count()

        # Generate recommendations
        checker.generate_recommendations()

        # Save results
        checker.save_results()

        if not args.quiet:
            checker.print_summary()

        print(
            f"GPU Check Complete: {adjusted_gpus} GPU(s) available for use",
            file=sys.stderr,
        )
        return 0

    except Exception as e:
        print(f"💥 GPU check failed: {e}", file=sys.stderr)
        checker.gpu_info["error"] = str(e)
        checker.save_results()
        return 1


if __name__ == "__main__":
    sys.exit(main())
