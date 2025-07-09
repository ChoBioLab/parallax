# bin/validate_gpu_params.py
#!/usr/bin/env python3
import json
import os
import sys

import torch


def validate_gpu_configuration():
    """Validate GPU configuration and adjust parameters"""

    # Read GPU info from CHECK_GPU_AVAILABILITY
    try:
        with open("gpu_info.txt", "r") as f:
            gpu_info = json.load(f)
    except FileNotFoundError:
        gpu_info = {"available_gpus": 0, "adjusted_gpus": 0, "cuda_available": False}

    # Get requested parameters
    requested_gpus = (
        int(os.environ.get("CUDA_VISIBLE_DEVICES", "0").count(",") + 1)
        if os.environ.get("CUDA_VISIBLE_DEVICES")
        else 0
    )

    # Validate PyTorch CUDA availability
    torch_cuda_available = torch.cuda.is_available()
    torch_gpu_count = torch.cuda.device_count()

    config = {
        "use_gpu": torch_cuda_available and gpu_info["adjusted_gpus"] > 0,
        "num_gpus": min(gpu_info["adjusted_gpus"], torch_gpu_count)
        if torch_cuda_available
        else 0,
        "device": "cuda"
        if torch_cuda_available and gpu_info["adjusted_gpus"] > 0
        else "cpu",
        "memory_fraction": float(os.environ.get("GPU_MEMORY_FRACTION", "0.8")),
    }

    # Set environment variables for downstream processes
    if config["use_gpu"]:
        os.environ["CUDA_VISIBLE_DEVICES"] = ",".join(
            str(i) for i in range(config["num_gpus"])
        )
        # Configure PyTorch memory allocation
        torch.cuda.set_per_process_memory_fraction(config["memory_fraction"])
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    print(f"GPU Configuration: {config}", file=sys.stderr)
    return config


if __name__ == "__main__":
    config = validate_gpu_configuration()
    with open("gpu_config.json", "w") as f:
        json.dump(config, f, indent=2)
