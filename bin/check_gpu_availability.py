#!/usr/bin/env python3
import json
import subprocess
import sys


def check_gpu_availability():
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        gpus = []
        for line in result.stdout.strip().split("\n"):
            if line.strip():
                idx, name, memory = line.split(", ")
                gpus.append(
                    {"index": int(idx), "name": name.strip(), "memory_mb": int(memory)}
                )

        gpu_info = {"available_gpus": len(gpus), "gpus": gpus, "cuda_available": True}

    except (subprocess.CalledProcessError, FileNotFoundError):
        gpu_info = {"available_gpus": 0, "gpus": [], "cuda_available": False}

    with open("gpu_info.txt", "w") as f:
        json.dump(gpu_info, f, indent=2)

    print(f"GPU Check: {gpu_info['available_gpus']} GPUs available")
    return gpu_info


if __name__ == "__main__":
    check_gpu_availability()
