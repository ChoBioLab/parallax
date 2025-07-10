process VALIDATE_GPU {
    tag "gpu_check"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pytorch:2.0.1--py3.11_cuda11.8_cudnn8.7.0_0' :
        'pytorch/pytorch:2.0.1-cuda11.8-cudnn8-runtime' }"

    output:
    path "gpu_info.json", emit: info
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    import json
    import torch

    def validate_gpu_environment():
        gpu_info = {
            'cuda_available': torch.cuda.is_available(),
            'device_count': torch.cuda.device_count() if torch.cuda.is_available() else 0,
            'cuda_version': torch.version.cuda,
            'pytorch_version': torch.__version__
        }

        if not gpu_info['cuda_available']:
            print("ERROR: CUDA is not available. This workflow requires GPU support.")
            sys.exit(1)

        if gpu_info['device_count'] == 0:
            print("ERROR: No GPU devices found.")
            sys.exit(1)

        # Get detailed GPU info
        for i in range(gpu_info['device_count']):
            props = torch.cuda.get_device_properties(i)
            gpu_info[f'gpu_{i}'] = {
                'name': props.name,
                'memory_total': props.total_memory,
                'memory_total_gb': props.total_memory / 1e9,
                'compute_capability': f"{props.major}.{props.minor}"
            }

        print(f"✓ GPU validation successful")
        print(f"✓ Found {gpu_info['device_count']} GPU(s)")
        print(f"✓ Primary GPU: {gpu_info['gpu_0']['name']}")
        print(f"✓ GPU Memory: {gpu_info['gpu_0']['memory_total_gb']:.1f} GB")

        return gpu_info

    # Run validation
    gpu_info = validate_gpu_environment()

    # Save GPU info
    with open('gpu_info.json', 'w') as f:
        json.dump(gpu_info, f, indent=2)

    # Create versions file
    versions_content = f'''"{task.process}":
    pytorch: "{torch.__version__}"
    cuda: "{torch.version.cuda}"
    gpu_name: "{gpu_info['gpu_0']['name']}"
    gpu_memory_gb: "{gpu_info['gpu_0']['memory_total_gb']:.1f}"
'''

    with open('versions.yml', 'w') as f:
        f.write(versions_content)
    """
}

