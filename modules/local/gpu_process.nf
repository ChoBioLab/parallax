process GPU_PROCESS {
    tag "$meta.id"
    label 'process_gpu'

    conda "${projectDir}/environment.yml"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}_output.*"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import sys
    import os
    import json
    import torch
    from pathlib import Path

    def validate_gpu_setup():
        \"\"\"Validate GPU is available and functional - exit if not\"\"\"
        if not torch.cuda.is_available():
            print("ERROR: CUDA is not available. This workflow requires GPU.")
            sys.exit(1)

        try:
            device = torch.device('cuda')
            test_tensor = torch.randn(100, 100).to(device)
            _ = test_tensor @ test_tensor.T
            torch.cuda.synchronize()

            print(f"✓ GPU validated: {torch.cuda.get_device_name()}")
            print(f"✓ CUDA version: {torch.version.cuda}")
            print(f"✓ GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

            return device
        except Exception as e:
            print(f"ERROR: GPU validation failed: {str(e)}")
            sys.exit(1)

    def main():
        # Validate GPU setup (will exit if GPU not available)
        device = validate_gpu_setup()

        # Set GPU as default device
        torch.cuda.set_device(0)

        # Your GPU-accelerated processing logic here
        print(f"Processing ${meta.id} on GPU")

        # Example: Load data to GPU
        # data = torch.load('${input_file}').to(device)

        # Create output
        output_file = Path("${prefix}_output.json")
        result = {
            'sample_id': '${meta.id}',
            'device': str(device),
            'gpu_name': torch.cuda.get_device_name(),
            'cuda_version': torch.version.cuda,
            'status': 'completed'
        }

        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)

    if __name__ == "__main__":
        main()

    # Create versions file
    python_version = sys.version.split()[0]
    torch_version = torch.__version__
    cuda_version = torch.version.cuda

    versions_content = f'''"{task.process}":
    python: "{python_version}"
    pytorch: "{torch_version}"
    cuda: "{cuda_version}"
    gpu_name: "{torch.cuda.get_device_name()}"
'''

    with open('versions.yml', 'w') as f:
        f.write(versions_content)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "ERROR: GPU required for this workflow" >&2
    exit 1
    """
}

