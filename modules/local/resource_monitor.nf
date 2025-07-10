process RESOURCE_MONITOR {
    tag "monitor_${stage}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--slim' :
        'python:3.11-slim' }"

    input:
    val(stage)

    output:
    path "resource_usage_${stage}.json", emit: report

    script:
    """
    #!/usr/bin/env python3
    import json
    import psutil
    import subprocess
    from datetime import datetime

    def get_system_resources():
        resources = {
            'timestamp': datetime.now().isoformat(),
            'stage': '${stage}',
            'cpu_count': psutil.cpu_count(),
            'cpu_percent': psutil.cpu_percent(interval=1),
            'memory_total_gb': psutil.virtual_memory().total / (1024**3),
            'memory_available_gb': psutil.virtual_memory().available / (1024**3),
            'memory_percent': psutil.virtual_memory().percent
        }

        # GPU information
        try:
            result = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total,memory.used,utilization.gpu', 
                                   '--format=csv,noheader,nounits'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                gpu_lines = result.stdout.strip().split('\\n')
                resources['gpus'] = []
                for i, line in enumerate(gpu_lines):
                    if line.strip():
                        parts = line.split(', ')
                        if len(parts) >= 4:
                            resources['gpus'].append({
                                'id': i,
                                'name': parts[0],
                                'memory_total_mb': int(parts[1]),
                                'memory_used_mb': int(parts[2]),
                                'utilization_percent': int(parts[3])
                            })
        except:
            resources['gpus'] = []

        return resources

    resources = get_system_resources()

    with open('resource_usage_${stage}.json', 'w') as f:
        json.dump(resources, f, indent=2)

    print(f"Resource monitoring for stage: ${stage}")
    print(f"CPU: {resources['cpu_percent']}%")
    print(f"Memory: {resources['memory_percent']}%")
    if resources['gpus']:
        for gpu in resources['gpus']:
            print(f"GPU {gpu['id']}: {gpu['utilization_percent']}% util, {gpu['memory_used_mb']}/{gpu['memory_total_mb']} MB")
    """
}

