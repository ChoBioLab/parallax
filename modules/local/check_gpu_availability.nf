process CHECK_GPU_AVAILABILITY {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--slim' :
        'python:3.11-slim' }"

    output:
    path "gpu_info.txt", emit: gpu_info
    path "versions.yml", emit: versions

    script:
    """
    python ${projectDir}/bin/check_gpu_availability.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

