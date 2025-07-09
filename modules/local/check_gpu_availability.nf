process CHECK_GPU_AVAILABILITY {
    label 'process_single'

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

