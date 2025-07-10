process ERROR_HANDLER {
    tag "$meta.id"
    label 'process_single'

    publishDir "${params.outdir}/errors", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(error_type), val(error_message)

    output:
    path "error_report_${meta.id}.json", emit: report
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import sys
    from datetime import datetime

    error_report = {
        'sample_id': '${meta.id}',
        'error_type': '${error_type}',
        'error_message': '''${error_message}''',
        'timestamp': datetime.now().isoformat(),
        'pipeline_stage': 'validation',
        'meta': ${meta}
    }

    with open('error_report_${meta.id}.json', 'w') as f:
        json.dump(error_report, f, indent=2)

    print(f"ERROR: {error_report['error_type']} for sample {error_report['sample_id']}")
    print(f"Message: {error_report['error_message']}")

    # Create versions file
    python_version = sys.version.split()[0]
    versions_content = f'''"{task.process}":
    python: "{python_version}"
'''

    with open('versions.yml', 'w') as f:
        f.write(versions_content)
    """
}

