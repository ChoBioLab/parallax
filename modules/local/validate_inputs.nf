process VALIDATE_INPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.9.6--pyhdfd78af_0' :
        'scverse/scanpy:1.9.6' }"

    input:
    tuple val(meta), path(zarr_path)

    output:
    tuple val(meta), path(zarr_path), emit: validated
    path "validation_report.json", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import sys
    import json
    import os
    from pathlib import Path

    def validate_zarr_input(zarr_path, meta):
        validation_results = {
            'sample_id': meta.get('id', 'unknown'),
            'zarr_path': str(zarr_path),
            'valid': True,
            'errors': [],
            'warnings': []
        }

        try:
            # Check if zarr package is available
            try:
                import zarr
                import pandas as pd
            except ImportError as e:
                validation_results['errors'].append(f"Required package not available: {str(e)}")
                validation_results['valid'] = False
                return validation_results

            # Check if path exists
            zarr_file = Path(zarr_path)
            if not zarr_file.exists():
                validation_results['errors'].append(f"Zarr file does not exist: {zarr_path}")
                validation_results['valid'] = False
                return validation_results

            # Try to open zarr store
            try:
                store = zarr.open(str(zarr_file), mode='r')
            except Exception as e:
                validation_results['errors'].append(f"Cannot open zarr file: {str(e)}")
                validation_results['valid'] = False
                return validation_results

            # Check required groups
            required_groups = ['X', 'obs', 'var']
            for group in required_groups:
                if group not in store:
                    validation_results['errors'].append(f"Missing required group: {group}")
                    validation_results['valid'] = False

            if not validation_results['valid']:
                return validation_results

            # Validate dimensions
            try:
                n_obs = store['X'].shape[0]
                n_vars = store['X'].shape[1]

                if 'obs' in store and hasattr(store['obs'], '__len__'):
                    obs_len = len(store['obs'])
                    if n_obs != obs_len:
                        validation_results['errors'].append(f"Dimension mismatch: X has {n_obs} observations but obs has {obs_len}")
                        validation_results['valid'] = False

                if 'var' in store and hasattr(store['var'], '__len__'):
                    var_len = len(store['var'])
                    if n_vars != var_len:
                        validation_results['errors'].append(f"Dimension mismatch: X has {n_vars} variables but var has {var_len}")
                        validation_results['valid'] = False

                # Add warnings for low counts
                if n_obs < 100:
                    validation_results['warnings'].append(f"Low cell count: {n_obs} cells")

                if n_vars < 1000:
                    validation_results['warnings'].append(f"Low gene count: {n_vars} genes")

                validation_results['stats'] = {
                    'n_obs': n_obs,
                    'n_vars': n_vars,
                    'has_spatial': 'spatial' in store
                }

            except Exception as e:
                validation_results['errors'].append(f"Error validating dimensions: {str(e)}")
                validation_results['valid'] = False

        except Exception as e:
            validation_results['errors'].append(f"Unexpected error during validation: {str(e)}")
            validation_results['valid'] = False

        return validation_results

    # Parse meta from Nextflow
    import ast
    meta_str = '''${meta}'''
    try:
        meta = ast.literal_eval(meta_str)
    except:
        meta = {'id': 'unknown'}

    zarr_path = "${zarr_path}"

    # Run validation
    results = validate_zarr_input(zarr_path, meta)

    # Write results
    with open('validation_report.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Print results
    if results['valid']:
        print(f"✓ Validation passed for {results['sample_id']}")
        if results['warnings']:
            print("Warnings:")
            for warning in results['warnings']:
                print(f"  - {warning}")
    else:
        print(f"✗ Validation failed for {results['sample_id']}")
        for error in results['errors']:
            print(f"  ERROR: {error}")
        sys.exit(1)

    # Create versions file
    python_version = sys.version.split()[0]

    try:
        import zarr
        zarr_version = zarr.__version__
    except:
        zarr_version = "unknown"

    versions_content = f'''"{task.process}":
    python: "{python_version}"
    zarr: "{zarr_version}"
'''

    with open('versions.yml', 'w') as f:
        f.write(versions_content)
    """
}
