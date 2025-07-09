process VALIDATE_INPUTS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(zarr_path)

    output:
    tuple val(meta), path(zarr_path), emit: validated
    path "validation_report.json", emit: report
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    import json
    import zarr
    import pandas as pd
    from pathlib import Path

    def validate_zarr_input(zarr_path, meta):
        validation_results = {
            'sample_id': meta['id'],
            'zarr_path': str(zarr_path),
            'valid': True,
            'errors': [],
            'warnings': []
        }

        try:
            if not Path(zarr_path).exists():
                validation_results['errors'].append(f"Zarr file does not exist: {zarr_path}")
                validation_results['valid'] = False
                return validation_results

            store = zarr.open(zarr_path, mode='r')

            required_groups = ['X', 'obs', 'var']
            for group in required_groups:
                if group not in store:
                    validation_results['errors'].append(f"Missing required group: {group}")
                    validation_results['valid'] = False

            if not validation_results['valid']:
                return validation_results

            n_obs = store['X'].shape[0]
            n_vars = store['X'].shape[1]

            if n_obs != len(store['obs']):
                validation_results['errors'].append(f"Dimension mismatch: X has {n_obs} observations but obs has {len(store['obs'])}")
                validation_results['valid'] = False

            if n_vars != len(store['var']):
                validation_results['errors'].append(f"Dimension mismatch: X has {n_vars} variables but var has {len(store['var'])}")
                validation_results['valid'] = False

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
            validation_results['errors'].append(f"Error reading zarr file: {str(e)}")
            validation_results['valid'] = False

        return validation_results

    # Main validation
    meta = ${meta}
    zarr_path = "${zarr_path}"

    results = validate_zarr_input(zarr_path, meta)

    with open('validation_report.json', 'w') as f:
        json.dump(results, f, indent=2)

    if results['valid']:
        print(f"✓ Validation passed for {meta['id']}")
        if results['warnings']:
            print("Warnings:")
            for warning in results['warnings']:
                print(f"  - {warning}")
    else:
        print(f"✗ Validation failed for {meta['id']}")
        for error in results['errors']:
            print(f"  ERROR: {error}")
        sys.exit(1)

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
    END_VERSIONS
    """
}

