process VALIDATE_INPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_path)

    output:
    tuple val(meta), path(zarr_path), emit: validated
    path "validation_report.json", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import sys
    import json
    from pathlib import Path

    def validate_spatialdata_input(zarr_path, meta):
        validation_results = {
            'sample_id': meta.get('id', 'unknown'),
            'zarr_path': str(zarr_path),
            'valid': True,
            'errors': [],
            'warnings': []
        }

        try:
            import spatialdata as sd
            import anndata as ad
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

        try:
            # Load SpatialData object
            sdata = sd.read_zarr(str(zarr_file))

            # Check if there's a table
            if not hasattr(sdata, 'tables') or len(sdata.tables) == 0:
                validation_results['errors'].append("No tables found in SpatialData object")
                validation_results['valid'] = False
                return validation_results

            # Get the main table (usually 'table')
            table_key = 'table' if 'table' in sdata.tables else list(sdata.tables.keys())[0]
            adata = sdata.tables[table_key]

            # Validate AnnData structure
            if adata.X is None:
                validation_results['errors'].append("Expression matrix (X) is None")
                validation_results['valid'] = False
            if adata.obs is None or len(adata.obs) == 0:
                validation_results['errors'].append("Cell metadata (obs) is empty")
                validation_results['valid'] = False
            if adata.var is None or len(adata.var) == 0:
                validation_results['errors'].append("Gene metadata (var) is empty")
                validation_results['valid'] = False

            if validation_results['valid']:
                # Get basic stats
                n_obs, n_vars = adata.shape

                validation_results['stats'] = {
                    'n_obs': n_obs,
                    'n_vars': n_vars,
                    'table_key': table_key,
                    'has_spatial_info': len(sdata.shapes) > 0 or len(sdata.points) > 0
                }

                # Add warnings for low counts
                if n_obs < 100:
                    validation_results['warnings'].append(f"Low cell count: {n_obs} cells")
                if n_vars < 100:
                    validation_results['warnings'].append(f"Low gene count: {n_vars} genes")

        except Exception as e:
            validation_results['errors'].append(f"Error reading SpatialData: {str(e)}")
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
    results = validate_spatialdata_input(zarr_path, meta)

    # Write results
    with open('validation_report.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Print results
    if results['valid']:
        print(f"✓ Validation passed for {results['sample_id']}")
        if 'stats' in results:
            stats = results['stats']
            print(f"  Table: {stats['table_key']}")
            print(f"  Cells: {stats['n_obs']}, Genes: {stats['n_vars']}")
            print(f"  Has spatial info: {stats['has_spatial_info']}")
        if results['warnings']:
            print("  Warnings:")
            for warning in results['warnings']:
                print(f"    - {warning}")
    else:
        print(f"✗ Validation failed for {results['sample_id']}")
        for error in results['errors']:
            print(f"  ERROR: {error}")
        sys.exit(1)

    # Create versions file
    python_version = sys.version.split()[0]
    try:
        import spatialdata as sd
        sd_version = sd.__version__
    except:
        sd_version = "unknown"

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{python_version}"\\n')
        f.write(f'    spatialdata: "{sd_version}"\\n')
    """
}
