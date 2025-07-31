process RESOLVI_VISUALIZE {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${projectDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'parallax.sif' :
        'parallax:latest' }"

    containerOptions '--writable-tmpfs'
    
    input:
    tuple val(meta), path(adata)
    val marker_genes
    
    output:
    tuple val(meta), path("plots/*"), emit: plots
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_arg = marker_genes ? marker_genes.join(',') : ""
    """
    #!/usr/bin/env python3

    # Import shared environment setup
    import sys
    sys.path.insert(0, '${projectDir}/bin')

    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import logging
    import warnings
    import subprocess
    import sys
    import os
    import tempfile
    from pathlib import Path
    
    # Set up container-safe numba caching
    if os.path.exists('/.singularity.d') or os.environ.get('SINGULARITY_CONTAINER'):
        # Inside Singularity - use writable temp location
        numba_cache = tempfile.mkdtemp(prefix='numba_cache_')
        os.environ['NUMBA_CACHE_DIR'] = numba_cache
        # Keep JIT enabled but use safe cache location
    else:
        # Regular environment - use default caching
        pass
    
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Suppress warnings
    warnings.filterwarnings("ignore")
    plt.ioff()  # Turn off interactive plotting
    
    logger.info("=== Starting ResolVI Visualization ===")
    
    # Create plots directory
    os.makedirs("plots", exist_ok=True)
    
    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad("${adata}")
    
    # Define marker genes
    marker_genes_list = "${genes_arg}".split(',') if "${genes_arg}" else []
    
    # Default marker genes if none provided
    if not marker_genes_list or marker_genes_list == ['']:
        marker_genes_list = [
            "CYP2A6", "ALDOB", "CLU", "APOE", "MST1", "CXCL2", "PROX1", "ADH1B", "SDS", "FBLN1",
            "MMP2", "CCDC80", "COL1A1", "COL3A1", "PDGFRA", "PODN", "PDGFRB", "VCAN", "MXRA8", "THY1"
        ]
    
    # Filter available genes
    available_markers = [gene for gene in marker_genes_list if gene in adata.var_names]
    logger.info(f"Found {len(available_markers)} of {len(marker_genes_list)} marker genes in dataset")
    
    # Custom spatial plotting function
    def xenium_spatial_plot(adata, color, layer=None, title=None, figsize=(8, 7), 
                           point_size=1, save_path=None, colormap="viridis", categorical=False):
        coords = adata.obsm["X_spatial"]
        
        if color in adata.obs_keys():
            values = adata.obs[color]
            if not categorical:
                categorical = True
        else:
            if layer is not None and layer in adata.layers:
                expr_matrix = adata.layers[layer]
            else:
                expr_matrix = adata.X
            
            if color not in adata.var_names:
                raise ValueError(f"Gene '{color}' not found in adata.var_names.")
            
            gene_loc = adata.var_names.get_loc(color)
            if hasattr(expr_matrix, "toarray"):
                values = expr_matrix[:, gene_loc].toarray().flatten()
            else:
                values = expr_matrix[:, gene_loc]
        
        fig, ax = plt.subplots(figsize=figsize)
        
        if categorical:
            values_series = pd.Series(values) if not isinstance(values, pd.Series) else values
            
            if pd.api.types.is_categorical_dtype(values_series.dtype):
                values_series = values_series.astype(str)
                categories = sorted(values_series.unique().tolist())
            else:
                try:
                    categories = sorted(values_series.unique().tolist(), key=lambda x: (isinstance(x, str), str(x)))
                except TypeError:
                    categories = values_series.unique().tolist()
            
            n_categories = len(categories)
            
            if n_categories == 0:
                ax.scatter(coords[:, 0], coords[:, 1], s=point_size, color="gray", alpha=0.7)
            else:
                if n_categories <= 10:
                    cmap_colors = plt.cm.get_cmap("tab10").colors
                elif n_categories <= 20:
                    cmap_colors = plt.cm.get_cmap("tab20").colors
                else:
                    cmap_colors = [plt.cm.get_cmap("viridis")(i) for i in np.linspace(0, 1, n_categories)]
                
                color_dict = {cat: cmap_colors[i % len(cmap_colors)] for i, cat in enumerate(categories)}
                point_colors = [color_dict.get(str(val), "lightgrey") for val in values_series]
                
                ax.scatter(coords[:, 0], coords[:, 1], s=point_size, c=point_colors, alpha=0.7)
                
                if n_categories <= 20 and n_categories > 0:
                    legend_elements = [
                        plt.Line2D([0], [0], marker="o", color="w", label=str(cat),
                                 markerfacecolor=color_dict.get(cat, "lightgrey"),
                                 markersize=max(5, np.sqrt(point_size * 10)), linestyle="None")
                        for cat in categories
                    ]
                    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1, 0.5),
                             frameon=False, fontsize="small")
        else:
            scatter = ax.scatter(coords[:, 0], coords[:, 1], s=point_size, c=values, 
                               cmap=colormap, alpha=0.7)
            plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.05)
        
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title if title else color)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            return fig, ax
    
    # Generate UMAP comparison plots
    if "X_umap" in adata.obsm:
        logger.info("Generating UMAP comparison visualizations...")
        try:
            fig, axes = plt.subplots(1, 2, figsize=(16, 7))
            
            if "annotation" in adata.obs:
                sc.pl.umap(adata, color="annotation", ax=axes[0], show=False, frameon=False)
                axes[0].set_title("Original Annotation")
            
            if "resolvi_predicted" in adata.obs:
                sc.pl.umap(adata, color="resolvi_predicted", ax=axes[1], show=False, frameon=False)
                axes[1].set_title("ResolVI Predicted")
            
            plt.tight_layout()
            plt.savefig("plots/umap_annotation_comparison.png", dpi=300, bbox_inches="tight")
            plt.close()
            logger.info("Saved UMAP annotation comparison")
            
            if "annotation_match" in adata.obs:
                sc.pl.umap(adata, color="annotation_match", show=False)
                plt.savefig("plots/umap_annotation_match.png", dpi=300, bbox_inches="tight")
                plt.close()
                logger.info("Saved UMAP annotation match")
                
        except Exception as e:
            logger.error(f"Error generating UMAP plots: {e}")
    
    # Generate confusion matrix
    if "annotation" in adata.obs and "resolvi_predicted" in adata.obs:
        logger.info("Generating confusion matrix...")
        try:
            confusion = pd.crosstab(
                adata.obs["annotation"],
                adata.obs["resolvi_predicted"],
                normalize="index",
            )
            plt.figure(figsize=(12, 10))
            sns.heatmap(
                confusion,
                cmap="viridis",
                annot=True,
                fmt=".2f",
                cbar_kws={"label": "Proportion"},
            )
            plt.title("Original vs ResolVI Cell Type Assignment")
            plt.ylabel("Original Annotation")
            plt.xlabel("ResolVI Predicted")
            plt.tight_layout()
            plt.savefig("plots/annotation_confusion_matrix.png", dpi=300, bbox_inches="tight")
            plt.close()
            logger.info("Saved confusion matrix")
        except Exception as e:
            logger.error(f"Error generating confusion matrix: {e}")
    
    # Generate spatial plots
    logger.info("Generating spatial plots...")
    for category in ["annotation", "resolvi_predicted", "annotation_match"]:
        if category in adata.obs:
            try:
                xenium_spatial_plot(
                    adata,
                    color=category,
                    title=f"Spatial distribution of {category}",
                    save_path=f"plots/spatial_{category}.png",
                    point_size=1,
                    categorical=True,
                )
                logger.info(f"Generated spatial plot for {category}")
            except Exception as e:
                logger.error(f"Error plotting category '{category}': {e}")
    
    # Generate spatial plots for noise components
    logger.info("Generating spatial plots for noise components...")
    for component in ["true_proportion", "diffusion_proportion", "background_proportion"]:
        if component in adata.obs:
            try:
                xenium_spatial_plot(
                    adata,
                    color=component,
                    title=f"Spatial distribution of {component}",
                    save_path=f"plots/spatial_{component}.png",
                    point_size=1,
                    colormap="viridis",
                    categorical=False,
                )
                logger.info(f"Generated spatial plot for {component}")
            except Exception as e:
                logger.error(f"Error plotting component '{component}': {e}")
    
    # Generate gene expression comparison plots
    logger.info("Generating gene expression comparison plots...")
    num_genes_to_plot = min(20, len(available_markers))
    
    for gene in available_markers[:num_genes_to_plot]:
        try:
            fig, axs = plt.subplots(1, 2, figsize=(16, 7))
            coords = adata.obsm["X_spatial"]
            gene_loc = adata.var_names.get_loc(gene)
            
            # Original counts
            if "counts" in adata.layers:
                if hasattr(adata.layers["counts"], "toarray"):
                    values_orig = adata.layers["counts"][:, gene_loc].toarray().flatten()
                else:
                    values_orig = adata.layers["counts"][:, gene_loc]
            else:
                if hasattr(adata.X, "toarray"):
                    values_orig = adata.X[:, gene_loc].toarray().flatten()
                else:
                    values_orig = adata.X[:, gene_loc]
            
            scatter1 = axs[0].scatter(coords[:, 0], coords[:, 1], s=1.0, c=values_orig, 
                                    cmap="viridis", alpha=0.7)
            axs[0].set_title(f"{gene} - Original")
            axs[0].set_xticks([])
            axs[0].set_yticks([])
            plt.colorbar(scatter1, ax=axs[0], shrink=0.5)
            
            # Corrected counts
            if "resolvi_corrected" in adata.layers:
                if hasattr(adata.layers["resolvi_corrected"], "toarray"):
                    values_corr = adata.layers["resolvi_corrected"][:, gene_loc].toarray().flatten()
                else:
                    values_corr = adata.layers["resolvi_corrected"][:, gene_loc]
                
                scatter2 = axs[1].scatter(coords[:, 0], coords[:, 1], s=1.0, c=values_corr, 
                                        cmap="viridis", alpha=0.7)
                axs[1].set_title(f"{gene} - Corrected (ResolVI)")
                plt.colorbar(scatter2, ax=axs[1], shrink=0.5)
            else:
                axs[1].set_title(f"{gene} - Corrected (N/A)")
                axs[1].text(0.5, 0.5, "Corrected counts not available", ha="center", va="center", 
                           transform=axs[1].transAxes)
            
            axs[1].set_xticks([])
            axs[1].set_yticks([])
            plt.tight_layout()
            
            plt.savefig(f"plots/gene_{gene}_comparison.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            logger.info(f"Generated comparison plot for gene {gene}")
            
        except Exception as e:
            logger.error(f"Error plotting gene '{gene}': {e}")
    
    # Generate overview spatial plot
    try:
        if "resolvi_predicted" in adata.obs:
            xenium_spatial_plot(
                adata,
                color="resolvi_predicted",
                title="Spatial Overview - ResolVI Predictions",
                save_path="plots/spatial_overview.png",
                point_size=1,
                categorical=True,
            )
        elif "annotation" in adata.obs:
            xenium_spatial_plot(
                adata,
                color="annotation",
                title="Spatial Overview - Annotations",
                save_path="plots/spatial_overview.png",
                point_size=1,
                categorical=True,
            )
        logger.info("Generated spatial overview plot")
    except Exception as e:
        logger.error(f"Error generating spatial overview: {e}")
    
    logger.info("=== ResolVI Visualization Completed Successfully ===")
    
    # Write versions file
    def get_version(package):
        try:
            result = subprocess.run([sys.executable, "-c", f"import {package}; print({package}.__version__)"], 
                                  capture_output=True, text=True)
            return result.stdout.strip()
        except:
            return "unknown"
    
    versions_content = f'''"${task.process}":
        scanpy: {get_version('scanpy')}
        matplotlib: {get_version('matplotlib')}
        seaborn: {get_version('seaborn')}
        pandas: {get_version('pandas')}
        numpy: {get_version('numpy')}'''
    
    with open("versions.yml", "w") as f:
        f.write(versions_content)
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p plots
    touch plots/umap_annotation_comparison.png
    touch plots/umap_annotation_match.png
    touch plots/annotation_confusion_matrix.png
    touch plots/spatial_overview.png
    touch plots/spatial_annotation.png
    touch plots/spatial_resolvi_predicted.png
    
    cat << 'END_VERSIONS' > versions.yml
    "${task.process}":
        scanpy: 1.9.6
        matplotlib: 3.7.1
        seaborn: 0.12.0
        pandas: 2.0.3
        numpy: 1.24.3
    END_VERSIONS
    """
}

