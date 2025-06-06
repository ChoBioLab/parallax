process RESOLVI_VISUALIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::scanpy bioconda::matplotlib bioconda::seaborn"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.9.3--pyhd8ed1ab_0' :
        'quay.io/biocontainers/scanpy:1.9.3--pyhd8ed1ab_0' }"

    input:
    tuple val(meta), path(adata)
    val marker_genes

    output:
    tuple val(meta), path("plots/"), emit: plots
    path "versions.yml", emit: versions

    script:
    def genes_arg = marker_genes ? "--marker_genes '${marker_genes.join(',')}'" : ""
    """
    #!/usr/bin/env python3

    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import os
    import argparse
    import logging
    import warnings

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker_genes", type=str, default="", help="Comma-separated marker genes")
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Suppress warnings
    warnings.filterwarnings("ignore")
    plt.ioff()  # Turn off interactive mode

    # Create output directory
    os.makedirs("plots", exist_ok=True)

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad("${adata}")

    # Parse marker genes
    marker_genes = []
    if args.marker_genes:
        marker_genes = [gene.strip() for gene in args.marker_genes.split(',') if gene.strip()]

    # Filter marker genes to those present in dataset
    if marker_genes:
        available_markers = [gene for gene in marker_genes if gene in adata.var_names]
        unavailable_markers = [gene for gene in marker_genes if gene not in adata.var_names]
        if unavailable_markers:
            logger.warning(f"Genes not found in AnnData: {', '.join(unavailable_markers)}")
        logger.info(f"Found {len(available_markers)} of {len(marker_genes)} marker genes in dataset")
    else:
        available_markers = []

    # Custom spatial plotting function
    def xenium_spatial_plot(adata, color, layer=None, title=None, figsize=(8, 7), 
                           point_size=1, save_path=None, colormap="viridis", categorical=False):
        """Custom spatial plotting function for Xenium data using Matplotlib."""
        logger.info(f"Plotting spatial distribution for {color}")
        coords = adata.obsm["X_spatial"]

        if color in adata.obs_keys():
            values = adata.obs[color]
            categorical = True
        else:
            if layer is not None and layer in adata.layers:
                expr_matrix = adata.layers[layer]
            else:
                expr_matrix = adata.X

            if color not in adata.var_names:
                logger.error(f"Gene '{color}' not found in adata.var_names.")
                return None, None

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
            scatter = ax.scatter(coords[:, 0], coords[:, 1], s=point_size, c=values, cmap=colormap, alpha=0.7)
            plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.05)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title if title else color)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            logger.info(f"Saved plot to {save_path}")
        else:
            return fig, ax

    # Generate UMAP comparison visualizations
    logger.info("Generating UMAP comparison visualizations...")
    try:
        sc.pl.umap(adata, color=["annotation", "resolvi_predicted"], show=False, wspace=0.4)
        save_path = "plots/umap_annotation_comparison.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved UMAP annotation comparison to {save_path}")

        sc.pl.umap(adata, color="annotation_match", show=False)
        save_path = "plots/umap_annotation_match.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved UMAP annotation match to {save_path}")
    except Exception as e:
        logger.error(f"Error generating UMAP comparison plots: {e}")

    # Create confusion matrix
    logger.info("Generating confusion matrix for annotations...")
    try:
        confusion = pd.crosstab(
            adata.obs["annotation"],
            adata.obs["resolvi_predicted"],
            normalize="index",
        )
        plt.figure(figsize=(12, 10))
        sns.heatmap(confusion, cmap="viridis", annot=True, fmt=".2f", cbar_kws={"label": "Proportion"})
        plt.title("Original vs ResolVI Cell Type Assignment")
        plt.ylabel("Original Annotation")
        plt.xlabel("ResolVI Predicted")
        plt.tight_layout()
        save_path = "plots/annotation_confusion_matrix.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved annotation confusion matrix to {save_path}")
    except Exception as e:
        logger.error(f"Error generating confusion matrix: {e}")

    # Generate spatial plots for annotations and noise components
    logger.info("Generating spatial plots...")
    for category in ["annotation", "resolvi_predicted", "annotation_match"]:
        try:
            xenium_spatial_plot(
                adata,
                color=category,
                title=f"Spatial distribution of {category}",
                save_path=f"plots/spatial_{category}.png",
                point_size=1,
                categorical=True,
            )
        except Exception as e:
            logger.error(f"Error plotting category '{category}': {e}")

    # Generate spatial plots for noise components
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
            except Exception as e:
                logger.error(f"Error plotting component '{component}': {e}")

    # Generate gene expression comparison plots
    if available_markers:
        logger.info("Generating gene expression comparison plots...")
        num_genes_to_plot = min(20, len(available_markers))
        logger.info(f"Will plot for the first {num_genes_to_plot} available marker genes.")

        for gene in available_markers[:num_genes_to_plot]:
            try:
                fig, axs = plt.subplots(1, 2, figsize=(16, 7))
                coords = adata.obsm["X_spatial"]
                gene_loc = adata.var_names.get_loc(gene)

                # Original counts
                if hasattr(adata.layers["counts"], "toarray"):
                    values_orig = adata.layers["counts"][:, gene_loc].toarray().flatten()
                else:
                    values_orig = adata.layers["counts"][:, gene_loc]

                scatter1 = axs[0].scatter(coords[:, 0], coords[:, 1], s=1.0, c=values_orig, 
                                        cmap="viridis", alpha=0.7)
                axs[0].set_title(f"{gene} - Original (counts layer)")
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
                    axs[1].text(0.5, 0.5, "Corrected counts not available", 
                              ha="center", va="center", transform=axs[1].transAxes)

                axs[1].set_xticks([])
                axs[1].set_yticks([])

                plt.tight_layout()
                save_path = f"plots/gene_{gene}_comparison.png"
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                plt.close(fig)
                logger.info(f"Generated comparison plot for gene {gene}")
            except Exception as e:
                logger.error(f"Error plotting gene '{gene}': {e}")

    logger.info("Visualization completed successfully!")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
        f.write(f'    matplotlib: {plt.matplotlib.__version__}\\n')
        f.write(f'    seaborn: {sns.__version__}\\n')

    ${genes_arg}
    """
}
