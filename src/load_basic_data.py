import pandas as pd
import streamlit as st

@st.cache_data
def load_data(insertion_LFCs_file: str, gene_level_LFCs_file: str, insertion_annotations_file: str, timepoint_file: str) -> tuple:
    """
    Load and process data from CSV files.

    Args:
        insertion_LFCs_file (str): Path to the insertion LFCs CSV file.
        gene_level_LFCs_file (str): Path to the gene-level LFCs CSV file.
        insertion_annotations_file (str): Path to the insertion annotations CSV file.
        timepoint_file (str): Path to the timepoint CSV file.
    Returns:
        tuple: A tuple containing three pandas DataFrames:
            - insertion_LFCs: DataFrame with insertion LFCs data.
            - gene_level_LFCs: DataFrame with gene-level LFCs data.
            - insertion_annotations: DataFrame with insertion annotations data.
            - timepoints: DataFrame with timepoints data.
    """

    insertion_LFCs = pd.read_csv(insertion_LFCs_file, index_col=[0, 1, 2, 3], header=[0, 1]).reorder_levels([1, 0], axis=1)
    gene_level_LFCs = pd.read_csv(gene_level_LFCs_file, index_col=0, header=0)
    insertion_annotations = pd.read_csv(insertion_annotations_file, index_col=[0, 1, 2, 3], header=0)
    timepoints = pd.read_csv(timepoint_file, index_col=0, header=0).mean(axis=1).sort_values()

    return insertion_LFCs, gene_level_LFCs, insertion_annotations, timepoints


@st.cache_data
def load_additional_info(gene_description_file: str, gene_essentiality_file: str, genome_region_file: str) -> tuple:
    """
    Load and process additional information from various files.

    Args:
        gene_description_file (str): Path to the gene description TSV file.
        gene_essentiality_file (str): Path to the gene essentiality Excel file.
        genome_region_file (str): Path to the genome region TSV file.

    Returns:
        tuple: A tuple containing two pandas DataFrames:
            - merged_gene_info: DataFrame with merged gene information.
            - genome_regions: DataFrame with genome region information.
    """
    # Load gene descriptions
    gene_info = pd.read_csv(gene_description_file, sep="\t")[
        ["gene_systematic_id", "gene_name", "gene_product", "synonyms"]]
    
    # Load gene essentiality data
    essentiality_columns = [
        "Systematic ID", "Deletion mutant phenotype description", 
        "Phenotypic classification used for analysis", "Gene dispensability. This study", 
        "Category", "One or multi basic phenotypes"
    ]
    essentiality_from_hayles = pd.read_excel(gene_essentiality_file, sheet_name="All genes")[essentiality_columns]

    # Fill NaN gene names with systematic IDs
    gene_info["gene_name"] = gene_info["gene_name"].fillna(gene_info["gene_systematic_id"])

    # Merge gene info with essentiality data
    merged_gene_info = gene_info.merge(
        essentiality_from_hayles,
        how="left",
        left_on="gene_systematic_id",
        right_on="Systematic ID"
    ).drop(columns=["Systematic ID"]).set_index("gene_systematic_id").rename_axis("Systematic ID")

    # Load genome regions
    genome_regions = pd.read_csv(genome_region_file, sep="\t", header=0)

    return merged_gene_info, genome_regions

