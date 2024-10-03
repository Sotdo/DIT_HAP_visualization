from pathlib import Path
import streamlit as st

from src.load_basic_data import load_all_data

from src.link_with_known_information import display_basic_information
from src.extract_DIT_HAP_data import get_insertions_in_genes, get_DR_in_genes
from src.plot_insertion import plot_layer

# Set page configuration
st.set_page_config(layout="wide")

# Load data
data_folder = Path("./data")
references_folder = Path("./references")

Groups = {
    "Group1":  "HD_DIT_HAP",
}

(
    individual_insertion_data,
    weighted_data,
    weighted_data_stat,
    DIT_HAP_features,
    insertion_annotations,
    merged_gene_info,
    protein_domains_based_on_AF,
    genome_regions,
    DNA_region_for_domain,
) = load_all_data(data_folder, references_folder, Groups)

# Sidebar
st.sidebar.title("DIT HAP")

# Get query gene
st.sidebar.markdown("---")
st.sidebar.subheader("Enter the query gene")
query_gene = st.sidebar.text_input(
    "Enter the gene name or systematic ID", "SPAC1002.09c")

# Display basic information
basic_informations = st.container()
sysID = display_basic_information(
    basic_informations, query_gene, merged_gene_info)

# Extract data
DR_feature_in_current_gene = get_DR_in_genes(
    sysID, DIT_HAP_features, insertion_annotations)

GM_stats, merged_weighted_stats = get_insertions_in_genes(
    sysID,
    insertion_annotations,
    individual_insertion_data,
    weighted_data,
    weighted_data_stat,
    DIT_HAP_features,
)

vconcated_plot = plot_layer(
    sysID,
    DR_feature_in_current_gene,
    DNA_region_for_domain,
    genome_regions,
    GM_stats,
    individual_insertion_data,
    weighted_data,
    weighted_data_stat,
    Groups
)

st.altair_chart(vconcated_plot.resolve_scale(
    size="independent"), use_container_width=True, theme=None)
