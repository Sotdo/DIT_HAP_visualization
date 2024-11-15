# %%
from pathlib import Path
import sys
import altair as alt
import streamlit as st
sys.path.append("/data/c/yangyusheng_optimized/DIT_HAP_visualization/Streamlit_DEseq2")
from src.load_basic_data import load_data, load_additional_info
from src.extract_DIT_HAP_data import get_insertions_in_genes
from src.plot_insertion import combine_plots
from src.link_with_known_information import display_basic_information
from src.utils_functions import get_gene_list_from_text_area, gene_information
# Set page configuration

# Load data
insertion_LFCs_file = "./data/0_raw/insertions_LFC.csv"
gene_level_LFCs_file = "./data/0_raw/GWMs.csv"
insertion_annotations_file = "./data/0_raw/DIT_HAP_20241001.annotated.csv"
timepoint_file = "./data/0_raw/samples_timepoints.csv"

gene_description_file = "./references/pombase_annotation/20241001/gene_IDs_names_products.tsv"
gene_essentiality_file = "./references/Hayles_2013_OB_merged_categories.xlsx"
genome_region_file = "./references/Genome_regions_CDS_intron_IGR_annotated.bed"


(
    insertion_LFCs,
    gene_level_LFCs,
    insertion_annotations,
    timepoints  
) = load_data(insertion_LFCs_file, gene_level_LFCs_file, insertion_annotations_file, timepoint_file)

(
    merged_gene_info,
    genome_regions,
) = load_additional_info(gene_description_file, gene_essentiality_file, genome_region_file)


# %%

# Sidebar
st.title("DIT HAP curve plot")
st.divider()

gene_info_file = Path(
    "./references/pombase_annotation/20241001/gene_IDs_names_products.tsv")
gene_id_to_name, coding_genes = gene_information(gene_info_file)
name_to_id = {v: k for k, v in gene_id_to_name.items()}

# Get query gene
st.sidebar.subheader("Enter the list of the query genes")
query_sysIDs, _ = get_gene_list_from_text_area(st.sidebar, "input genes", gene_id_to_name, name_to_id)

for query in query_sysIDs:
    
    igene = st.container()
    sysID, gene_col = display_basic_information(igene, query, merged_gene_info)
    if sysID in gene_level_LFCs.index:
        insertion_GMs, gene_level_GMs, insertion_last_tp = get_insertions_in_genes(sysID, insertion_annotations, insertion_LFCs, gene_level_LFCs, timepoints)
        combined_plot = combine_plots(insertion_GMs, gene_level_GMs, insertion_last_tp)
        gene_col.altair_chart(combined_plot, use_container_width=True, theme=None)
    else:
        st.warning(f"No data found for {query}")
    igene.divider()