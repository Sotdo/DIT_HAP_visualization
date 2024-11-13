import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import altair as alt
from src.utils_functions import get_gene_list_from_text_area, gene_information, load_cluster_info
from src.depletion_curve import get_depletion_curve_data, plot_depletion_curve
# Sidebar
st.title("The similarity analysis of the depletion curve")
st.divider()
st.cache_data.clear()

gene_info_file = Path(
    "./references/pombase_annotation/20241001/gene_IDs_names_products.tsv")
gene_id_to_name, coding_genes = gene_information(gene_info_file)
name_to_id = {v: k for k, v in gene_id_to_name.items()}

LFCs = get_depletion_curve_data("./data/0_raw/GWMs.csv", "./data/0_raw/samples_timepoints.csv")
coding_genes_in_DIT_HAP, cluster_info_grouped = load_cluster_info("./data/0_raw/clustered_GWMs_customed_distance_1000_renamed.csv")
gene_sets = {"Coding genes in DIT_HAP": coding_genes_in_DIT_HAP, **cluster_info_grouped}
gene_sets_with_size = {f"{k} ({len(v)})": v for k, v in gene_sets.items()}

use_customed_gene_list = st.sidebar.toggle("Use customed gene list", value=False)
if use_customed_gene_list:
    st.sidebar.subheader("Enter the list of the query genes")
    query_sysIDs, _ = get_gene_list_from_text_area(st.sidebar, "input genes", gene_id_to_name, name_to_id)
else:
    gene_set = st.sidebar.selectbox(
        "Please select the gene set",
        gene_sets_with_size.keys(),
        index=11
    )
    query_sysIDs = gene_sets_with_size[gene_set]
with st.spinner("Plot the depletion curve..."):
    query_LFCs = LFCs.query("Gene in @query_sysIDs").copy()
    query_LFCs["Gene"] = query_LFCs["Gene"].map(gene_id_to_name)
    line_point_plot = plot_depletion_curve(query_LFCs)
    st.altair_chart(line_point_plot, use_container_width=False, theme=None)




