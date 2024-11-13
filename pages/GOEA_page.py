import pandas as pd
from src.enrichment_analysis import load_GO_data, GOEA, format_GOEA_results, display_GOEA_results, parse_string_enrichment, display_string_enrichment
import streamlit as st
import altair as alt
import requests ## python -m pip install requests 
import json
from src.utils_functions import gene_information, get_gene_list_from_text_area, load_cluster_info
import io
from pathlib import Path

gene_info_file = Path(
    "./references/pombase_annotation/20241001/gene_IDs_names_products.tsv")
gene_id_to_name, coding_genes = gene_information(gene_info_file)
name_to_id = {v: k for k, v in gene_id_to_name.items()}

st.title("Enrichment analysis")
st.divider()
st.cache_data.clear()
# get query genes
input_container = st.container()
query_container, bg_container = input_container.columns(2)
query_container.header("Enter the list of the query genes")
query_genes, _ = get_gene_list_from_text_area(query_container, "query_genes", gene_id_to_name, name_to_id)

bg_container.header("Enter the list of the background genes")
bg_genes, _ = get_gene_list_from_text_area(bg_container, "bg_genes", gene_id_to_name, name_to_id)

# load GWMs
coding_genes_in_DIT_HAP, cluster_info_grouped = load_cluster_info("./data/0_raw/clustered_GWMs_customed_distance_1000_renamed.csv")
gene_sets = {"All coding genes": coding_genes, "Coding genes in DIT_HAP": coding_genes_in_DIT_HAP, **cluster_info_grouped}
gene_sets_with_size = {f"{k} ({len(v)})": v for k, v in gene_sets.items()}

st.sidebar.header("Configuration:")
st.sidebar.markdown("- **Background gene setting:**")
use_customed_bg = st.sidebar.toggle("Use customed background genes", value=False)
if use_customed_bg:
    bg_genes = bg_genes
else:
    bg_gene_set = st.sidebar.selectbox(
        "Please select the background gene set",
        gene_sets_with_size.keys(),
    )
    bg_genes = gene_sets_with_size[bg_gene_set]

st.sidebar.markdown("- **Query gene setting:**")
use_customed_query = st.sidebar.toggle("Use customed query genes", value=False)
if use_customed_query:
    query_genes = query_genes
else:
    query_gene_set = st.sidebar.selectbox(
        "Please select the query gene set",
        gene_sets_with_size.keys(),
    )
    query_genes = gene_sets_with_size[query_gene_set]

ontology_tab = st.tabs(["GO enrichment", "FYPO enrichment", "STRING enrichment"])
with ontology_tab[0]:
    with st.spinner("Performing ontology enrichment analysis..."):
        obo_file = Path(
                "./references/pombase_annotation/20241001/go-basic.obo")
        gaf_file = Path(
            "./references/pombase_annotation/20241001/go_style_gaf.tsv")
        godag, ns2assoc = load_GO_data(obo_file, gaf_file)
        go_sig = GOEA(query_genes, bg_genes, godag, ns2assoc)
        go_sig_results, reorder_columns = format_GOEA_results(go_sig, gene_id_to_name)
        go_sig_results = display_GOEA_results(ontology_tab[0], go_sig_results, reorder_columns)
with ontology_tab[1]:
    with st.spinner("Performing FYPO enrichment analysis..."):
        FYPO_file = Path(
                "./references/pombase_annotation/20241001/fypo-simple.obo")
        FYPO_gaf_file = Path(
            "./references/pombase_annotation/phaf_go_style_gaf.tsv")
        FYPOdag, ns2FYPOassoc = load_GO_data(FYPO_file, FYPO_gaf_file)
        fy_sig = GOEA(query_genes, bg_genes, FYPOdag, ns2FYPOassoc)
        fy_sig_results, reorder_columns = format_GOEA_results(fy_sig, gene_id_to_name)
        fy_sig_results = display_GOEA_results(ontology_tab[1], fy_sig_results, reorder_columns)
with ontology_tab[2]:
    with st.spinner("Performing STRING enrichment analysis..."):
        string_sig = parse_string_enrichment(query_genes, bg_genes)
        string_sig_results = display_string_enrichment(ontology_tab[2], string_sig)


# with ontology_tab[3]:
#     with st.spinner("Generating STRING network..."):
#         st.html("./references/string.html")

