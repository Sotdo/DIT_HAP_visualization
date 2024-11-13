from pathlib import Path

import pandas as pd

import streamlit as st
def display_basic_information(container, query, merged_gene_info):

    sysIDs_names = merged_gene_info[["gene_name"]].copy()
    names_sysIDs = merged_gene_info[["gene_name"]].copy().reset_index(drop=False).set_index("gene_name")
    synonyms = merged_gene_info["synonyms"].str.split(",").explode().str.strip().to_frame("synonyms").reset_index(drop=False)
    synonyms_sysIDs = synonyms.set_index("synonyms")
    if query in sysIDs_names.index:
        sysID = query
        gene_name = sysIDs_names.loc[query, "gene_name"]
    elif query in names_sysIDs.index:
        sysID = names_sysIDs.loc[query, "Systematic ID"]
        gene_name = query
    elif query in synonyms_sysIDs.index:
        sysID = synonyms_sysIDs.loc[query, "Systematic ID"]
        gene_name = sysIDs_names.loc[sysID, "gene_name"]
        container.text(f"The gene name {query} has been updated to {sysID} / {gene_name}")
    else:
        container.text(f"No gene information found for {query}")
    
    pombase_link = f"https://www.pombase.org/gene/{sysID}"
    # info_col, pombase_col = container.columns([6,4])
    info_col = st.container()
    if sysID == gene_name:
        info_col.header(f"Gene: [{sysID}]({pombase_link})")
    else:
        info_col.header(f"Gene: [{sysID} / {gene_name}]({pombase_link})")
    info_dict = {
        "Product": merged_gene_info.loc[sysID, 'gene_product'],
        "Essentiality (Hayles 2013)": merged_gene_info.loc[sysID, 'Gene dispensability. This study'],
        "Deletion phenotype": merged_gene_info.loc[sysID, 'Deletion mutant phenotype description'],
        "Category": merged_gene_info.loc[sysID, 'Category']
    }

    info_name_col, info_value_col = info_col.columns([2,8])
    for info_name, info_value in info_dict.items():
        info_name_col.markdown(f"**{info_name}:**")
        info_value_col.markdown(info_value)

    # pombase_col.write(
    #     f'<iframe width="100%" height="800" src="{pombase_link}"></iframe>',
    #     unsafe_allow_html=True,
    # )

    return sysID, info_col