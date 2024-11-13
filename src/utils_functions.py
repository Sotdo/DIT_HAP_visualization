import pandas as pd
import streamlit as st

@st.cache_data
def gene_information(gene_info_file):

    gene_info = pd.read_csv(gene_info_file, index_col=0, sep="\t")
    # fill the NA in gene_name with index
    gene_info["gene_name"] = gene_info["gene_name"].fillna(gene_info.index.to_series())
    gene_id_to_name = gene_info["gene_name"].to_dict()
    coding_genes = gene_info.query("gene_type == 'protein coding gene'").index.tolist()
    
    return gene_id_to_name, coding_genes

def transform_query_genes_to_sysIDs(query_genes, gene_id_to_name, name_to_id):

    query_sysIDs = []
    missing_genes = []
    for query in query_genes:
        if query in gene_id_to_name:
            query_sysIDs.append(query)
        elif query in name_to_id:
            query_sysIDs.append(name_to_id[query])
        else:
            missing_genes.append(query)
    print_missing_genes = "\n".join(missing_genes)
    print_info = f"There are {len(query_genes)} genes in the list.\n{len(query_sysIDs)} were found.\nThe following {len(missing_genes)} genes were not found:\n{print_missing_genes}"

    return query_sysIDs, missing_genes, print_info

def get_gene_list_from_text_area(container, key, gene_id_to_name, name_to_id):
    genes = container.text_area(
        "Enter the genes (each gene by one row)", "SPAC1002.09c\nSPAC3G9.12", key=key)
    if "," in genes:
        genes = genes.replace(",", "\n")
    genes = genes.strip().split("\n")
    genes = [gene.strip() for gene in genes]
    query_sysIDs, missing_genes, print_info = transform_query_genes_to_sysIDs(genes, gene_id_to_name, name_to_id)
    container.text(print_info)
    return query_sysIDs, genes

@st.cache_data
def load_cluster_info(cluster_info_file):
    cluster_info = pd.read_csv(cluster_info_file)
    coding_genes_in_DIT_HAP = cluster_info["Systematic ID"].tolist()
    cluster_info_grouped = cluster_info.groupby("revised_cluster")["Systematic ID"].apply(list).to_dict()
    return coding_genes_in_DIT_HAP, cluster_info_grouped
    