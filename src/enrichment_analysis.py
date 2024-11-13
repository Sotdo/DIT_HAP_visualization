
import numpy as np
import pandas as pd
import altair as alt
import streamlit as st
import requests
import io
from requests.exceptions import ConnectionError, RequestException
import time

from goatools.obo_parser import GODag
from goatools.rpt.rpt_lev_depth import RptLevDepth
from goatools.anno.gaf_reader import GafReader
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.mapslim import mapslim
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
from goatools.godag.go_tasks import get_go2parents

from src.utils_functions import *

@st.cache_resource
def load_GO_data(obo_file, gaf_file):

    godag = GODag(str(obo_file), optional_attrs=["defn", "relationship"])
    objanno = GafReader(gaf_file, godag=godag)
    ns2assoc = objanno.get_ns2assc()

    return godag, ns2assoc

def GOEA(query_genes, bg_genes, godag, ns2assoc):
    goeaobj = GOEnrichmentStudyNS(
        bg_genes,  # List of mouse protein-coding genes
        ns2assoc,  # geneid/GO associations
        godag,  # Ontologies
        propagate_counts=False,
        alpha=0.05,  # default significance cut-off
        methods=["fdr_bh"],
    )  # defult multipletest correction method
    
    goea_results = goeaobj.run_study(query_genes, prt=None)
    goea_results_all = [
        r for r in goea_results if r.enrichment == "e"
    ]
    goea_results_sig = [
        r for r in goea_results if (r.p_fdr_bh < 0.05) and (r.enrichment == "e")
    ]

    return goea_results_sig


def format_GOEA_results(goea_results_sig, gene_id_to_name=None):
    enriched_variables = [
        'GO',
        'NS',
        'enrichment',
        'name',
        'p_fdr_bh',
        'p_uncorrected',
        'study_count',
        'pop_count',
        'study_n',
        'pop_n',
        'ratio_in_study',
        'ratio_in_pop',
        'study_items',
        'pop_items',
    ]
    enrichment_results = pd.DataFrame()
    for idx, enrichment in enumerate(goea_results_sig):
        for enriched_var in enriched_variables:
            var_value = enrichment.__getattribute__(enriched_var)
            if isinstance(var_value, set):
                if gene_id_to_name is not None:
                    var_value = [gene_id_to_name[gene_id] for gene_id in var_value]
                    var_value = ", ".join(var_value)
                else:
                    var_value = ", ".join(var_value)
            elif isinstance(var_value, tuple):
                var_value = "/".join([str(v) for v in var_value])
            enrichment_results.loc[idx, enriched_var] = var_value
        coverage_frac = enrichment_results.loc[idx, "study_count"] / enrichment_results.loc[idx, "pop_count"]
        enrichment_results.loc[idx, "coverage_frac"] = round(coverage_frac, 3)
        missing_items = enrichment.__getattribute__("pop_items") - enrichment.__getattribute__("study_items")
        missing_items = [gene_id_to_name[gene_id] for gene_id in missing_items]
        enrichment_results.loc[idx, "missing_items"] = ", ".join(missing_items)
    reorder_columns = enriched_variables[:-2] + ["coverage_frac"] + ["study_items", "missing_items", "pop_items"]
    return enrichment_results, reorder_columns

def display_GOEA_results(container, enrichment_results, reorder_columns):
    if enrichment_results.empty:
        message = "No significant terms found"
        container.warning(message)
    else:
        enrichment_results = enrichment_results.sort_values(by=["p_fdr_bh"], ascending=True)
        enrichment_results.reset_index(drop=True, inplace=True)
        enrichment_results = enrichment_results[reorder_columns]
        enrichment_results.rename(columns={"GO": "Term ID"}, inplace=True)
        message = "The following terms were found to be significantly enriched:"
        container.success(message)
        with st.expander("Enrichment result table", expanded=True):
            st.dataframe(enrichment_results, hide_index=True, use_container_width=True)
        sig_chart = plot_GOEA_results(enrichment_results)
        with st.expander("Enrichment result plots", expanded=True):
            st.altair_chart(sig_chart, use_container_width=True)
    
    return enrichment_results

def plot_GOEA_results(enrichment_results):
    charts = []
    for ns, ns_sig in enrichment_results.groupby("NS"):
        chart = (
            alt.Chart(ns_sig)
            .mark_circle()
            .encode(
                alt.X("study_count:Q", axis=alt.Axis(grid=True)),
                alt.Y(
                    "name:N",
                    axis=alt.Axis(grid=True, labelLimit=500, title=""),
                    sort=alt.EncodingSortField(
                        field="study_count", order="descending"),
                ),
                size="coverage_frac:Q",
                color=alt.Color("p_fdr_bh:Q").scale(
                    scheme="yelloworangered", reverse=True),
                tooltip=ns_sig.columns.tolist()
            )
        ).properties(
            title=f"Enrichment results for {ns}",
            width=120
        )
        charts.append(chart)

    return alt.hconcat(*charts)

def parse_string_enrichment(query_genes, bg_genes, max_retries=3, retry_delay=5):

    output_format = "tsv"
    # Try to get the current STRING version
    for attempt in range(max_retries):
        try:
            string_version = requests.post(f"https://string-db.org/api/{output_format}/version").text
            string_api_url = "https" + str(string_version.split("https")[1].strip()) + "/api"
            break
        except (ConnectionError, RequestException) as e:
            if attempt == max_retries - 1:
                st.error(f"Failed to connect to STRING database after {max_retries} attempts. Error: {str(e)}")
                return pd.DataFrame()  # Return an empty DataFrame
            time.sleep(retry_delay)

    method = "get_string_ids"
    params = {
        "identifiers" : "\r".join(bg_genes),
        "species" : 4896,
        "limit" : 1,
        "echo_query" : 1,
        "caller_identity" : "DIT_HAP_visualization"
    }
    
    # Try to get STRING IDs
    for attempt in range(max_retries):
        try:
            request_url = "/".join([string_api_url, output_format, method])
            results = requests.post(request_url, data=params)
            results.raise_for_status()  # Raise an exception for bad status codes
            break
        except (ConnectionError, RequestException) as e:
            if attempt == max_retries - 1:
                st.error(f"Failed to get STRING IDs after {max_retries} attempts. Error: {str(e)}")
                return pd.DataFrame()  # Return an empty DataFrame
            time.sleep(retry_delay)

    bg_string_ids = []
    for line in results.text.strip().split("\n")[1:]:
        l = line.split("\t")
        string_identifier = l[2]
        bg_string_ids.append(string_identifier)

    method = "enrichment"
    params = {
        "identifiers" : "%0d".join(query_genes),
        "background_string_identifiers" : "%0d".join(bg_string_ids),
        "species" : 4896,
        "caller_identity" : "DIT_HAP_visualization",
    }
    
    # Try to get enrichment results
    for attempt in range(max_retries):
        try:
            request_url = "/".join([string_api_url, output_format, method])
            response = requests.post(request_url, data=params)
            response.raise_for_status()  # Raise an exception for bad status codes
            break
        except (ConnectionError, RequestException) as e:
            if attempt == max_retries - 1:
                st.error(f"Failed to get enrichment results after {max_retries} attempts. Error: {str(e)}")
                return pd.DataFrame()  # Return an empty DataFrame
            time.sleep(retry_delay)

    data = response.text
    dataframe = pd.read_csv(io.StringIO(data), sep="\t")
    return dataframe

def display_string_enrichment(container, dataframe):
    if dataframe.empty:
        st.warning("No significant STRING enrichment found")
    else:
        kept_column_names = ["category", "term", "description", "p_value", "fdr", "number_of_genes", "number_of_genes_in_background", "inputGenes", "preferredNames"]
        dataframe = dataframe.loc[:, kept_column_names].copy()

        category_description = {
            "Process": "Biological Process (Gene Ontology)",
            "Component": "Cellular Component (Gene Ontology)",
            "Function": "Molecular Function (Gene Ontology)",
            "PMID": "Reference Publications (PubMed)",
            "NetworkNeighborAL": "Local Network Cluster (STRING)",
            "KEGG": "KEGG Pathways",
            "RCTM": "Reactome Pathways",
            "COMPARTMENTS": "Subcellular Localization (COMPARTMENTS)",
            "Keyword": "Annotated Keywords (UniProt)",
            "InterPro": "Protein Domains and Features (InterPro)",
            "SMART": "Protein Domains and Features (SMART)"
        }

        category_order = list(category_description.keys())
        dataframe['category_order'] = dataframe['category'].map({v: i for i, v in enumerate(category_order)})
        dataframe.sort_values(by='category_order', inplace=True)
        dataframe.drop('category_order', axis=1, inplace=True)
        dataframe["category"] = dataframe["category"].map(category_description)

        for category, category_df in dataframe.groupby("category", sort=False):
            container.subheader(category)
            container.dataframe(category_df, hide_index=True, use_container_width=True)

    return dataframe
