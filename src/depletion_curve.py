import pandas as pd
import numpy as np
import streamlit as st
import altair as alt

@st.cache_data
def get_depletion_curve_data(GWMs_file, timepoints_file):
    timepoints = pd.read_csv(timepoints_file, index_col=0)
    Gs = timepoints.mean(axis=1).rename("G").round(3)

    GWMs = pd.read_csv(GWMs_file, index_col=0)
    LFCs = GWMs[["YES0", "YES1", "YES2", "YES3", "YES4"]].copy().stack().rename("LFC").round(3)
    pvalues = GWMs.filter(like="pvalue").copy()
    pvalues.columns = pvalues.columns.str.replace("_pvalue", "")
    pvalues = pvalues.stack().rename("pvalue")
    LFCs = pd.concat([LFCs, pvalues], axis=1).rename_axis(["Gene", "Timepoint"]).reset_index()
    LFCs = LFCs.merge(Gs, left_on="Timepoint", right_index=True, how="left")
    LFCs["Confidence"] = -np.log10(LFCs["pvalue"].apply(lambda x: 1e-10 if x <= 1e-10 else 1 if x > 1-1e-10 else x))

    return LFCs

def plot_depletion_curve(sub_LFCs):

    gene_selector = alt.selection_multi(fields=["Gene"], bind="legend")

    base = alt.Chart(sub_LFCs)

    line_plot = base.mark_line(size=4, opacity=0.2).encode(
        x=alt.X("G:Q", title="Generations"),
        y=alt.Y("LFC:Q", title="Log2 Fold Change", scale=alt.Scale(domain=(-3, 10))),
        color=alt.Color("Gene:N", title="Gene"),
        tooltip=["Gene", "G", "LFC", "pvalue"],
        opacity=alt.condition(gene_selector, alt.value(0.2), alt.value(0)),
    )
    point_plot = base.mark_circle(opacity=0.2).encode(
        x=alt.X("G:Q", title="Generations"),
        y=alt.Y("LFC:Q", title="Log2 Fold Change", scale=alt.Scale(domain=(-3, 10))),
        color=alt.Color("Gene:N", title="Gene", legend=alt.Legend(columns=10, symbolLimit=100, orient="bottom")),
        opacity=alt.condition(gene_selector, alt.value(0.2), alt.value(0)),
        size=alt.Size("Confidence:Q", title="-log10(pvalue)", legend=alt.Legend(
            orient="top"), scale=alt.Scale(type='sqrt')),
        tooltip=["Gene", "G", "LFC", "pvalue"]
    )

    line_point_plot = (line_plot+point_plot).properties(
        width=600,
        height=800
    ).add_params(gene_selector)

    return line_point_plot