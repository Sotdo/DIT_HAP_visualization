import pandas as pd
import numpy as np
import streamlit as st
import altair as alt


def insertion_data_plot(insertion_GMs, gene_level_GMs, point_selector):

    insertion_GMs["Confidence"] = insertion_GMs["weights"].apply(lambda x: 0.1 if x < 0.1 else x)
    gene_level_GMs["Confidence"] = gene_level_GMs["weights"].apply(lambda x: 0.1 if x < 0.1 else x)

    insertion_GMs_plot = alt.Chart(insertion_GMs.reset_index()).mark_point(opacity=0.7).encode(
        x=alt.X('G:Q', title="Generations"),
        y=alt.Y('M:Q', scale=alt.Scale(domain=(-3, 9)), title="M values"),
        size=alt.Size("Confidence:Q", title="-log10(Padj)", legend=alt.Legend(
            orient="top"), scale=alt.Scale(type='sqrt')),
        order="Timepoint:N",
        tooltip=insertion_GMs.reset_index().columns.tolist()
    ).transform_filter(point_selector).add_params(point_selector)

    gene_level_GMs_plot = alt.Chart(gene_level_GMs.reset_index()).mark_circle(opacity=0.7).encode(
        x=alt.X('G:Q', title="Generations"),
        y=alt.Y('M:Q', scale=alt.Scale(domain=(-3, 9)), title="M values"),
        size=alt.Size("Confidence:Q", title="-log10(Padj)", legend=alt.Legend(
            orient="top")),
        tooltip=gene_level_GMs.reset_index().columns.tolist()
    ) + alt.Chart(gene_level_GMs.reset_index()).mark_line(width=2, opacity=0.8, color="gray").encode(
        x=alt.X('G:Q'),
        y=alt.Y('M:Q'),
    )

    return insertion_GMs_plot+gene_level_GMs_plot


def gene_feature_across_gene_plot(insertion_last_tp, point_selector):

    insertion_last_tp["Confidence"] = insertion_last_tp["weights"].apply(lambda x: 0.1 if x < 0.1 else x)

    gene_feature_plot = alt.Chart(insertion_last_tp.reset_index()).mark_circle(opacity=0.8).encode(
        x=alt.X("Fraction_to_start_codon:Q", scale=alt.Scale(domain=(0, 1))),
        y=alt.Y("M:Q", title="LFC",
                scale=alt.Scale(domain=(-3, 10))),
        size=alt.Size("Confidence:Q", title="-log10(Padj)", scale=alt.Scale(type='sqrt')),
        color=alt.condition(point_selector, 'Insertion_direction:N', alt.value(
            'lightgray'), legend=alt.Legend(orient="top")),
        # color = alt.Color('Insertion_direction:N', legend=alt.Legend(orient="top"), scale=alt.Scale(scheme='set1')),
        tooltip=insertion_last_tp.reset_index().columns.tolist()
    ).add_params(point_selector)  # .transform_filter(point_selector)
    return gene_feature_plot

def combine_plots(insertion_GMs, gene_level_GMs, insertion_last_tp):
    point_selector = alt.selection_multi(fields=["#Chr", "Coordinate", "Strand", "Target"])
    insertion_curve_plot = insertion_data_plot(insertion_GMs, gene_level_GMs, point_selector)
    YES4_plot = gene_feature_across_gene_plot(insertion_last_tp, point_selector)

    return alt.hconcat(insertion_curve_plot, YES4_plot).resolve_scale(x='independent', y="shared", size="shared")

