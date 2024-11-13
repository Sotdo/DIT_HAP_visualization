import numpy as np
import pandas as pd
import streamlit as st


def get_insertions_in_genes(sysID, insertion_annotations, insertion_LFCs, gene_level_LFCs, timepoints):

    insertions_in_current_genes = insertion_annotations.query(
        "(`Systematic ID` == @sysID) and (Distance_to_stop_codon > 4)"
    ).index.intersection(insertion_LFCs.index)

    tps_without_init = timepoints.index.tolist()[1:]
    last_tp = tps_without_init[-1]

    insertion_Ms = insertion_LFCs.loc[insertions_in_current_genes]["log2FoldChange"].stack().rename("M")
    insertion_pvalues = insertion_LFCs.loc[insertions_in_current_genes]["padj"].stack().rename("Padj")
    insertion_GMs = pd.concat([insertion_Ms, insertion_pvalues], axis=1).rename_axis(["#Chr", "Coordinate", "Strand", "Target", "Timepoint"], axis=0)
    insertion_GMs["G"] = timepoints.loc[insertion_GMs.index.get_level_values(-1)].values
    insertion_GMs["weights"] = -np.log10(insertion_GMs["Padj"].apply(lambda x: 1e-10 if x <= 1e-10 else 1 if x > 1-1e-10 else x))

    insertion_info_cols = [
        "Type", "Distance_to_start_codon", "Distance_to_stop_codon",
        "Fraction_to_start_codon", "Fraction_to_stop_codon", "Residue_affected", 
        "Residue_frame", "Insertion_direction"
    ]

    insertion_GMs[insertion_info_cols] = insertion_annotations.loc[insertion_GMs.index.droplevel(-1), insertion_info_cols].values

    gene_level_Ms = gene_level_LFCs.loc[sysID, tps_without_init].rename("M")
    gene_level_pvalues = gene_level_LFCs.filter(like="_pvalue")
    gene_level_pvalues.columns = [f"{col.split('_pvalue')[0]}" for col in gene_level_pvalues.columns]
    gene_level_pvalues = gene_level_pvalues.loc[sysID, tps_without_init].rename("Padj")
    gene_level_GMs = pd.concat([gene_level_Ms, gene_level_pvalues], axis=1)
    gene_level_GMs["G"] = timepoints.loc[gene_level_GMs.index.get_level_values(-1)].values
    gene_level_GMs["weights"] = -np.log10(gene_level_GMs["Padj"].apply(lambda x: 1e-10 if x <= 1e-10 else 1 if x > 1-1e-10 else x))

    insertion_last_tp = insertion_GMs.query("Timepoint == @last_tp")

    return insertion_GMs, gene_level_GMs, insertion_last_tp
