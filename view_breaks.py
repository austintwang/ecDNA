import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_facet_breaks(df, freq_col, seq_order, result_path):
    print(df) ####
    sns.set(style="whitegrid", font="Roboto")
    sns.relplot(
        data=df, 
        x="pos_from",
        y="pos_to",
        row="seq_to", 
        col="seq_from",
        hue=freq_col, 
        row_order=reversed(seq_order), 
        col_order=seq_order, 
        hue_norm=(0,5), 
        height=5, 
        aspect=1, 
        s=30
    )
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def view_breaks(in_path, min_supports, seq_sets, results_dir):
    with open(in_path, "rb") as in_file:
        breaks_df, breaks, freqs_from, freqs_to = pickle.load(in_file)

    # seq_order = sorted(list(set(breaks_df["seq_from"]) | set(breaks_df["seq_to"])))
    # print(seq_order) ####
    facet_dir = os.path.join(results_dir, "genome_pair")
    os.makedirs(facet_dir, exist_ok=True)

    for k, v in seq_sets.items():
        for min_support in min_supports:
            df_from = breaks_df[breaks_df["freq_from"] >= min_support]
            df_to = breaks_df[breaks_df["freq_to"] >= min_support]
            df_pair = breaks_df[breaks_df["freq_pair"] >= min_support]

            result_path_from = os.path.join(facet_dir, f"{k}_from_s_{min_support}.svg")
            result_path_to = os.path.join(facet_dir, f"{k}_to_s_{min_support}.svg")
            result_path_pair = os.path.join(facet_dir, f"{k}_pair_s_{min_support}.svg") 

            plot_facet_breaks(df_from, "freq_from", v, result_path_from)
            plot_facet_breaks(df_to, "freq_to", v, result_path_to)
            plot_facet_breaks(df_pair, "freq_pair", v, result_path_pair)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_breaks_proc.pickle")

    results_dir = "/oak/stanford/groups/wjg/atwang/ecdna/results"
    os.makedirs(results_dir, exist_ok=True)

    min_supports = [1, 2]
    seq_sets = {
        # "3_8": ["NC_000003.12, +", "NC_000003.12, -", "NC_000008.11, +", "NC_000008.11, -"],
        "3_22P": ["NC_000003.12, +", "NC_000003.12, -", "NW_021160026.1, +", "NW_021160026.1, -"]
    }
    view_breaks(in_path, min_supports, seq_sets, results_dir)