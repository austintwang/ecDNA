import pickle
import os
import numpy as np
import pandas as pd
import seaborn as sns

def plot_facet_breaks(df, freq_col, seq_order, result_path):
    sns.set(style="whitegrid", font="Roboto")
    sns.relplot(data=df, x="seq_from", y="seq_to", hue=freq_col, seq_order=seq_order, hue_norm=(0,10), height=3, aspect=1, s=9)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def view_breaks(in_path, min_supports, results_dir):
    with open(in_path, "rb") as in_file:
        breaks_df, breaks, freqs_from, freqs_to = pickle.load(in_file)

    seq_order = sorted(list(set(breaks_df["seq_from"]) | set(breaks_df["seq_to"])))
    facet_dir = os.path.join(results_dir, "genome_pair")
    os.makedirs(facet_dir, exist_ok=True)

    for min_support in min_supports:
        df_from = breaks_df["freq_from"] >= min_support
        df_to = breaks_df["freq_to"] >= min_support
        df_pair = breaks_df["freq_pair"] >= min_support

        result_path_from = os.path.join(facet_dir, f"from_s_{min_support}.png")
        result_path_to = os.path.join(facet_dir, f"to_s_{min_support}.png")
        result_path_pair = os.path.join(facet_dir, f"pair_s_{min_support}.png") 

        plot_facet_breaks(df_from, "freq_from", seq_order, result_path_from)
        plot_facet_breaks(df_to, "freq_to", seq_order, result_path_to)
        plot_facet_breaks(df_pair, "freq_pair", seq_order, result_path_pair)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_breaks_proc.pickle")

    results_dir = "/oak/stanford/groups/wjg/atwang/ecdna/results"
    os.makedirs(results_dir, exist_ok=True)

    min_supports = [2]
    view_breaks(in_path, min_supports, results_dir)